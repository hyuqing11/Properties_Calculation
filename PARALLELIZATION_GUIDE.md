# Parallelization Strategy for HPC

## Quick Decision Guide

| Your Scenario | Best Choice | Expected Speedup |
|---------------|-------------|------------------|
| Single node (2-64 cores) | **Multiprocessing** | 2-8x |
| Multiple nodes, small jobs | **Multiprocessing** | 2-8x |
| Multiple nodes, large jobs | **MPI** | 10-100x |
| Many dump files to process | **MPI** or **Job arrays** | Nx files |

---

## Strategy 1: Multiprocessing (Recommended for most users)

### What to parallelize:
- Loop over atom types in `Simulator._compute_vacf_and_pdos()`
- Loop over atom types in `Simulator._compute_dynamic_structure()`
- Loop over q-vectors in `ComputeDynamicProperties.calculate_intermediate_scattering()`

### Implementation approach:

```python
from multiprocessing import Pool
import functools

def compute_pdos_for_atom_type(j, atom_pos_dict, atom_vel_dict, pos, vel, parameters, latt, omega):
    """Wrapper function for parallel PDOS computation."""
    if j == parameters['num_types'] + 1:
        pdos_cal = ComputeDynamicProperties(pos, vel, parameters, latt)
    else:
        pdos_cal = ComputeDynamicProperties(atom_pos_dict[j], atom_vel_dict[j], parameters, latt)

    vacf_non, vacf_output, pdos = pdos_cal.pdos(omega)
    return j, vacf_non, vacf_output, pdos

# In Simulator._compute_vacf_and_pdos():
omega = np.arange(0, self.parameters['max_omega'], self.parameters['d_omega'])
nu = omega / (2 * np.pi)

# Create partial function with fixed arguments
compute_func = functools.partial(
    compute_pdos_for_atom_type,
    atom_pos_dict=self.atom_pos_dict,
    atom_vel_dict=self.atom_vel_dict,
    pos=self.pos,
    vel=self.vel,
    parameters=self.parameters,
    latt=self.latt,
    omega=omega
)

# Parallel computation
n_processes = min(len(self.parameters['compute_type']), os.cpu_count())
with Pool(processes=n_processes) as pool:
    results = pool.map(compute_func, self.parameters['compute_type'])

# Write results
for j, vacf_non, vacf_output, pdos in results:
    self._write_vacf_and_pdos(j, nu, pdos, t, vacf_non, vacf_output, folder)
```

### HPC batch script:

```bash
#!/bin/bash
#SBATCH --job-name=properties_calc
#SBATCH --nodes=1                    # Single node
#SBATCH --ntasks=1                   # Single task
#SBATCH --cpus-per-task=16          # 16 cores for multiprocessing
#SBATCH --mem=64G
#SBATCH --time=02:00:00

module load python/3.11
export OMP_NUM_THREADS=1            # Disable NumPy threading
export MKL_NUM_THREADS=1

python main.py input.json
```

---

## Strategy 2: MPI (For large-scale HPC)

### What to parallelize:
- Distribute atom types across MPI ranks
- Distribute q-vectors across MPI ranks
- Parallelize time-correlation calculations with MPI reductions

### Implementation approach:

```python
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Distribute atom types across ranks
compute_types = self.parameters['compute_type']
my_types = compute_types[rank::size]  # Round-robin distribution

for j in my_types:
    if j == self.parameters['num_types'] + 1:
        pdos_cal = ComputeDynamicProperties(self.pos, self.vel, self.parameters, self.latt)
    else:
        pdos_cal = ComputeDynamicProperties(self.atom_pos_dict[j], self.atom_vel_dict[j],
                                           self.parameters, self.latt)
    vacf_non, vacf_output, pdos = pdos_cal.pdos(omega)

    # Each rank writes its own results
    self._write_vacf_and_pdos(j, nu, pdos, t, vacf_non, vacf_output, folder)

# Synchronize before exit
comm.Barrier()
```

### For even better scaling - parallelize correlation loops:

```python
# In ComputeDynamicProperties.pdos()
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

M = self.parameters['num_frame'] - self.parameters['Nc']
vel_array = np.array(self.atom_velocity)

# Each rank computes a subset of correlation steps
my_nc_values = range(rank, self.parameters['Nc'], size)
local_vacf = np.zeros(self.parameters['Nc'])

for nc in my_nc_values:
    correlations = np.sum(vel_array[0:M+1] * vel_array[nc:M+1+nc], axis=(1, 2))
    local_vacf[nc] = np.sum(correlations)

# Reduce across all ranks
vacf = np.zeros(self.parameters['Nc'])
comm.Allreduce(local_vacf, vacf, op=MPI.SUM)
```

### HPC batch script:

```bash
#!/bin/bash
#SBATCH --job-name=properties_mpi
#SBATCH --nodes=4                    # Multiple nodes
#SBATCH --ntasks=128                 # 128 MPI ranks
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00

module load python/3.11
module load openmpi/4.1.2

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

mpirun -np 128 python main.py input.json
```

---

## Strategy 3: Hybrid MPI + Threading (Advanced)

For maximum performance on modern HPC clusters with many-core nodes:

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8         # 8 MPI ranks per node
#SBATCH --cpus-per-task=8           # 8 threads per rank
#SBATCH --mem=128G

module load python/3.11
module load openmpi/4.1.2

export OMP_NUM_THREADS=8            # NumPy will use 8 threads
export MKL_NUM_THREADS=8

mpirun -np 32 python main.py input.json
```

- MPI parallelizes over atom types and q-vectors
- NumPy/BLAS threads parallelize within vector operations
- Best for 64-256 core jobs

---

## Performance Expectations

Based on typical molecular dynamics workloads:

### File Reading (already optimized):
- **Current:** ~3-5x faster than original
- **Multiprocessing:** No benefit (I/O bound)
- **MPI:** No benefit (I/O bound)

### PDOS Calculation:
- **Current vectorization:** ~10-50x faster than original
- **+ Multiprocessing (8 cores):** 3-6x additional speedup
- **+ MPI (128 ranks):** 10-80x additional speedup
- **Total potential:** 100-4000x faster than original code

### Intermediate Scattering:
- **Current vectorization:** ~50-100x faster than original
- **+ Multiprocessing (8 cores):** 4-7x additional speedup
- **+ MPI (128 ranks):** 20-100x additional speedup
- **Total potential:** 1000-10000x faster than original code

---

## When NOT to Parallelize

1. **Small systems** (<1000 atoms, <100 frames):
   - Overhead > benefit
   - Just use current optimizations

2. **Single atom type calculations**:
   - Nothing to parallelize at top level
   - Must parallelize inner loops (more complex)

3. **Limited HPC allocation**:
   - Current optimizations are sufficient
   - Save resources for other jobs

---

## Recommended Next Steps

1. **Start simple:** Add multiprocessing to atom type loop (1-2 hours work)
2. **Test scaling:** Run with 1, 2, 4, 8, 16 cores and measure speedup
3. **If good scaling:** Consider MPI for production runs
4. **Profile:** Use `cProfile` or HPC profiling tools to find remaining bottlenecks

---

## Quick Benchmark Command

```bash
# Test current optimizations
time python main.py input.json

# Test with multiprocessing (after implementation)
export NUM_PROCESSES=8
time python main.py input.json

# Test with MPI (after implementation)
mpirun -np 16 python main.py input.json
```

Expected results for 10k atoms, 1k frames, 4 atom types:
- **Serial (original):** 2-4 hours
- **Optimized (current):** 10-30 minutes
- **+ Multiprocessing (8 cores):** 3-8 minutes
- **+ MPI (64 ranks):** 30-90 seconds
