# Multiprocessing Usage Guide

This guide explains how to use the multiprocessing features added to the Properties Calculation package.

## Overview

Multiprocessing has been added to parallelize computations across different atom types. This provides significant speedup when computing properties for multiple atom types simultaneously.

## Affected Calculations

The following property calculations now support multiprocessing:

1. **VACF and PDOS** (property_type = 0)
2. **Van Hove Self Correlation** (property_type = 1)
3. **Intermediate Scattering Function and Dynamic Structure** (property_type = 2)

## Configuration

### Enabling Multiprocessing

Add the `num_processes` parameter to your `input.json` file:

```json
{
  "num_atoms": 10000,
  "num_frame": 1000,
  "num_types": 2,
  "num_processes": -1,
  ...
}
```

### `num_processes` Options:

| Value | Behavior |
|-------|----------|
| `-1` | Use all available CPU cores (default) |
| `1` | Disable multiprocessing (sequential execution) |
| `N > 1` | Use exactly N processes (capped at available cores) |

### Examples:

**Use all available cores:**
```json
"num_processes": -1
```

**Use 8 cores:**
```json
"num_processes": 8
```

**Disable multiprocessing (sequential):**
```json
"num_processes": 1
```

**If not specified:** Defaults to using all available cores

---

## Performance Guidelines

### When Multiprocessing Helps:

✅ **Multiple atom types**: Computing properties for 2+ atom types
✅ **Large systems**: 5000+ atoms, 500+ frames
✅ **Multi-core systems**: HPC nodes with 8+ cores
✅ **Computationally intensive**: PDOS, intermediate scattering

### When to Use Sequential:

❌ **Single atom type**: No parallelism benefit
❌ **Small systems**: <1000 atoms, <100 frames (overhead > benefit)
❌ **Limited resources**: Systems with few cores or memory constraints

---

## Expected Speedup

Based on typical molecular dynamics workloads:

| Configuration | Speedup |
|---------------|---------|
| 2 atom types, 4 cores | 1.8-2x |
| 2 atom types, 8 cores | 1.9-2x |
| 4 atom types, 8 cores | 3.5-4x |
| 4 atom types, 16 cores | 3.8-4x |

*Note: Speedup is limited by the number of atom types and available cores*

---

## HPC Usage

### SLURM Batch Script Example:

```bash
#!/bin/bash
#SBATCH --job-name=props_mp
#SBATCH --nodes=1                    # Single node
#SBATCH --ntasks=1                   # Single task
#SBATCH --cpus-per-task=16          # 16 cores for multiprocessing
#SBATCH --mem=64G
#SBATCH --time=02:00:00

module load python/3.11

# Important: Disable NumPy internal threading to avoid oversubscription
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Run with multiprocessing enabled
python main.py input.json
```

### Key Environment Variables:

Set these to prevent NumPy from using internal threading, which would conflict with multiprocessing:

```bash
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
```

---

## Example Input Files

### For Multiprocessing (Multiple Atom Types):

```json
{
  "property_type": 0,
  "num_atoms": 10000,
  "num_frame": 1000,
  "num_types": 4,
  "dim": 3,
  "dt": 0.001,
  "compute_type": [1, 2, 3, 4],
  "num_processes": -1,
  "max_omega": 20,
  "d_omega": 0.1,
  "Nc": 500
}
```

### For Sequential (Single Atom Type):

```json
{
  "property_type": 0,
  "num_atoms": 10000,
  "num_frame": 1000,
  "num_types": 2,
  "dim": 3,
  "dt": 0.001,
  "compute_type": [1],
  "num_processes": 1,
  "max_omega": 20,
  "d_omega": 0.1,
  "Nc": 500
}
```

---

## Benchmarking

To benchmark the performance improvement:

### Sequential Run:

```bash
# Set num_processes to 1 in input.json
time python main.py input.json
```

### Parallel Run (8 cores):

```bash
# Set num_processes to 8 in input.json
time python main.py input.json
```

### Parallel Run (all cores):

```bash
# Set num_processes to -1 in input.json
time python main.py input.json
```

Compare the execution times to measure speedup.

---

## Troubleshooting

### "Using sequential computation" message

**Causes:**
- `num_processes` set to 1
- Only one atom type in `compute_type`
- Multiprocessing is automatically disabled

**Solution:** Check your `input.json` configuration.

### Memory Issues

**Problem:** System runs out of memory with multiprocessing

**Solution:**
```json
"num_processes": 4  // Reduce number of processes
```

### No Speedup Observed

**Possible causes:**
1. System is too small (overhead dominates)
2. NumPy threading conflicts (set OMP_NUM_THREADS=1)
3. I/O bound (file reading/writing dominates)
4. Only one atom type being computed

---

## Technical Details

### Implementation:

- Uses Python's `multiprocessing.Pool` for parallel execution
- Each atom type computation runs in a separate process
- Results are collected and written sequentially (to avoid file conflicts)

### Memory Considerations:

Each process receives a copy of the atom position/velocity data. Memory usage:

```
Total Memory ≈ Base Memory + (num_processes × Per-Process Memory)
```

For large systems (>100k atoms), limit `num_processes` to avoid memory exhaustion.

### Thread Safety:

File writing is done sequentially after parallel computation to ensure thread safety and prevent corrupted output files.

---

## Comparison with MPI

| Feature | Multiprocessing | MPI |
|---------|----------------|-----|
| **Ease of use** | ✅ Simple | ⚠️ Complex |
| **Setup** | ✅ No dependencies | ⚠️ Requires mpi4py |
| **Best for** | 1-2 nodes | Multiple nodes |
| **Scalability** | Up to ~64 cores | Hundreds of cores |
| **Memory** | Shared on node | Distributed |

**Recommendation:**
- Start with multiprocessing for easy performance gains
- Consider MPI for large-scale HPC deployments (see PARALLELIZATION_GUIDE.md)

---

## Questions or Issues?

If you encounter problems or have questions about using multiprocessing:

1. Check this guide for configuration issues
2. Review PARALLELIZATION_GUIDE.md for advanced usage
3. Verify your input.json configuration
4. Check SLURM logs for error messages
