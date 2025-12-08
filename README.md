# Material Properties Calculation

This repository contains code for calculating various material properties using data obtained from molecular dynamics simulations (LAMMPS). The code has been optimized for large atom files and supports multiprocessing for parallel computation.

## 🚀 Performance Features

- **Optimized file I/O**: 3-5x faster reading of large LAMMPS dump files
- **Vectorized computations**: 10-100x speedup for PDOS and dynamic structure calculations
- **Multiprocessing support**: Parallel computation across atom types (see `MULTIPROCESSING_USAGE.md`)
- **HPC ready**: Batch scripts and configuration for SLURM clusters (see `PARALLELIZATION_GUIDE.md`)

## Supported Calculations

Currently, the following calculations are supported:

## 1. Plot Strings in OVITO (property_type = 6)
This option allows you to visualize the string configuration using OVITO software. Strings are accessed by connecting mobile atoms i and j if:
$$|\vec{r}_i(t)-\vec{r}_j(0)|<\delta$$
or
$$|\vec{r}_j(t)-\vec{r}_i(0)|<\delta$$

### Input Data Required (input.json):
- `num_atoms`: Total number of atoms in the system
- `num_frame`: Total number of steps
- `num_types`: Number of types of atoms in the system
- `dim`: The dimension you computed
- `dt`: The time step
- `gap`: `gap * dt` is the t in the equation above
- `steps_read`: The folder that you save for the STRINGS at time interval t
- `initial_read`: The folder that you save for the STRINGS at a specific initial time
- `time_index`: The time steps you would like to plot the strings

### Required Files:
- `STRINGS`: Includes information about the strings
- `conf.dump_all`: Includes atom positions and velocities

## 2. Van Hove Self Correlation (property_type = 1)
The Van Hove self-correlation function quantifies the probability of finding an ion within a volume element dr, centered at r and time t, given that the ion was at the origin at time t=0. The formula is given by:
$$G_s(\vec{r},t)=\left<\frac{1}{N} \sum_{i=1}^{N} \delta(\vec{r}-|\vec{r}_i(t)-\vec{r}_i(0)|)\right>$$

### Input Data Required (input.json):
- `num_atoms`: Total number of atoms in the system
- `num_frame`: Total number of steps
- `num_types`: Number of types of atoms in the system
- `dim`: The dimension you computed
- `dt`: The time step
- `rDel`: The value of dr
- `rCutoff`: The cut-off distance
- `time_series`: The array that includes the different time intervals you would like to compute \( G_s \)
- `gap`: For one run, the start of computing \( G_s \) is delayed by gap steps (different initial time)
- `ave_num`: The number of \( G_s \) calculations you are averaging
### Required Files:
- `conf.dump_all`: Includes atom positions and velocities

## 3. VACF and PDOS (property_type = 0)
This option calculates the Velocity AutoCorrelation Function (VACF) and the Phonon Density of States (PDOS) for the material.

- **VACF (Velocity AutoCorrelation Function):**
  VACF provides information about particle motion and diffusion. It is defined as:
  $$C(t) = \frac{1}{N} \sum_{i=1}^{N} \vec{v}_i(t) \cdot \vec{v}_i(0)$$
  Where $\vec{v}_i(t)$ is the velocity of particle i at time t, and N is the total number of particles.

- **PDOS (Phonon Density of States):**
  PDOS gives insights into the material's vibrational modes. It is calculated using Fourier transform techniques and provides the distribution of vibrational frequencies in the material.

### Input Data Required (input.json):
- `num_atoms`: Total number of atoms in the system
- `num_frame`: Total number of steps
- `num_types`: Number of types of atoms in the system
- `dim`: The dimension you computed
- `dt`: The time step
- `max_omega`: The largest frequency to compute PDOS
- `d_omega`: The frequency step
- `Nc`: The correlation time for VACF
- `compute_type`: The array that specifies which atom type you are going to compute the VACF and PDOS for. If `compute_type` equals `num_types + 1`, then compute VACF for all the atoms.

### Required Files:
- `conf.dump_all`: Includes atom positions and velocities

## 4. Intermediate Scattering Function and Dynamic Structure Factor (property_type=2)

This option calculates the **intermediate scattering function** F(q,t), the **dynamic structure factor** S(q,ω), and their **momentum-integrated** forms. These quantities characterize the dynamics of density fluctuations in materials and are directly comparable to neutron/X-ray scattering experiments.

### Physical Quantities:

- **Intermediate Scattering Function F(q,t):**
  Describes the time evolution of density fluctuations at wavevector q:
  $$F(\vec{q},t) = \frac{1}{N} \left\langle \sum_{i=1}^{N} \sum_{j=1}^{N} e^{i\vec{q} \cdot [\vec{r}_i(0) - \vec{r}_j(t)]} \right\rangle$$

  The real part is computed as:
  $$F(\vec{q},t) = \frac{1}{N} \left\langle \left[\sum_{i} \cos(\vec{q} \cdot \vec{r}_i(0))\right] \times \left[\sum_{j} \cos(\vec{q} \cdot \vec{r}_j(t))\right] + \left[\sum_{i} \sin(\vec{q} \cdot \vec{r}_i(0))\right] \times \left[\sum_{j} \sin(\vec{q} \cdot \vec{r}_j(t))\right] \right\rangle$$

- **Dynamic Structure Factor S(q,ω):**
  The Fourier transform of F(q,t) with respect to time:
  $$S(\vec{q},\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} F(\vec{q},t) e^{-i\omega t} dt$$

  This reveals the frequency spectrum of density fluctuations at each wavevector.

- **Integrated Dynamic Structure Factor:**
  Momentum-integrated S(q,ω) over specified q-ranges:
  $$S_{int}(\omega, [q_{min}, q_{max}]) = \sum_{q \in [q_{min},q_{max}]} S(q,\omega)$$

  Useful for comparison with angle-integrated scattering experiments.

### Input Data Required (input.json):

**Basic parameters:**
- `num_atoms`: Total number of atoms in the system
- `num_frame`: Total number of time steps
- `num_types`: Number of atom types in the system
- `dim`: Dimension (typically 3)
- `dt`: Time step (in ps)
- `Nc`: Number of correlation steps to compute F(q,t)
- `compute_type`: Array specifying which atom types to compute (e.g., `[1, 2]`)

**Q-vector parameters:**
- `vectors`: Number of q-vectors to compute (e.g., `50`)
- `q_dir`: Direction of q-vectors in reciprocal space (e.g., `[1, 0, 0]` for x-direction)
- `uCell`: Unit cell dimensions (e.g., `[4, 4, 4]`)

**Frequency parameters:**
- `max_omega`: Maximum angular frequency ω (in rad/ps)
- `d_omega`: Frequency step Δω (in rad/ps)

**Output control:**
- `write_parameters`: Array `[wr1, wr2, wr3]` controlling which outputs to write:
  - `wr1`: Write intermediate scattering function F(q,t) (True/False)
  - `wr2`: Write dynamic structure factor S(q,ω) (True/False)
  - `wr3`: Write integrated dynamic structure factor (True/False)

**Integration ranges (optional):**
- `integration_list`: List of q-ranges for integration (e.g., `[[0.5, 2.0], [2.0, 5.0]]`)

### Example Input (input.json):

```json
{
  "property_type": 2,
  "num_atoms": 10000,
  "num_frame": 2000,
  "num_types": 2,
  "dim": 3,
  "dt": 0.001,
  "Nc": 1000,
  "compute_type": [1, 2],
  "vectors": 50,
  "q_dir": [1, 0, 0],
  "uCell": [4, 4, 4],
  "max_omega": 30,
  "d_omega": 0.1,
  "write_parameters": [true, true, true],
  "integration_list": [[1.0, 3.0], [3.0, 6.0]],
  "num_processes": -1
}
```

### Output Files:

1. **Intermediate Scattering Function** (if `wr1` is true):
   - Files: `Intermediate_scattering_{atom_type}_{q_vector_index}.txt`
   - Format: Two columns `[time, F(q,t)]`
   - One file per atom type and q-vector

2. **Dynamic Structure Factor** (if `wr2` is true):
   - Files: `dynamic_structure_{atom_type}_{q_vector_index}.txt`
   - Format: Two columns `[frequency, S(q,ω)]`
   - One file per atom type and q-vector

3. **Integrated Dynamic Structure Factor** (if `wr3` is true):
   - Files: `Integration_dynamic_structure_{atom_type}_{integration_range_index}.txt`
   - Format: Two columns `[frequency, S_integrated(ω)]`
   - One file per atom type and integration range

### Physical Interpretation:

**F(q,t) decay:**
- Fast decay → rapid structural relaxation
- Slow decay → slow dynamics (e.g., glasses, supercooled liquids)
- Exponential decay → single relaxation time
- Stretched exponential → distribution of relaxation times

**S(q,ω) features:**
- **Sharp peak**: Well-defined phonon or collective excitation
- **Broad peak**: Strongly damped excitations
- **Peak position**: Characteristic frequency of density fluctuations
- **Integrated intensity**: Total scattering strength (sum rule: ∫S(q,ω)dω = N)

**Applications:**
- Phonon dispersion relations
- Sound velocities (from low-q limit of S(q,ω))
- Diffusive dynamics (from low-ω behavior)
- Comparison with inelastic neutron/X-ray scattering

### Required Files:
- `conf.dump_all`: Includes atom positions and velocities

### Notes:
- The implementation uses time-averaging over multiple time origins for better statistics
- A Hann window function is applied before Fourier transform to reduce spectral leakage
- The code exploits the symmetry F(q,t) = F(q,-t) for efficiency
- For detailed equation verification, see `EQUATIONS_VERIFICATION.md`

## 5. Other Properties (Work in Progress)
Additional material properties are currently under development and will be added to the repository soon.

---

## ⚡ Using Multiprocessing for Faster Computation

This package supports parallel computation using Python's multiprocessing to speed up calculations when computing properties for multiple atom types.

### Quick Start

Add `num_processes` to your `input.json`:

```json
{
  "num_atoms": 10000,
  "num_frame": 1000,
  "num_types": 4,
  "compute_type": [1, 2, 3, 4],
  "num_processes": -1,
  ...
}
```

| `num_processes` | Behavior |
|-----------------|----------|
| `-1` | Use all CPU cores (recommended) |
| `1` | Sequential (no parallelization) |
| `N` | Use N cores |

### Performance Example

For a system with 10,000 atoms, 1,000 frames, and 4 atom types:
- **Sequential:** 30 minutes
- **8 cores:** 8-10 minutes (3-4x speedup)
- **16 cores:** 7-9 minutes (3.5-4x speedup)

### Documentation

- **Quick guide**: See `MULTIPROCESSING_USAGE.md`
- **Advanced HPC usage**: See `PARALLELIZATION_GUIDE.md`

---

Stay tuned for updates and improvements!
