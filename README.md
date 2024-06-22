# Material Properties Calculation

This repository contains code for calculating various material properties using data obtained from molecular dynamics simulations (LAMMPS). Currently, the following calculations are supported:

## 1. Plot Strings in OVITO (property_type = 6)
This option allows you to visualize the string configuration using OVITO software. Strings are accessed by connecting mobile atoms i and j if:
$$|\vec{r}_i(t)-\vec{r}_j(t)|<\delta$$
or
$$|\vec{r}_j(t)-\vec{r}_i(t)|<\delta$$

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

## 4. Imtermeidate Scattering Function (property_type=2)
This opetion calculates the imtermediate scattering function and dynamics structure factors for materials

## 5. Other Properties (Work in Progress)
Additional material properties are currently under development and will be added to the repository soon.

Stay tuned for updates and improvements!
