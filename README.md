# Material Properties Calculation

This repository contains code for calculating various material properties using data obtained from molecular dynamics simulations (LAMMPS). Currently, the following calculations are supported:

1. **Plot Strings in OVITO (property_type = 6):**
   This option allows you to visualize the string configuration using OVITO software.

2. **Van Hove Self Correlation (property_type = 1):**
   The Van Hove self-correlation function quantifies the probability of finding an ion within a volume element dr, centered at r and time t, given that the ion was at the origin at time t=0. The formula is given by:
   $$G_s(\vec{r},t)=<\frac{1}{N} \sum_{i=1}^{N} \delta(\vec{r}-|\vec{r}_i(t)-\vec{r}_i(0)|)>$$

3. **VACF and PDOS (property_type = 0):**
   This option calculates the Velocity AutoCorrelation Function (VACF) and the Phonon Density of States (PDOS) for the material.
   
   - **VACF (Velocity AutoCorrelation Function):**
     VACF provides information about particle motion and diffusion. It is defined as:
     $$C(t) = \frac{1}{N} \sum_{i=1}^{N} \vec{v}_i(t) \cdot \vec{v}_i(0)$$
     Where $\vec{v}_i(t)$ is the velocity of particle i at time t, and N is the total number of particles.

   - **PDOS (Phonon Density of States):**
     PDOS gives insights into the material's vibrational modes. It is calculated using Fourier transform techniques and provides the distribution of vibrational frequencies in the material.

4. **Other Properties (Work in Progress):**
   Additional material properties are currently under development and will be added to the repository soon.

Stay tuned for updates and improvements!


