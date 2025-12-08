import numpy as np
from MathFunctions import MathFunctions
from Post_Position import Post_Position
class ComputeDynamicProperties:
    def __init__(self,atom_positions,atom_velocity, parameters ,lattice):
        self.parameters = parameters
        self.lattice = lattice
        self.atom_positions = atom_positions
        self.atom_velocity = atom_velocity


    def spatial_correlation_time(self,num1,j,kk):
        shells = int(self.parameters["rCutOff"] / self.parameters['rDel'] + 1)
        siOutput = np.zeros([num1, shells])
        pos_org = self.atom_positions[kk * self.parameters['gap']]
        displacement_cal = Post_Position()
        dr = displacement_cal.cal_displacement(pos_org, self.atom_positions[kk * self.parameters['gap'] + j], self.lattice)
        dr_mag = np.sqrt(np.sum(dr * dr, axis=1))

        # Vectorized version - avoid explicit loop
        # Create boolean mask for valid shells
        valid_mask = dr_mag < self.parameters["rCutOff"]
        shellNums = np.floor(dr_mag / self.parameters['rDel']).astype(int)

        # Only process atoms within cutoff and with shellNum > 0
        valid_indices = np.where(valid_mask & (shellNums > 0))[0]

        # Use advanced indexing to increment counts
        for i in valid_indices:
            siOutput[i, shellNums[i]] += 1

        mean_Shells = siOutput.mean(0)
        return mean_Shells, shells


    def calculate_van_hove_function(self):
        sz_pos = np.shape(self.atom_positions)
        Gr = []
        for kk in range(self.parameters['ave_num']):
            Gr.append([])
        for kk in range(self.parameters['ave_num']):
            for jj in self.parameters['time_series']:
                meanShells, shells = self.spatial_correlation_time(sz_pos[1],jj,kk)
                Gr[kk].append(meanShells)

        Gr_mean = np.mean(np.array(Gr), axis=0)
        r = self.parameters['rDel'] *(np.array(range(shells))+0.5)
        return Gr_mean,shells,r

    def pdos(self,omega):
        # v_all: all the velocity data
        # Nc: number of correlation steps
        # dt: time interval between two frames, in units of ps
        # omega: phonon angular frequency points you want to consider
        M = self.parameters['num_frame'] - self.parameters['Nc']

        # Convert velocity list to numpy array once (avoid repeated conversions)
        vel_array = np.array(self.atom_velocity)

        # Vectorized VACF calculation - much faster than nested loops
        # Using broadcasting to compute all correlations at once
        vacf = np.zeros(self.parameters['Nc'])

        for nc in range(self.parameters['Nc']):
            # Vectorized dot product over all time origins m
            # Shape: (M+1, num_atoms, 3) dot (M+1, num_atoms, 3) -> (M+1,)
            correlations = np.sum(vel_array[0:M+1] * vel_array[nc:M+1+nc], axis=(1, 2))
            vacf[nc] = np.sum(correlations)

        # Store unnormalized version (avoid deep copy)
        vacf_non = vacf.copy()

        # Normalize the VACF
        vacf = vacf / vacf[0]
        vacf_output = vacf.copy()

        ff_cal = MathFunctions()
        pdos = ff_cal.compute_fourier_transform(vacf, self.parameters['Nc'], omega, self.parameters['dt'])

        return vacf_non, vacf_output, pdos

    def compute_q_vectors(self):
        q = np.zeros((self.parameters['vectors'], 3))
        #print(self.parameters['q_dir'])
        #print(np.where(self.parameters['q_dir'] != 0))
        non_zero_index = np.where(np.array(self.parameters['q_dir']) != 0)[0][0]
        for i in range(self.parameters['vectors']):
            q[i] = np.array(self.parameters['q_dir'] )* (i + 1) / self.parameters['uCell'][non_zero_index] * 2 * np.pi / \
                   self.lattice[non_zero_index]
        return q


    def calculate_intermediate_scattering(self):
        q = self.compute_q_vectors()
        M = self.parameters['num_frame'] - self.parameters['Nc']
        fd = np.zeros((self.parameters['vectors'], self.parameters['Nc']))

        # Vectorized computation - eliminates inner two nested loops
        for kk in range(self.parameters['vectors']):
            # Compute cos and sin for all frames at once
            q_dot_r = np.sum(q[kk] * self.atom_positions, axis=2)
            c = np.sum(np.cos(q_dot_r), axis=1)
            s = np.sum(np.sin(q_dot_r), axis=1)

            # Vectorized correlation calculation over all nc and m
            # Instead of nested loops, use array slicing and broadcasting
            for nc in range(self.parameters['Nc']):
                # Compute all m values at once using vectorized operations
                # c[0:M+1] and c[nc:M+1+nc] are arrays of shape (M+1,)
                fd[kk, nc] = np.sum(c[0:M+1] * c[nc:M+1+nc] + s[0:M+1] * s[nc:M+1+nc])

        num_atoms = np.shape(self.atom_positions)
        fd_scale = fd / ((M + 1) * num_atoms[1])
        return fd_scale

    def calculate_dynamic_structure(self, omega, fd_scale):
        fft = MathFunctions()
        Sv = np.zeros((self.parameters['vectors'],len(omega)))
        for i in range(self.parameters['vectors']):
            Sv[i] = fft.compute_fourier_transform(fd_scale[i], self.parameters['Nc'], omega,self.parameters['dt'])
        return Sv

    def Integrate_dynamic_structure(self, Sv):

        q = self.compute_q_vectors()
        non_zero_index = np.where(np.array(self.parameters['q_dir']) != 0)[0][0]
        num_integration = len(self.parameters['integration_list'])
        num_Sv = np.shape(Sv)
        S_int_record = np.zeros((num_integration, num_Sv[1]))
        for n, (q_min, q_max) in enumerate(self.parameters['integration_list']):
            for i in range(self.parameters['vectors']):
                if (q[i,non_zero_index] > q_min) & (q[i,non_zero_index] < q_max):
                    S_int_record[n] += Sv[i]
        return S_int_record









