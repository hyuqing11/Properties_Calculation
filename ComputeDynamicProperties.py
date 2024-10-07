import numpy as np
from MathFunctions import MathFunctions
import copy
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
        for i in range(num1):
            if dr_mag[i] < self.parameters["rCutOff"]:
                shellNum = int(np.floor(dr_mag[i] / self.parameters['rDel']))
                if shellNum > 0:
                    rCount = siOutput[i, shellNum]
                    rCount = rCount + 1
                    siOutput[i, shellNum] = rCount

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
        vacf = np.zeros(self.parameters['Nc'] )
        for nc in (range(self.parameters['Nc'] )):

            for m in range(M + 1):  # loop over the time origins
                delta = np.sum(np.array(self.atom_velocity[m + 0]) * np.array(self.atom_velocity[m + nc]))
                vacf[nc] = vacf[nc] + delta

        vacf_non = copy.deepcopy(vacf)
        vacf = vacf / vacf[0]  # normalize the VACF
        vacf_output = copy.deepcopy(vacf)
        ff_cal = MathFunctions()
        pdos = ff_cal.compute_fourier_transform(vacf, self.parameters['Nc'] , omega, self.parameters['dt'])  # copy the VACF before modifying it

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
        fd = np.zeros((self.parameters['vectors'],self.parameters['Nc']))
        for kk in range(self.parameters['vectors']):
            c = np.sum(np.cos(np.sum(q[kk] * self.atom_positions, axis=2)), axis=1)
            s = np.sum(np.sin(np.sum(q[kk] * self.atom_positions, axis=2)), axis=1)
            for nc in range(self.parameters['Nc']):
                for m in range(M+1):
                    delta = (c[m + 0] * c[m + nc] + s[m + 0] * s[m + nc])
                    fd[kk,nc] = fd[kk,nc] + delta
        num_atoms = np.shape(self.atom_positions)
        fd_scale = fd/((M+1)*num_atoms[1])
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












