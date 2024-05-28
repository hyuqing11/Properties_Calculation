import numpy as np
from MathFunctions import MathFunctions
import copy
from Post_Position import Post_Position
class ComputeDynamicProperties:
    def __init__(self,atom_positions,atom_velocity, parameters ,lattice,atom_type):
        self.parameters = parameters
        self.lattice = lattice
        self.atom_positions = atom_positions
        self.atom_velocity = atom_velocity
        self.atom_types = atom_type

    def spatial_correlation_time(self,num1,j,kk):
        shells = int(self.parameters["rCutOff"] / self.parameters('rDel') + 1)
        siOutput = np.zeros([num1, shells])
        pos_org = self.atom_positions[kk * self.parameters['Nc']]
        displacement_cal = Post_Position()
        dr = displacement_cal.cal_displacement(pos_org, self.atom_positions[kk * self.parameters['Nc'] + j], self.lattice)
        dr_mag = np.sqrt(np.sum(dr * dr, axis=1))
        for i in range(num1):
            if dr_mag[i] < self.parameters["rCutOff"]:
                shellNum = int(np.floor(dr_mag[i] / self.parameters('rDel')))
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
        #sz_time = np.size(specific_time)
        for kk in range(self.parameters['ave_num']):
            for jj in self.parameters['time_series']:
                meanShells, shells = self.spatial_correlation_time(self.parameters['rDel'],sz_pos[1],self.parameters['rCutOff'],jj,kk)
                Gr[kk].append(meanShells)

        Gr_mean = np.mean(np.array(Gr), axis=0)
        r = self.parameters['rDel'] *(np.array(range(shells))+0.5)
        return Gr_mean,shells,r

    def pdos(self,omega):
        # v_all: all the velocity data
        # Nc: number of correlation steps
        # dt: time interval between two frames, in units of ps
        # omega: phonon angular frequency points you want to consider
        M = self.parameters['endsteps'] - self.parameters['Nc']  # number of time origins for time average
        vacf = np.zeros(self.parameters['Nc'] )  # the velocity autocorrelation function (VACF)
        for nc in (range(self.parameters['Nc'] )):  # loop over the correlation steps

            for m in range(M + 1):  # loop over the time origins
                delta = np.sum(np.array(self.atom_velocity[m + 0]) * np.array(self.atom_velocity[m + nc]))
                vacf[nc] = vacf[nc] + delta

        vacf_non = copy.deepcopy(vacf)
        vacf = vacf / vacf[0]  # normalize the VACF
        vacf_output = copy.deepcopy(vacf)
        ff_cal = MathFunctions()
        pdos = ff_cal.compute_fourier_transform(vacf, self.parameters['Nc'] , omega, self.parameters['dt'])  # copy the VACF before modifying it

        return vacf_non, vacf_output, pdos



