from ReadConfiguration import ReadConfiguration
from Post_Position import Post_Position
import numpy as np
class ComputePropensity:
    def __init__(self,current_folder,parameters):
        self.current_folder = current_folder
        self.parameters = parameters

    def compute_propensity(self):
        dis_sqr = np.zeros((self.parameters['num_independent_runs'], self.parameters['num_frame'],
                            self.parameters['num_compute_atoms'], self.parameters['dim']))

        for run_id in range(self.parameters['num_independent_runs']):
            folder = f"{self.current_folder}/{run_id + 1}/"
            pos_filename = f"{folder}conf.dump_all"

            # Read atom position data
            pos_current, latt = self.read_atom_positions(pos_filename)
            pos_org = pos_current[0]

            # Compute displacement squared for each frame
            for frame_id in range(self.parameters['num_frame']):
                dis_sqr[run_id, frame_id] = self.compute_displacement_squared(pos_current[frame_id], pos_org, latt)

        return dis_sqr, pos_org, latt

    def compute_displacement_squared(self, pos_current, pos_org, latt):
        displacement_cal = Post_Position()
        dr = displacement_cal.cal_displacement(pos_current, pos_org, latt)
        return dr ** 2

    def read_atom_positions(self, pos_filename):
        config_reader = ReadConfiguration(pos_filename, self.parameters['num_frame'],
                                          self.parameters['num_types'], self.parameters['num_atoms'],
                                          self.parameters['dim'])
        pos, _, latt_matrix, atom_pos_dict, _ = config_reader.Read_lammps()

        latt = latt_matrix[:, 1] - latt_matrix[:, 0]
        pos_current = atom_pos_dict[self.parameters['compute_type']]
        return pos_current, latt

    def compute_delta_propensity(self,dis_sqr_sum):
        x_mean = np.mean(dis_sqr_sum, axis=0)

        x_mean_t = np.mean(x_mean, axis=1)
        x_mean_t_shape = x_mean_t.reshape(800, 1)
        delta_x = x_mean - x_mean_t_shape
        delta_x_t = np.mean(delta_x ** 2, axis=1)
        return delta_x, delta_x_t

    def spatial_correlation_time(self, pos, lattice, delta_w):

        shells = int(self.parameters['rCut'] / self.parameters['rDel']) + 1
        si_output = np.zeros((self.parameters['num_compute_atoms'], shells))
        displacement_cal = Post_Position()

        for i in range(self.parameters['num_compute_atoms'] - 1):
            for j in range(i + 1, self.parameters['num_compute_atoms']):
                dij, shell_num = self.compute_interatomic_distance(pos[i], pos[j], lattice, displacement_cal)
                if dij < self.parameters['rCut'] and shell_num > 0:
                    si_output[i, shell_num] += delta_w[i] * delta_w[j]

        mean_shell = si_output.mean(axis=0)
        return mean_shell, shells

    def compute_interatomic_distance(self, atom_i, atom_j, lattice, displacement_cal):
        dx, dy, dz = atom_j - atom_i
        dx = displacement_cal.wrap(dx, lattice[0])
        dy = displacement_cal.wrap(dy, lattice[1])
        dz = displacement_cal.wrap(dz, lattice[2])

        dij = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        shell_num = int(np.floor(dij / self.parameters['rDel']))
        return dij, shell_num

    def compute_spatial_corr(self):
        dis_sqr, pos_org, latt = self.compute_propensity()
        dis_sqr_sum = np.sum(dis_sqr, axis=3)
        delta_w, delta_w_mean = self.compute_delta_propensity(dis_sqr_sum)

        Cd_record = []
        for frame_id in range(self.parameters['num_frame']):
            mean_shell, _ = self.spatial_correlation_time(pos_org, latt, delta_w[frame_id])
            Cd_record.append(mean_shell / delta_w_mean[frame_id])

        return Cd_record





