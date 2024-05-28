import os
import numpy as np
from Read_String_Data import Read_String_Data
from String_Visualization import String_Visualization
from ComputeDynamicProperties import ComputeDynamicProperties
from WriteFile import WriteFile


class Simulator:
    def __init__(self, parameters, pos, vel, latt, atom_pos_dict, atom_vel_dict, latt_matrix):
        self.parameters = parameters
        self.pos = pos
        self.vel = vel
        self.latt = latt
        self.atom_pos_dict = atom_pos_dict
        self.atom_vel_dict = atom_vel_dict
        self.latt_matrix = latt_matrix

    def simulator_choice(self, folder):
        property_type = self.parameters.get('property_type')

        if property_type == 4:
            self._compute_propensity()
        elif property_type == 6:
            self._visualize_string(folder)
        elif property_type == 1:
            self._compute_van_hove()
        elif property_type == 0:
            self._compute_vacf_and_pdos()
        else:
            print(f"Unknown property type: {property_type}")

    def _compute_propensity(self):
        print("Compute propensity")

    def _visualize_string(self, folder):
        print("Visualize the string")
        string_reader = Read_String_Data(folder, self.parameters)
        string_libraries = string_reader.read_array_string()
        visual = String_Visualization(self.pos, self.parameters, folder, self.latt_matrix, string_libraries)
        visual.string_visual()
        print("Visualization completed")

    def _compute_van_hove(self):
        print("Compute van Hove self-correlation")
        for j in self.parameters['compute_type']:
            van_hov_cal = ComputeDynamicProperties(self.atom_pos_dict[j], self.atom_vel_dict[j], self.parameters,
                                                   self.latt_matrix)
            Gr_mean, shells, r = van_hov_cal.calculate_van_hove_function(self.parameters)

            sz_time = np.size(self.parameters['time_series'])
            for i in range(sz_time):
                output_filename = f'atom_Gs_{j}_{self.parameters["time_series"][i]}.txt'
                van_hov_write = WriteFile(self.folder)
                van_hov_write.write_results(output_filename, r, Gr_mean[i])

    def _compute_vacf_and_pdos(self):
        print('Compute VACF and Pdos')
        nu = self.parameters['omega'] / (2 * np.pi)
        dt = self.parameters['dt']
        t = np.arange(self.parameters['Nc']) * dt

        for j in self.parameters['compute_type']:
            if j == self.parameters['compute_type'] + 1:
                pdos_cal = ComputeDynamicProperties(self.pos, self.vel, self.parameters, self.latt)
            else:
                pdos_cal = ComputeDynamicProperties(self.atom_pos_dict[j], self.atom_vel_dict[j], self.parameters,
                                                    self.latt)

            vacf_non, vacf_output, pdos = pdos_cal.pdos()

            self._write_vacf_and_pdos(j, nu, pdos, t, vacf_non, vacf_output)

    def _write_vacf_and_pdos(self, atom_type, nu, pdos, t, vacf_non, vacf_output):
        vacf_write = WriteFile(self.folder)

        output_filename = f'atom{atom_type}_pdos_{self.parameters["Nc"]}.txt'
        vacf_write.write_properties(output_filename, nu, pdos)

        output_filename = f'atom{atom_type}_vacf_non_{self.parameters["Nc"]}.txt'
        vacf_write.write_properties(output_filename, t, vacf_non)

        output_filename = f'atom{atom_type}_vacf_scale_{self.parameters["Nc"]}.txt'
        vacf_write.write_properties(output_filename, t, vacf_output)

