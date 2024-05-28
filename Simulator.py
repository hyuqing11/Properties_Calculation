from Read_String_Data import Read_String_Data
from String_Visualization import String_Visualization
import numpy as np
from ComputeDynamicProperties import ComputeDynamicProperties
class Simulator:
    def __init__(self, parameters,pos,vel,latt,atom_pos_dict,atom_vel_dict,latt_matrix):
        self.pos = pos
        self.vel = vel
        self.latt = latt
        self.atom_pos_dict = atom_pos_dict
        self.atom_vel_dict = atom_vel_dict
        self.parameters = parameters
        self.latt_matrix = latt_matrix

    def simulator_choice(self,folder):
        if self.parameters['property_type'] == 4:
            print("compute propensity")

        if self.parameters['property_type'] == 6:
            print("vitualize the string")
            string_reader = Read_String_Data(folder,self.parameters)
            #time, string_length, string_atoms = string_reader.read_string()
            #visual = String_Visualization(string_length,string_atoms,time,self.pos,self.parameters,folder,self.latt_matrix)
            #string_length, string_atoms, time, pos, parameters, folder, latt_matrix
            #visual.string_visual(4,0,154)
            string_libraries = string_reader.read_array_string()
            visual = String_Visualization(self.pos,self.parameters,folder,self.latt_matrix,string_libraries)
            visual.string_visual()
            print("completed")


        if self.parameters['property_type'] == 1:
            print("compute van hove self correlation")
            sz_ctype = np.size(self.parameters['compute_type'])
            for j in self.parameters['compute_type']:
                van_hov_cal = ComputeDynamicProperties(self.atom_pos_dict[j],
                                                       self.atom_vel_dict[j], self.paramters,
                                                       self.latt_matrix)
                Gr_mean, shells, r = van_hov_cal.calculate_van_hove_function(self.parameters)
        if self.parameters['property_type'] == 0:
            print('compute VACF and Pdos')
            for j in self.parameters['compute_type']:
                if compute_type[j] == num_types + 1:
                    pdos_cal = ComputeDynamicProperties(self.pos, self.vel, self.parameters, self.latt)
                    vacf_non, vacf_output, pdos = pdos_cal.pdos()
                else:
                    pdos_cal = ComputeDynamicProperties(self.atom_pos_dict[j], atom_vel_dict[j],
                                                        self.parameters, self.latt)
                    vacf_non, vacf_output, pdos = pdos_cal.pdos()







