from WriteFile import WriteFile
import numpy as np
import sys
class String_Visualization:
    def __init__(self,pos,parameters,folder,latt_matrix,string_dictionaries):
        self.pos = pos
        self.parameters = parameters
        self.folder = folder
        self.latt_matrix = latt_matrix
        self.string_dictionaries = string_dictionaries


    def build_string_library(self,n,j):
        string_library = {}
        for tt, atom_indices in zip(self.time[n][j], self.string_atoms[n][j]):
            if len(atom_indices) >= 3:
                if tt in string_library:
                    string_library[tt].append(atom_indices)
                else:
                    string_library[tt] = [atom_indices]
        return string_library

    '''def string_visual(self,i,j,tt):
        file = '/strint_visual_' + str(i) + '_' + str(j)
        n = int((i-self.parameters['start_steps'])/self.parameters['gap'])
        string_library = self.build_string_library(n,j)
        current_pos = self.pos[tt * i + j]
        write_string = WriteFile(self.folder)
        write_string.write_string_visual(current_pos,file, string_library[tt],self.latt_matrix)
        print("checkpoint")'''
    def string_visual(self):
        for i in self.parameters['steps_read']:
            for j in self.parameters['initial_read']:
                filepath = f'{self.folder}/{i}/{j}/STRINGS'
                for tt in self.parameters['time_index']:
                    file = '/string_visual_' + str(i) + '_' + str(j) + '_' + str(tt)
                    current_pos = self.pos[tt * i + j]
                    write_string = WriteFile(self.folder)
                    if tt in self.string_dictionaries[filepath]:
                        write_string.write_string_visual(current_pos, file, self.string_dictionaries[filepath][tt], self.latt_matrix)
                    else:
                        print(f"no string forming in {filepath}  at time {tt:5d}")



