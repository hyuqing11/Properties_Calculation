import numpy as np
class WriteFile:
    def __init__(self,folder):
        self.folder = folder
    def write_string_visual(self,pos,file, atom_indices,latt_matrix):
        current_folder = self.folder
        filename = current_folder + '/'  + file
        n = 0
        n1 = 0
        total_lines = sum(len(sublist) for sublist in atom_indices)
        with open(filename, 'w') as fw:
            first_row = f'ITEM: TIMESTEP \n 0 \nITEM: NUMBER OF ATOMS \n {total_lines} \nITEM: BOX BOUNDS pp pp pp \n ' \
                        f'{latt_matrix[0][0]} {latt_matrix[0][1]} \n {latt_matrix[1][0]} {latt_matrix[1][1]} \n ' \
                        f'{latt_matrix[2][0]} {latt_matrix[2][1]} \nITEM: ATOMS id type x y z\n'
            fw.write(first_row)
            for atom_index in atom_indices:
                sz_index = np.size(atom_index)
                for i in range(sz_index):
                    fw.write(
                        f'{n1 + 1:5d} {n+1:5d} \t {pos[atom_index[i] - 1][0]:.8f} \t {pos[atom_index[i] - 1][1]:.8f} \t {pos[atom_index[i] - 1][2]:.8f}\n')
                    n1 += 1
                n += 1





