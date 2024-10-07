import numpy as np
class WriteFile:
    def __init__(self, folder, parameters):
        self.folder = folder
        self.parameters = parameters

    def write_string_visual(self, pos, file, atom_indices, latt_matrix):
        current_folder = self.folder
        filename = f"{current_folder}/{file}"
        total_lines = sum(len(sublist) for sublist in atom_indices)

        with open(filename, 'w') as fw:
            first_row = (
                "ITEM: TIMESTEP \n0\nITEM: NUMBER OF ATOMS \n"
                f"{total_lines}\nITEM: BOX BOUNDS pp pp pp \n"
                f"{latt_matrix[0][0]} {latt_matrix[0][1]}\n"
                f"{latt_matrix[1][0]} {latt_matrix[1][1]}\n"
                f"{latt_matrix[2][0]} {latt_matrix[2][1]}\n"
                "ITEM: ATOMS id type x y z\n"
            )
            fw.write(first_row)

            n1 = 0
            for n, atom_index in enumerate(atom_indices, start=1):
                for idx in atom_index:
                    atom_id = idx - 1
                    x, y, z = pos[atom_id]
                    fw.write(f"{n1 + 1:5d} {n:5d} \t {x:.8f} \t {y:.8f} \t {z:.8f}\n")
                    n1 += 1

    def write_results(self, filename, x, y):
        output_filename = f"{self.folder}/{filename}"

        with open(output_filename, "w") as fw:
            for xi, yi in zip(x, y):
                fw.write(f"{xi:.8f} \t {yi:.16f}\n")


    def write_cd(self,Cd, folder,time):
        shells = int(self.parameters['rCut'] / self.parameters['rDel']) + 1
        r = np.array(range(shells)) * self.parameters['rDel']
        sz1, sz2 = np.shape(Cd)
        for i in range(sz2):
            filename = 'Cd_' + str(r[i]) + '.txt'
            with open(folder + filename, 'w') as fw1:
                for j in range(sz1):
                    fw1.write(f'{time[j]:.8f} \t {Cd[j][i]:.8f}\n')


    def write_facilitation(self,r, mean_Shells_record,mean_Shells1_record):
        sz_gap = len(mean_Shells1_record)
        size_r = np.size(r)
        gap_array = range(self.parameters['start_gap'], self.parameters['end_gap'], self.parameters['gap'])
        for i in range(sz_gap):
            with open(self.folder + 'Prob_ion_' + str(gap_array[i]), 'w')as fw:
                for j in range(size_r):
                    fw.write(f'{r[j]:.8f} \t {mean_Shells_record[i][j]:.8f} \t {mean_Shells1_record[i][j]:.8f} \n')
