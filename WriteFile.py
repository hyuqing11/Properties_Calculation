
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
