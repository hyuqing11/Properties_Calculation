import numpy as np
from typing import Tuple, Dict, List
from MathFunctions import MathFunctions

class ReadConfiguration:
    def __init__(self, filename: str,parameters):
        self.filename = filename
        self.parameters = parameters
        self.num_frame = parameters['num_frame']
        self.num_types = parameters['num_types']
        self.num_atoms = parameters['num_atoms']
        self.dim = parameters['dim']

    def read_lammps(self) -> Tuple[
        np.ndarray, np.ndarray, np.ndarray, Dict[int, List[List[np.ndarray]]], Dict[int, List[List[np.ndarray]]]]:
        """Reads LAMMPS dump file and extracts positions, velocities, and lattice information.

        Optimized version: Uses batch reading and vectorized operations for better performance
        with large atom files.

        Returns:
            pos: Positions array of shape (num_frame, num_atoms, dim).
            vel: Velocities array of shape (num_frame, num_atoms, dim).
            latt_matrix: Lattice matrix of shape (3, 2).
            atom_pos_dict: Dictionary of atom positions indexed by atom type.
            atom_vel_dict: Dictionary of atom velocities indexed by atom type.
        """
        pos = np.zeros((self.num_frame, self.num_atoms, self.dim))
        vel = np.zeros((self.num_frame, self.num_atoms, self.dim))
        latt_matrix = np.zeros((3, 2))

        # Pre-allocate array for atom types
        atom_types = np.zeros((self.num_frame, self.num_atoms), dtype=int)

        with open(self.filename, "r") as fin:
            for frame in range(self.num_frame):
                # Skip headers (5 lines)
                for _ in range(5):
                    next(fin)

                # Read lattice matrix (3 lines)
                for jj in range(3):
                    line = fin.readline().split()
                    latt_matrix[jj] = [float(value) for value in line]

                next(fin)  # Skip the blank line

                # Read all atom data for this frame at once (batch reading)
                atom_lines = []
                for _ in range(self.num_atoms):
                    atom_lines.append(fin.readline())

                # Join lines and parse as a single block using NumPy
                atom_data_str = ''.join(atom_lines)
                atom_data = np.fromstring(atom_data_str, sep=' ').reshape(self.num_atoms, -1)

                # Extract positions, velocities, and atom types using array slicing
                pos[frame] = atom_data[:, 1:4]
                vel[frame] = atom_data[:, 4:7]
                atom_types[frame] = atom_data[:, -1].astype(int)

        # Build atom dictionaries using pre-computed atom types (more efficient)
        atom_pos_dict = {i: [] for i in range(1, self.num_types + 1)}
        atom_vel_dict = {i: [] for i in range(1, self.num_types + 1)}

        for frame in range(self.num_frame):
            for atom_type in range(1, self.num_types + 1):
                # Use boolean indexing to select atoms of specific type
                type_mask = atom_types[frame] == atom_type
                atom_pos_dict[atom_type].append(pos[frame, type_mask])
                atom_vel_dict[atom_type].append(vel[frame, type_mask])

        return pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict

    def displacement(self,data1,data2):
        dr = np.subtract(data1, data2)
        dr = dr - np.round(dr)
        return dr
    def compute_velocity(self,pos, vel, latt_matrix, atom_vel_dict,elements,atom_counts):
        dr_car = np.zeros((self.num_frame - 1, self.num_atoms, 3))
        dr = np.zeros((self.num_frame - 1,self.num_atoms, 3))
        if len(latt_matrix) > 1:
            print(
                "Warning: The lattice constant changes over time, and the velocity calculation will be computed based on the previous time crystal framework")
            npt = 1
        else:
            npt = 0
        for j in range(1,self.num_frame):
            if npt:
                latt_prev = latt_matrix[j - 1]
                latt_curr = latt_matrix[j]
            else:
                latt_prev = latt_matrix[0]
                latt_curr = latt_matrix[0]
            trans = np.transpose(
            np.dot(np.transpose(np.linalg.inv((latt_prev))), (np.dot(np.transpose(latt_curr), np.transpose(pos[j])))))
            dr[j - 1] = self.displacement(trans, pos[j - 1])
            dr_car[j - 1] = np.transpose(np.dot(np.transpose(latt_prev), np.transpose(dr[j - 1])))
            vel[j - 1] = dr_car[j - 1] / self.parameters['dt']
            start_i = 0
            for i in range(len(elements)):
                end_i = start_i + atom_counts[i]
                atom_vel_dict[i + 1].append(vel[j-1, start_i:end_i])
                start_i = end_i
        #if len(latt_matrix) == 1:

    def read_vasp(self,sysName):
        with open(self.filename,'r') as f:
            contcar = f.readlines()
        pos = np.zeros((self.num_frame,self.num_atoms,self.dim))
        vel = np.zeros((self.num_frame,self.num_atoms,self.dim))
        elements = contcar[5].strip().split()
        atom_counts = list(map(int,contcar[6].strip().split()))
        atom_pos_dict = {i: [] for i in range(1,len(elements)+1)}
        atom_vel_dict = {i: [] for i in range(1,len(elements)+1)}
        latt_matrix = []
        n = 0

        for line_i, line in enumerate(contcar):
            if sysName in line:
                matrix1 = self.convert_lines_to_data(contcar[line_i + 2: line_i + 2 + self.dim])
                latt_matrix.append(matrix1)
            if "Direct" in line:
                pos[n] = self.convert_lines_to_data(contcar[line_i + 1: line_i + 1 + self.num_atoms ])
                start_i = 0
                for i in range(len(elements)):
                    end_i = start_i  + atom_counts[i]
                    atom_pos_dict[i+1].append(pos[n, start_i:end_i])
                    start_i = end_i
                n += 1
                if n + 1 > self.num_frame:
                    break
        if  self.parameters['compute_velocity']:
            self.compute_velocity(pos, vel, latt_matrix, atom_vel_dict,elements,atom_counts)
        return pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict


    def convert_lines_to_data(self,lines):
        """Convert lines to numpy array efficiently using vectorized parsing."""
        # Join all lines and use numpy's faster string parsing
        lines_str = ''.join(lines)
        # Count columns from first line
        first_line = lines[0].split()
        num_cols = len(first_line)
        # Parse all at once
        data = np.fromstring(lines_str, sep=' ').reshape(len(lines), num_cols)
        return data

    def _skip_lines(self, file, num_lines):
        for _ in range(num_lines):
            next(file)

