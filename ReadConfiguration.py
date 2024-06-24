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
        atom_pos_dict = {i: [] for i in range(1, self.num_types + 1)}
        atom_vel_dict = {i: [] for i in range(1, self.num_types + 1)}

        with open(self.filename, "r") as fin:
            for frame in range(self.num_frame):
                # Skip headers
                for _ in range(5):
                    next(fin)

                # Read lattice matrix
                for jj in range(3):
                    line = fin.readline().split()
                    latt_matrix[jj] = [float(value) for value in line]

                next(fin)  # Skip the blank line

                # Initialize storage for this frame
                for atom_type in range(1, self.num_types + 1):
                    atom_pos_dict[atom_type].append([])
                    atom_vel_dict[atom_type].append([])

                # Read positions and velocities
                for atom in range(self.num_atoms):
                    line = fin.readline().split()
                    values = [float(value) for value in line]
                    pos[frame, atom] = values[1:4]
                    vel[frame, atom] = values[4:7]
                    atom_type = int(values[-1])
                    atom_pos_dict[atom_type][-1].append(pos[frame, atom])
                    atom_vel_dict[atom_type][-1].append(vel[frame, atom])

        return pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict

    def displacement(self,data1,data2):
        dr = np.subtract(data1, data2)
        dr = dr - np.round(dr)
        return dr
    def compute_velocity(self,pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict,,elements,atom_counts):
        dr_car = np.zeros((self.num_frame - 1, self.num_atoms, 3))
        dr = np.zeros((self.num_frame - 1,self.num_atoms, 3))
        if len(latt_matrix) > 1:
            print(
                "Warning: The lattice constant changes over time, and the velocity calculation will be computed based on the previous time crystal framework")
            npt = 1
        else:
            npt = 0
        for j in range(1,self.num_frame):
            latt = latt_matrix[j - 1] if npt else latt_matrix[0]
            trans = np.transpose(
                np.dot(np.transpose(np.linalg.inv(latt)),
                       (np.dot(np.transpose(latt), np.transpose(pos[j])))))
            dr[j - 1] = self.displacement(trans, pos[j - 1])
            dr_car[j - 1] = np.transpose(np.dot(np.transpose(latt), np.transpose(dr[j - 1])))
            vel[j - 1] = dr_car[j - 1] / self.parameters['dt']
            for i in range(len(elements)):
                end_i = start_i + atom_counts[i]
                atom_vel_dict[i + 1].append(vel[j-1, start_i:end_i])
                start_i = end_i
        #if len(latt_matrix) == 1:

    def read_vasp(self,sysName,compute_velocity):
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
        if compute_velocity:
            self.compute_velocity(pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict,elements,atom_counts)
        return pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict


    def convert_lines_to_data(self,lines):
        data = []
        for line in lines:
            tem = line.split()
            row = list(map(float, tem))
            data.append(row)
        return data

    def _skip_lines(self, file, num_lines):
        for _ in range(num_lines):
            next(file)

