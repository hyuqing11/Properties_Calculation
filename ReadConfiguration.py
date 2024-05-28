import numpy as np
from typing import Tuple, Dict, List


class ReadConfiguration:
    def __init__(self, filename: str, num_frame: int, num_types: int, num_atoms: int, dim: int):
        self.filename = filename
        self.num_frame = num_frame
        self.num_types = num_types
        self.num_atoms = num_atoms
        self.dim = dim

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
