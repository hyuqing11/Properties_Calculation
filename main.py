import numpy as np
import sys
import json
import os
from ReadConfiguration import ReadConfiguration
from Read_InputFile import Read_InputFile
from Simulator import Simulator

def main():
    """
    Simulation options:
    - set vacf: 0
    - set van hove: 1
    - set dynamics structures: 2
    - set propensity plot: 3
    - set Cd: 4
    - set oct analysis: 5
    - set string analysis: 6
    """
    current_folder = '/Users/hyuqing/Downloads/research/YH2/YHx/YH1.98/msd/900/1/'
    input_filename = 'input.json'

    # Read and validate input parameters
    read_input = Read_InputFile(current_folder)
    parameters = read_input.read_and_validate_parameters(input_filename)

    # Read configuration from LAMMPS dump file
    pos_filename = os.path.join(current_folder, 'conf.dump_all')
    config_reader = ReadConfiguration(
        pos_filename,
        parameters['num_frame'],
        parameters['num_types'],
        parameters['num_atoms'],
        parameters['dim']
    )
    pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict = config_reader.read_lammps()
    latt = latt_matrix[:, 1] - latt_matrix[:, 0]

    # Run simulation based on input parameters
    simulator = Simulator(parameters, pos, vel, latt, atom_pos_dict, atom_vel_dict, latt_matrix)
    simulator.simulator_choice(current_folder)

    print("Simulation completed")

if __name__ == "__main__":
    main()
