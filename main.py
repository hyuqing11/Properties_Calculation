
import os
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
from ReadConfiguration import ReadConfiguration
from Read_InputFile import Read_InputFile
from Simulator import Simulator
from IntegrationDynamicStructure import IntegrationDynamicStructure
from MathFunctions import MathFunctions




def read_input_parameters(current_folder, input_filename):
    """
    Reads and validates input parameters from a file.

    Parameters:
    - current_folder: str, path to the current working directory
    - input_filename: str, name of the input file

    Returns:
    - parameters: dict, validated input parameters
    """
    read_input = Read_InputFile(current_folder)
    return read_input.read_and_validate_parameters(input_filename)


def read_configuration(parameters, current_folder):
    """
    Reads the configuration based on the system type.

    Parameters:
    - parameters: dict, input parameters
    - current_folder: str, path to the current working directory

    Returns:
    - pos: ndarray, positions of atoms
    - vel: ndarray, velocities of atoms
    - latt_matrix: ndarray, lattice matrix
    - atom_pos_dict: dict, dictionary of atomic positions
    - atom_vel_dict: dict, dictionary of atomic velocities
    - latt: ndarray, lattice parameters
    """
    pos_filename = None
    if parameters['system'] == 'lammps':
        pos_filename = os.path.join(current_folder, 'conf.dump_all')
    elif parameters['system'] == 'vasp':
        pos_filename = os.path.join(current_folder, 'XDATCAR')
    else:
        raise ValueError(f"Unknown system type: {parameters['system']}")
    config_reader = ReadConfiguration(pos_filename,parameters)

    if parameters['system'] == 'lammps':
        pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict = config_reader.read_lammps()
        latt = latt_matrix[:, 1] - latt_matrix[:, 0]
    elif parameters['system'] == 'vasp':
        pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict = config_reader.read_vasp(parameters['SysName'])
        sum_square = MathFunctions()
        latt = sum_square.sum_of_squares(latt_matrix[0])
        pos = pos * latt
        for i in range(len(atom_pos_dict)):
            atom_pos_dict[i + 1] = atom_pos_dict[i + 1] * latt

    return pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict, latt


def run_simulation(parameters, current_folder):
    """
    Runs the simulation based on input parameters.

    Parameters:
    - parameters: dict, input parameters
    - current_folder: str, path to the current working directory
    """
    if 'integrationOnly' in parameters and parameters['integrationOnly']:
        IntDyn = IntegrationDynamicStructure(current_folder, parameters)
        IntDyn._integration_dynamic_structure()
    else:
        pos, vel, latt_matrix, atom_pos_dict, atom_vel_dict, latt = read_configuration(parameters, current_folder)
        simulator = Simulator(parameters, pos, vel, latt, atom_pos_dict, atom_vel_dict, latt_matrix)
        simulator.simulator_choice(current_folder)


def main():
    """
    Main function to execute the simulation.

    Simulation options:
    - set vacf and pdos: 0
    - set van hove: 1
    - set dynamics structures: 2
    - set propensity plot: 3
    - set Cd: 4
    - set oct analysis: 5
    - set string analysis: 6
    - set dynamic facilition: 7
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting calculating the materials properties")
    current_folder = os.getcwd()
    input_filename = 'input.json'

    # Read and validate input parameters
    parameters = read_input_parameters(current_folder, input_filename)

    # Run the simulation based on input parameters
    run_simulation(parameters, current_folder)\

    logging.info("Simulation completed")


if __name__ == "__main__":
    main()