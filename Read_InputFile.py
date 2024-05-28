import json
import sys
from typing import Dict, Any
class Read_InputFile:
    def __init__(self,folder):
        self.folder = folder

    def read_and_validate_parameters(self,filename,property_type=0) -> dict:
        input_filename = self.folder + '/' + filename
        parameters = self.read_parameters_from_json(input_filename)
        property_type = parameters['property_type']
        self.check_parameters(parameters,property_type)
        return parameters

    def read_parameters_from_json(self,filename) -> Dict[str, Any]:
        """Reads parameters from a JSON file.

        Returns:
            A dictionary containing the parameters.

        Raises:
            SystemExit: If the file is not found or contains invalid JSON.
        """

        try:
            with open(filename, 'r') as file:
                parameters = json.load(file)
            return parameters
        except FileNotFoundError:
            print(f"Error: JSON file '{filename}' not found.")
            sys.exit(1)
        except json.JSONDecodeError:
            print(f"Error: Invalid JSON format in '{filename}'.")
            sys.exit(1)

    def check_parameters(self, parameters: Dict[str, Any], property_type: int):
        """Checks if all required parameters are present for a given property type.

        Args:
            parameters: A dictionary of parameters read from the JSON file.
            property_type: An integer representing the property type.

        Raises:
            SystemExit: If any required parameters are missing.
        """
        required_params = self._get_required_params(property_type)
        missing_params = [param for param in required_params if param not in parameters]

        if missing_params:
            print("Error: The following parameters are missing from the JSON file:")
            for param in missing_params:
                print(f"- {param}")
            sys.exit(1)

    def _get_required_params(self, property_type: int) -> list:
        """Returns a list of required parameter keys based on the property type.

        Args:
            property_type: An integer representing the property type.

        Returns:
            A list of required parameter keys.
        """
        if property_type == 1:
            return ['num_atoms', 'num_frame', 'num_types', 'dim', 'dt', 'rDel', 'rCutOff', 'Nc', 'time_series',
                    'output_filename', 'compute_type']
        elif property_type == 0:
            return ['num_atoms', 'num_frame', 'num_types', 'dim', 'dt', 'ave_num', 'max_omega', 'd_omega', 'Nc',
                    'output_filename', 'compute_type']
        elif property_type == 4:
            return ['num_atoms', 'num_frame', 'num_types', 'dim', 'dt', 'gap', 'compute_type', 'start_steps','cutoff_distance','rDel']
        elif property_type ==6:
            #return ['num_atoms', 'num_frame', 'num_types', 'dim', 'dt', 'gap', 'ave_num', 'start_steps','end_steps','time_index']
            return ['num_atoms', 'num_frame', 'num_types', 'dim', 'dt', 'gap', 'ave_num', 'steps_read', 'initial_read',
                    'time_index','minimum_length']

        # Add other property types and their required parameters here

        else:
            print(f"Error: Unsupported property type '{property_type}'.")
            sys.exit(1)

