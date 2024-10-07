import json
import sys
from typing import Dict, Any
import re
class Read_InputFile:
    def __init__(self,folder):
        self.folder = folder

    def read_and_validate_parameters(self,filename) -> dict:
        input_filename = self.folder + '/' + filename
        parameters = self.read_parameters_from_json(input_filename)
        property_type = parameters['property_type']
        self.check_parameters(parameters,property_type)
        return parameters

    def read_parameters_from_json(self,filename) -> Dict[str, Any]:
        """Reads parameters from a JSON file.
        """
        try:
            with open(filename, 'r') as file:
                content = file.read()
                content_no_comments = self.remove_comments(content)
                parameters = json.loads(content_no_comments)
                if parameters.get('property_type') in [0] and parameters.get('system') == 'vasp':
                    parameters['compute_velocity'] = 1
                else:
                    parameters['compute_velocity'] = 0
            return parameters
        except FileNotFoundError:
            print(f"Error: JSON file '{filename}' not found.")
            sys.exit(1)
        except json.JSONDecodeError:
            print(f"Error: Invalid JSON format in '{filename}'.")
            sys.exit(1)


    def remove_comments(self, json_str: str) -> str:
        """Removes any lines starting with // or # and inline comments from a JSON string.
        """
        # Split the content into lines
        lines = json_str.split('\n')

        # Remove lines that start with // or #
        lines = [line for line in lines if not line.strip().startswith(("//", "#"))]

        # Join lines back into a single string
        json_str = '\n'.join(lines)

        # Remove inline comments
        json_str = re.sub(r'\/\/.*', '', json_str)  # Remove // comments
        json_str = re.sub(r'\s*#.*', '', json_str)  # Remove # comments
        return json_str

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
        """
        if property_type == 1:
            return ['system','num_atoms', 'num_frame', 'num_types', 'dim', 'dt', 'rDel', 'rCutOff', 'gap', 'time_series',
                    'compute_type','ave_num']
        elif property_type == 0:
            return ['system','num_atoms', 'num_frame', 'num_types', 'dim', 'dt', 'max_omega', 'd_omega', 'Nc',
                    'compute_type']
        elif property_type == 4:
            return ['system','num_atoms', 'num_frame', 'num_types', 'compute_type', 'q_dir','uCell']
        elif property_type == 6:
            return ['system','num_atoms', 'num_frame', 'num_types', 'dim', 'dt', 'gap', 'ave_num', 'steps_read', 'initial_read',
                    'time_index','minimum_length']
        elif property_type ==2:
            return ['system','num_atoms', 'num_frame', 'num_types', 'dim', 'dt', 'max_omega', 'd_omega', 'Nc',
                    'compute_type', 'write_parameters','vectors','q_dir','integration_list']
        elif property_type == 3:
            return ['max_omega', 'd_omega']
        elif property_type == 7:
            return ['system','num_atoms','num_frame','num_types', 'dim','start_gap','end_gap','gap','compute_type','rDel']

        # Add other property types and their required parameters here

        else:
            print(f"Error: Unsupported property type '{property_type}'.")
            sys.exit(1)



