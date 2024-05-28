import os
import re
import numpy as np
from typing import List, Tuple
class Read_String_Data:
    def __init__(self,folder,parameters):
        self.folder = folder
        self.parameters = parameters
        '''self.start_time = start_time
        self.end_time = end_time
        self.gap = gap
        self.ave_num = ave_num'''
        self.pattern = r'Time= (\d+)\s+(\d+)\s+(.+)'

    def _initialize_storage(self) -> Tuple[List[List[List[int]]], List[List[List[int]]], List[List[List[np.ndarray]]]]:
        """Initialize nested lists for time, string length, and string atoms."""
        steps = (self.parameters['end_steps'] - self.parameters['start_steps']) // self.parameters['gap']
        return (
            [[[] for _ in range(self.parameters['ave_num'])] for _ in range(steps)],
            [[[] for _ in range(self.parameters['ave_num'])] for _ in range(steps)],
            [[[] for _ in range(self.parameters['ave_num'])] for _ in range(steps)]
        )

    def _parse_line(self, line: str, pattern: str) -> Tuple[int, int, np.ndarray]:
        """Parse a line using the given regex pattern."""
        match = re.match(pattern, line)
        if match:
            time = int(match.group(1))
            length = int(match.group(2))
            atoms = np.array([int(x) for x in match.group(3).split('-->') if x.strip().isdigit()])
            return time, length, atoms
        return None

    def _read_file(self, filepath: str) -> List[str]:
        """Read lines from a file if it exists."""
        if os.path.exists(filepath):
            with open(filepath) as file:
                return file.readlines()
        print(f'{filepath} does not exist')
        return []

    def read_string_file(self,filename):
        string_length = []
        string_atom  = []
        string_time = []
        lines = self._read_file(filename)
        for line in lines:
            parsed_data = self._parse_line(line.strip(), self.pattern)
            if parsed_data:
                t, length, atoms = parsed_data
                string_time.append(t)
                string_length.append(length)
                string_atom.append(atoms)

        return string_time,string_length,string_atom


    def build_string_library(self,string_time, string_atom):
        string_library = {}
        for tt, atom_indices in zip(string_time, string_atom):
            if len(atom_indices) >= 3:
                if tt in string_library:
                    string_library[tt].append(atom_indices)
                else:
                    string_library[tt] = [atom_indices]
        return string_library

    def read_array_string(self):
        string_dictionaries = {}
        for i in self.parameters['steps_read']:
            for j in self.parameters['initial_read']:
                filepath = f'{self.folder}/{i}/{j}/STRINGS'
                string_time,string_length,string_atom  = self.read_string_file(filepath)
                string_dictionaries[filepath] = self.build_string_library(string_time,string_atom)
        return string_dictionaries


    def read_string(self) -> Tuple[List[List[List[int]]], List[List[List[int]]], List[List[List[np.ndarray]]]]:
        pattern = r'Time= (\d+)\s+(\d+)\s+(.+)'
        time, string_lengths, string_atoms = self._initialize_storage()

        for n, i in enumerate(range(self.parameters['start_steps'], self.parameters['end_steps'], self.parameters['gap'])):
            for j in range(self.parameters['ave_num']):
                filepath = f'{self.folder}/{i}/{j}/STRINGS'
                lines = self._read_file(filepath)
                for line in lines:
                    parsed_data = self._parse_line(line.strip(), pattern)
                    if parsed_data:
                        t, length, atoms = parsed_data
                        time[n][j].append(t)
                        string_lengths[n][j].append(length)
                        string_atoms[n][j].append(atoms)

        return time, string_lengths, string_atoms
