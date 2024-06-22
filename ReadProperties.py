import numpy as np
import os
class ReadProperties:
    def __init__(self,folder):
        self.folder = folder

    def convert_lines_to_data(self,lines):
        data = []
        for line in lines:
            line = line.strip().split()
            row = list(map(float, line))
            data.append(row)
        return data

    def read_properties(self, filename):
        fileName = os.path.join(self.folder, filename)
        # Check if the file exists
        if not os.path.exists(fileName):
            raise FileNotFoundError(f"The file {fileName} does not exist.")

        with open(fileName) as fr:
            lines = fr.readlines()
            data = self.convert_lines_dat(lines)
        x = np.array(data)[:,0]
        y = np.array(data)[:,1:]

        return x, y