import numpy as np
from ReadProperties import ReadProperties
from WriteFile import WriteFile
class IntegrationDynamicStructure:
    def __init__(self, folder, parameters,lattice):
        self.folder = folder
        self.parameters = parameters
        self.lattice = lattice
    def _integration_dynamic_structure(self):
        read_dynstructure = ReadProperties(self.folder)

        # Compute q vectors and non-zero indices
        q, non_zero_index = self.compute_q_vectors()

        # Initialize the integration record array
        num_nu = len(read_dynstructure.read_properties(self._generate_filename(1))[0])
        num_integration = len(self.parameters['integration_list'])
        S_int_record = np.zeros((num_integration, num_nu))
        # Iterate over integration list and vectors to populate S_int_record
        for n, (q_min, q_max) in enumerate(self.parameters['integration_list']):
            for i in range(self.parameters['vectors']):
                filename = self._generate_filename(i)
                nu, Sv = read_dynstructure.read_properties(filename)

                if self._is_within_integration_range(q[i,non_zero_index], q_min,q_max):
                    S_int_record[n] += Sv[i]
            output_filename = self._generate_output_filename(n)
            self._write_results(output_filename, nu, S_int_record[n])

    def _generate_filename(self, index):
        return f'dynamic_structure_{self.parameters["compute_type"]}_{index}.txt'

    def _generate_output_filename(self, index):
        return f'Interation_dynamic_structure_{self.parameters["compute_type"]}_{index}.txt'

    def _write_results(self, filename, nu, S_int_record):
        intgdyn_write = WriteFile(self.folder, self.parameters)
        intgdyn_write.write_results(filename, nu, S_int_record)

    def _is_within_integration_range(self, q, q_min, q_max):
        return (q_min < q < q_max)

    def compute_q_vectors(self):
        q = np.zeros((self.parameters['vectors'], 3))
        non_zero_index = np.where(self.parameters['q_dir'] != 0)[0][0]
        for i in range(self.parameters['vectors']):
            q[i] = np.array(self.parameters['q_dir']) * (i + 1) / self.parameters['uCell'][
                non_zero_index] * 2 * np.pi / \
                   self.lattice[non_zero_index]
        return q, non_zero_index
