import numpy as np
import random
from Post_Position import Post_Position
class ComputeDynamicFacilitation:
    def __init__(self,atom_positions,atom_velocity, parameters ,lattice):
        self.parameters = parameters
        self.lattice = lattice
        self.atom_positions = atom_positions
        self.atom_velocity = atom_velocity
        self.displacement_calculator = Post_Position()

    def find_smallest_distance(self, atoms_pos: np.ndarray, sorted_drSqr2_index: np.ndarray,
                               sorted_drSqr1_index: np.ndarray, sorted_drSqr1: np.ndarray,
                               sorted_drSqr2: np.ndarray) -> tuple[list[float], list[float]]:
        """
        Find the smallest distance between mobile and non-mobile atoms.
        """
        num_h = int(len(sorted_drSqr2_index) * 0.05)
        mobile_in_t23 = self._get_mobile_atoms(sorted_drSqr2, sorted_drSqr2_index)
        mobile_in_t12 = self._get_mobile_atoms(sorted_drSqr1, sorted_drSqr1_index)

        # Get mobile atoms that only exist in t2-t3 and not in t1-t2
        mobile_in_t23_only = np.setdiff1d(mobile_in_t23, mobile_in_t12)

        non_mobile_in_t12 = sorted_drSqr1_index[num_h:]
        sz_atoms = len(mobile_in_t23_only)
        mobile_in_t23_random = random.choices(non_mobile_in_t12, k=sz_atoms)

        smallest_distance, random_distance_record = [], []
        for i in range(sz_atoms):
            min_dist, min_dist_random = self._calculate_distances(
                np.array(atoms_pos)[mobile_in_t23_only[i]],
                np.array(atoms_pos)[mobile_in_t12],
                np.array(atoms_pos)[mobile_in_t23_random[i]]
            )
            smallest_distance.append(min_dist)
            random_distance_record.append(min_dist_random)

        return smallest_distance, random_distance_record


    def _calculate_distances(self, atom1_pos: np.ndarray, atoms_group1: np.ndarray, atom2_pos: np.ndarray) -> tuple[
        float, float]:
        """
        Calculate the smallest distance for both the actual and random sets of mobile atoms.
        """
        distances = self.displacement_calculator.cal_displacement(atom1_pos, atoms_group1, self.lattice)
        distances_random = self.displacement_calculator.cal_displacement(atom2_pos, atoms_group1, self.lattice)

        min_distance = np.sqrt(np.min(np.sum(distances ** 2, axis=1)))
        min_distance_random = np.sqrt(np.min(np.sum(np.square(distances_random), axis=1)))

        return min_distance, min_distance_random

    def _get_mobile_atoms(self, sorted_drSqr: np.ndarray, sorted_drSqr_index: np.ndarray) -> np.ndarray:
        """
        Identify mobile atoms based on the top 5% of displacements.
        """
        num_h = int(len(sorted_drSqr_index) * 0.05)
        non_zero_mask = sorted_drSqr > 0  # Filter out zero displacements
        filtered_drSqr_index = sorted_drSqr_index[non_zero_mask]
        num_h = min(num_h, len(filtered_drSqr_index))
        return filtered_drSqr_index[:num_h]

    def collect_distance(self) -> tuple[list[list[float]], list[list[float]]]:

        gap_array = range(self.parameters['start_gap'],self.parameters['end_gap'],self.parameters['gap'])
        num_frames = self.parameters['num_frame']
        random.seed(self.parameters['seed'])

        smallest_distance, random_choice_distance = [[] for _ in gap_array], [[] for _ in gap_array]

        for j, gap in enumerate(gap_array):
            for i in range(gap, num_frames - gap):
                sorted_drSqr1, sorted_drSqr1_index = self._calculate_displacements(i - gap, i)
                sorted_drSqr2, sorted_drSqr2_index = self._calculate_displacements(i, i + gap)

                a, b = self.find_smallest_distance(self.atom_positions[i + gap], sorted_drSqr2_index,
                                                   sorted_drSqr1_index, sorted_drSqr1, sorted_drSqr2)
                smallest_distance[j].append(a)
                random_choice_distance[j].append(b)

        return smallest_distance, random_choice_distance

    def compute_probability_distribution(self)-> tuple[np.ndarray, list[np.ndarray], list[np.ndarray]]:
        """
                Compute the probability distribution of distances between mobile atoms over time.
        """
        smallest_distance, random_choice_distance = self.collect_distance()
        rDel = self.parameters['rDel']
        max_shells = int(87 / rDel)
        mean_shells_record = []
        mean_random_shells_record = []
        gap_array = range(self.parameters['start_gap'], self.parameters['end_gap'], self.parameters['gap'])

        sz_gap = len(gap_array)
        for k in range(sz_gap):
            num_samples = len(smallest_distance[k])
            shell_counts = np.zeros([num_samples, max_shells])
            random_shell_counts = np.zeros([num_samples, max_shells])
            for i in range(num_samples):
                for j in range(len(smallest_distance[k][i])):
                    shell_num = int(np.floor(smallest_distance[k][i][j] / rDel))
                    random_shell_num = int(np.floor(random_choice_distance[k][i][j] / rDel))
                    if shell_num > 0:
                        shell_counts[i, shell_num] += 1
                    if random_shell_num > 0:
                        random_shell_counts[i, random_shell_num] += 1

            mean_shells = shell_counts.mean(0)
            mean_random_shells = random_shell_counts.mean(0)
            mean_shells_record.append(mean_shells)
            mean_random_shells_record.append(mean_random_shells)

        r = rDel * (np.array(range(max_shells)) + 0.5)
        return r, mean_shells_record, mean_random_shells_record

    def _calculate_displacements(self, index1: int, index2: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate and sort the displacement vectors and their corresponding indices.
        """
        dr = self.displacement_calculator.cal_displacement(self.atom_positions[index2], self.atom_positions[index1],
                                                          self.lattice)
        drSqr = np.sqrt(np.sum(dr ** 2, axis=1))
        sorted_drSqr = np.sort(drSqr)[::-1]
        sorted_drSqr_index = np.argsort(drSqr)[::-1]
        return sorted_drSqr, sorted_drSqr_index

