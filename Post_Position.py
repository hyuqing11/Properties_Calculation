import numpy as np

class Post_Position:
    def __init__(self):
        pass

    def cal_displacement(self,data1: np.ndarray, data2: np.ndarray, lattice: np.ndarray) -> np.ndarray:
        """Calculate the displacement vector considering periodic boundary conditions."""
        dr = np.array(data1) - np.array(data2)
        dr = dr - np.round(dr / lattice) * lattice
        return dr

    def wrap(self,dr: float, region: float) -> float:
        """Apply perodic boundary condistion."""
        if dr > 0.5 * region:
            dr -= region
        elif dr < -0.5 * region:
            dr += region
        return dr


