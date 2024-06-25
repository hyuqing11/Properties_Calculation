import numpy as np
import math
class MathFunctions:
    def __init(self,f):
        pass
    def compute_fourier_transform(self,f,Nc, omega, dt):
        #f is the funtion that you would like to do the Fourier transform
        # ff is the fourier transform of f
        #f = np.transpose(f)
        f = f * (np.cos(np.pi * np.arange(Nc) / Nc) + 1) * 0.5  # window function
        f = f * np.append(np.ones(1), 2 * np.ones(Nc - 1))  # C(t) = C(-t)
        ff = np.zeros(len(omega))
        for n in range(len(omega)):  # Discrete cosine transform
            ff[n] = dt * np.sum(f * np.cos(omega[n] * np.arange(Nc) * dt))
        return ff


    def trapezoidal_integration(self,f,nf):
        s = 0.5 * (f[0] + f[nf - 1])
        for i in range(1, nf - 1):
            s = s + f[i]
        return s


    def sum_of_squares(self,matrix):
        return np.sqrt([sum(x**2 for x in row) for row in matrix])


