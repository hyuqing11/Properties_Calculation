import numpy as np
from MathFunctions import MathFunctions
import pytest
def test_trapezoidal_integration():
    math_functions = MathFunctions()
    f = np.array([1, 2, 3, 4, 5])
    nf = len(f)

    result = math_functions.trapezoidal_integration(f, nf)
    expected_output = 12.0

    assert result == pytest.approx(expected_output, rel=1e-3)


def test_sum_of_squares():
    math_functions = MathFunctions()
    # Define input matrix
    matrix = np.array([[1, 2], [3, 4], [5, 6]])  # Sample matrix

    # Expected output
    expected_output = np.array([np.sqrt(1 ** 2 + 2 ** 2), np.sqrt(3 ** 2 + 4 ** 2), np.sqrt(5 ** 2 + 6 ** 2)])

    # Call the method
    result = math_functions.sum_of_squares(matrix)

    # Compare result with expected output
    np.testing.assert_array_almost_equal(result, expected_output, decimal=3)



