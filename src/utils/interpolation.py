from scipy.interpolate import interp1d
from src.utils.deprecated import deprecated
import numpy as np

# TODO: review all docstrings and add newton divided method into the configuration

# @deprecated("use function :py:function:`src.utils.interpolation.linear_interpolation` instead")
def linear_interpolation_scipy(x, x_vals, y_vals):
    """
    Perform linear interpolation using scipy's interp1d function.

    Parameters:
        x (float or array-like): The x-value(s) at which to estimate y.
        x_vals (array-like): Array of known x-values.
        y_vals (array-like): Array of corresponding y-values for x_vals.

    Returns:
        float or array-like: The estimated y-value(s) at x using linear interpolation.
    """
    # Create a linear interpolation function
    interp_func = interp1d(x_vals, y_vals, kind='linear', fill_value='extrapolate')

    # Use the interpolation function to estimate y at the given x value(s)
    y = interp_func(x)

    return y


def linear_interpolation(x, x0, y0, x1, y1):
    """
    Perform linear interpolation to estimate y-value at given x.

    Parameters:
        x (float): The x-value at which to estimate y.
        x0 (float): x-value of the first known point.
        y0 (float): y-value of the first known point corresponding to x0.
        x1 (float): x-value of the second known point.
        y1 (float): y-value of the second known point corresponding to x1.

    Returns:
        float: The estimated y-value at x using linear interpolation.
    """
    # Check if the two x-values are the same (division by zero prevention)
    if x0 == x1:
        raise ValueError("x0 and x1 cannot be the same for interpolation")

    # Linear interpolation formula
    y = y0 + ((y1 - y0) / (x1 - x0)) * (x - x0)

    return y


def lagrange_interpolation(x_values, y_values, degree):
    """
    Perform Lagrange interpolation of specified degree.

    Parameters:
        x_values (list or array-like): List of x-values (data points).
        y_values (list or array-like): List of corresponding y-values.
        degree (int): Degree of the Lagrange interpolation polynomial.

    Returns:
        function: A function representing the Lagrange interpolation polynomial.
    """
    # Ensure the number of data points matches the degree + 1
    if len(x_values) != degree + 1 or len(y_values) != degree + 1:
        raise ValueError(f"Number of data points (len={len(x_values)}) should be equal to degree + 1 (degree={degree})")

    def lagrange_poly(x):
        result = 0.0
        for i in range(degree + 1):
            term = y_values[i]
            for j in range(degree + 1):
                if i != j:
                    term *= (x - x_values[j]) / (x_values[i] - x_values[j])
            result += term
        return result

    return lagrange_poly


def lagrange_interpolation_derivative(x_values, y_values, degree):
    """
    Computes the derivative of the Lagrange polynomial function of the specified degree.
    For example, it is used to compute velocity vectors using position data points

    Parameters:
        x_values (list or array-like): List of x-values (data points).
        y_values (list or array-like): List of corresponding y-values.
        degree (int): Degree of the Lagrange interpolation polynomial.

    Returns:
        function: A function representing the derivative of the Lagrange interpolation polynomial.
    """
    # Ensure the number of data points matches the degree + 1
    if len(x_values) != degree + 1 or len(y_values) != degree + 1:
        raise ValueError(f"Number of data points (len={len(x_values)}) should be equal to degree + 1 (degree={degree})")

    def lagrange_poly(x):
        result = 0.0
        for i in range(degree + 1):
            cum_sum = 0
            for j in range(degree + 1):
                if j != i:
                    cum_prod = 1 / (x_values[i] - x_values[j])
                    for m in range(degree + 1):
                        if m != i and m != j:
                            cum_prod *= (x - x_values[m]) / (x_values[i] - x_values[m])
                    cum_sum += cum_prod

            result += y_values[i] * cum_sum
        return result

    return lagrange_poly


# I had some trouble setting up the derivative of the lagrange interpolation, this function also seems to work,
# but it is untested
@deprecated("use function :py:function:`src.utils.interpolation.lagrange_interpolation_derivative` instead")
def lagrange_interpolation_derivative_alternative(x_values, y_values, degree):
    """
    Computes the derivative of the Lagrange polynomial function of the specified degree.
    For example, it is used to compute velocity vectors using position data points

    Parameters:
        x_values (list or array-like): List of x-values (data points).
        y_values (list or array-like): List of corresponding y-values.
        degree (int): Degree of the Lagrange interpolation polynomial.

    Returns:
        function: A function representing the derivative of the Lagrange interpolation polynomial.
    """
    # Ensure the number of data points matches
    if len(x_values) != len(y_values):
        raise ValueError("Number of x-values must be equal to the number of y-values.")

    # Compute the derivative of the Lagrange interpolation polynomial (velocity)
    def velocity_function(t):
        velocity = 0.0
        for j in range(len(x_values)):
            lagrange_derivative = 0.0
            for m in range(len(x_values)):
                if m != j:
                    derivative_term = 1.0 / (x_values[j] - x_values[m])
                    for n in range(len(x_values)):
                        if n != j and n != m:
                            derivative_term *= (t - x_values[n]) / (x_values[j] - x_values[n])
                    lagrange_derivative += derivative_term
            velocity += y_values[j] * lagrange_derivative
        return velocity

    return velocity_function


def newton_divided_differences(x_values, y_values, degree):
    """
    Perform Newton's Divided Differences interpolation for vector-valued functions.

    Parameters:
        x_values (numpy.ndarray): Vector of x-values (data points).
        y_values (numpy.ndarray): Matrix of corresponding y-values, where each row is a vector [y0, y1, ..., yn].
        degree (int): Desired degree of the Newton interpolation polynomial.

    Returns:
        function: A function representing the Newton interpolation polynomial.
    """
    # Ensure x_values is a column vector
    x_values = np.atleast_2d(x_values).T if x_values.ndim == 1 else x_values
    if x_values.shape[1] != 1:
        raise ValueError('x_values must be a column vector')

    # Ensure the number of data points in y_values matches the number of x_values
    if y_values.shape[0] != x_values.shape[0]:
        raise ValueError('Number of rows in y_values must match the length of x_values')

    # Ensure the number of data points is sufficient for the desired degree
    if len(x_values) < degree + 1:
        raise ValueError('Not enough data points for the desired degree of interpolation')

    # Use only the first degree + 1 points
    x_values = x_values[:degree + 1]
    y_values = y_values[:degree + 1, :]

    # Number of data points and dimensions
    n = len(x_values)
    num_dimensions = y_values.shape[1]

    # Initialize divided differences table
    divided_diffs = np.zeros((n, n, num_dimensions))
    divided_diffs[:, 0, :] = y_values

    # Compute the divided differences table
    for j in range(1, n):
        for i in range(n - j):
            for dim in range(num_dimensions):
                divided_diffs[i, j, dim] = (
                                                   divided_diffs[i + 1, j - 1, dim] - divided_diffs[i, j - 1, dim]
                                           ) / (x_values[i + j] - x_values[i])

    def newton_poly_func(x):
        result = np.zeros(num_dimensions)
        for dim in range(num_dimensions):
            term = divided_diffs[0, 0, dim]
            for j in range(1, n):
                term_j = divided_diffs[0, j, dim]
                for k in range(j):
                    term_j *= (x - x_values[k])
                term += term_j
            result[dim] = term
        return result

    return newton_poly_func


# Example Usage
if __name__ == "__main__":
    x_values_ = np.array([0, 1, 2, 3]).reshape(-1, 1)
    y_values_ = np.array([[1, 2], [2, 3], [3, 4], [4, 5]])
    degree = 2

    newton_poly = newton_divided_differences(x_values_, y_values_, degree)
    x_test = np.array([1.5])
    y_test = newton_poly(x_test)
    print(y_test)
