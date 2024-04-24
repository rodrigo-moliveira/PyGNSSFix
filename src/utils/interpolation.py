from scipy.interpolate import interp1d, lagrange


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
    TODO: need to check this function
        see also https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.lagrange.html
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
