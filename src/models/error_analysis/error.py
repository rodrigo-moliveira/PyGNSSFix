""" Module with utility functions for performing estimation error and residual analysis"""

import numpy as np
import scipy.stats as stats

from src.models.frames import cartesian2geodetic, latlon2dcm_e_enu, latlon2dcm_e_ned


def compute_error_static(est_data, true_data, true_pos, local="ENU"):
    """
    Computes estimation errors and covariances.

    Args:
        est_data (src.data_mng.csv.csv_data.CSVData): The estimated data object (position or velocity) in CSV format.
        true_data (numpy.ndarray or list): The true/reference vector to compare against the estimated data.
        true_pos (numpy.ndarray or list): The true/reference position used to compute the local transformation matrix.
        local (str): The local frame to use for computations. Options are 'ENU' (East-North-Up) or
                           'NED' (North-East-Down). Defaults to 'ENU'.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray, list[numpy.ndarray], list[numpy.ndarray]]: The computed outputs are:

            - error_ecef (numpy.ndarray): An Nx3 array representing the estimation error in the ECEF frame for each
              epoch (true value minus estimated value).
            - error_local (numpy.ndarray): An Nx3 array representing the estimation error in the local frame for each
              epoch (true value minus estimated value).
            - cov_ecef (list[numpy.ndarray]): A list of N 3x3 covariance matrices representing the estimation
              uncertainty in the ECEF frame for each epoch.
            - cov_local (list[numpy.ndarray]): A list of N 3x3 covariance matrices representing the estimation
              uncertainty in the local frame for each epoch.

        where N corresponds to the number of epochs.

    Raises:
        ValueError: an exception is raised if the `local` argument is not valid
    """
    [lat, lon, _] = cartesian2geodetic(*true_pos)

    if local == "ENU":
        rot_e_local = latlon2dcm_e_enu(lat, lon)  # rotation matrix from ECEF to ENU
    elif local == "NED":
        rot_e_local = latlon2dcm_e_ned(lat, lon)  # rotation matrix from ECEF to NED
    else:
        raise ValueError(f"Argument `local` must either be 'ENU' or 'NED' (local={local})")

    # Compute errors without rotation matrix (default identity matrix)
    errors_ecef = est_data.data.apply(lambda row: error_vector(row, true_data), axis=1)
    error_matrix_ecef = np.vstack(errors_ecef.values)

    # Compute errors with the specified rotation matrix
    errors_local = est_data.data.apply(lambda row: error_vector(row, true_data, rot_e_local), axis=1)
    error_matrix_local = np.vstack(errors_local.values)

    # Get ECEF Covariance per epoch
    cov_ecef = est_data.data.apply(lambda row: get_cov(row), axis=1)

    # Rotate covariance matrix from ECEF to local frame
    cov_local = est_data.data.apply(lambda row: get_cov(row, rot_e_local), axis=1)

    return error_matrix_ecef, error_matrix_local, list(cov_ecef), list(cov_local)


def compute_rms_error(error_matrix):
    """
    Computes the Root Mean Square (RMS) Error for the provided error matrix time series.
    The RMS error is calculated for each component (x, y, z), as well as for the 2D (x, y) and 3D (x, y, z)
    vectors.

    The RMS error for a component `c` is computed as:
        RMS_c = sqrt(1/N * sum((error_c[i])^2 for i in range(N)))
    where `N` is the number of estimation epochs, and `error_c[i]` is the error for the `i`-th epoch.

    Args:
        error_matrix (numpy.ndarray): An Nx3 array representing the estimation error for each epoch, where N is
                                      the number of estimation epochs. Each row corresponds to the error in the
                                      x, y, and z components for a single epoch.
    Returns:
        dict: A dictionary containing the computed RMS errors with the following keys:

            - "x" (float): RMS error for the x component.
            - "y" (float): RMS error for the y component.
            - "z" (float): RMS error for the z component.
            - "2D" (float): RMS error for the 2D vector (x, y).
            - "3D" (float): RMS error for the 3D vector (x, y, z).
    """

    accum_x = accum_y = accum_z = accum_2d = accum_3d = 0

    for error in error_matrix:
        accum_x += error[0] * error[0]
        accum_y += error[1] * error[1]
        accum_z += error[2] * error[2]

        accum_2d += error[0] * error[0] + error[1] * error[1]
        accum_3d += error[0] * error[0] + error[1] * error[1] + error[2] * error[2]

    # 1D RMS stats
    rms_x = np.sqrt(1 / len(error_matrix) * accum_x)
    rms_y = np.sqrt(1 / len(error_matrix) * accum_y)
    rms_z = np.sqrt(1 / len(error_matrix) * accum_z)

    # 2D RMS stats
    rms_2d = np.sqrt(1 / len(error_matrix) * accum_2d)

    # 3D RMS stats
    rms_3d = np.sqrt(1 / len(error_matrix) * accum_3d)

    return {"x": rms_x,
            "y": rms_y,
            "z": rms_z,
            "2D": rms_2d,
            "3D": rms_3d
            }


def error_vector(row, true_pos, rot_matrix=np.eye(3)):
    """
    Computes the estimation error vector in either the ECEF frame or a local reference frame.

    This function calculates the error vector by subtracting the true position from the estimated position
    provided in a pandas row. The resulting error vector is optionally transformed into a local reference
    frame using the provided rotation matrix. If the identity matrix is used as the rotation matrix, the
    error vector remains in the ECEF (Earth-Centered, Earth-Fixed) frame.

    Args:
        row (pandas.Series): A pandas Series representing a single row of data with the following columns:

            - 'x' (float): The estimated x-coordinate in the ECEF frame.
            - 'y' (float): The estimated y-coordinate in the ECEF frame.
            - 'z' (float): The estimated z-coordinate in the ECEF frame.
        true_pos (numpy.ndarray): A 3-element array-like object representing the true position
                                  [x_true, y_true, z_true] in the ECEF frame.
        rot_matrix (numpy.ndarray): A 3x3 rotation matrix to transform the error vector  into a local reference frame.
            If set to the identity matrix (default), the output error vector remains in the ECEF frame.

    Returns:
        numpy.ndarray:
            A 3-element array representing the estimation error vector. If `rot_matrix` is the identity
            matrix, the error vector is in the ECEF frame; otherwise, it is in the specified local
            reference frame, calculated as:

                error_vector = rot_matrix @ (pos - true_pos)

            where `pos` is the estimated position [x, y, z] from the row, and `@` denotes matrix multiplication.
    """

    pos = np.array([row['x'], row['y'], row['z']])
    return rot_matrix @ (pos - true_pos)


def get_cov(row, rot_matrix=np.eye(3)):
    """
    Computes the covariance matrix in either the ECEF frame or a local reference frame.

    This function constructs the 3x3 covariance matrix from the provided components in a pandas row, where the
    components are given in the ECEF (Earth-Centered, Earth-Fixed) frame. The function optionally transforms
    this covariance matrix into a local reference frame using the provided rotation matrix. If the identity
    matrix is used as the rotation matrix, the covariance matrix remains in the ECEF frame.

    Args:
        row (pandas.Series): A pandas Series representing a single row of data with the following covariance components:
            - 'cov_xx' (float): The xx component of the covariance matrix in the ECEF frame.
            - 'cov_yy' (float): The yy component of the covariance matrix in the ECEF frame.
            - 'cov_zz' (float): The zz component of the covariance matrix in the ECEF frame.
            - 'cov_xy' (float): The xy component of the covariance matrix in the ECEF frame.
            - 'cov_xz' (float): The xz component of the covariance matrix in the ECEF frame.
            - 'cov_yz' (float): The yz component of the covariance matrix in the ECEF frame.
        rot_matrix (numpy.ndarray, optional): A 3x3 rotation matrix to transform the covariance matrix into
                                              a local reference frame. If set to the identity matrix (default),
                                              the output covariance matrix remains in the ECEF frame.

    Returns:
        numpy.ndarray: A 3x3 covariance matrix. If `rot_matrix` is the identity matrix, the covariance matrix
                       is in the ECEF frame; otherwise, it is transformed into the specified local reference frame,
                       calculated as:

                       cov_local = rot_matrix @ cov_ecef @ rot_matrix.T

                       where `cov_ecef` is the covariance matrix in the ECEF frame constructed from the row data.
    """

    cov_ecef = np.array([[row['cov_xx'], row['cov_xy'], row['cov_xz']],
                         [row['cov_xy'], row['cov_yy'], row['cov_yz']],
                         [row['cov_xz'], row['cov_yz'], row['cov_zz']]])
    cov_local = rot_matrix @ cov_ecef @ rot_matrix.T
    return cov_local


def chi_squared_test(residuals, p, alpha):
    """
    This function performs the chi-square test, which is a common method to check if the residuals from
    a Least-Squares adjustment in GNSS positioning follow a white noise distribution.

    Formulate Hypotheses:
        Null Hypothesis (H0): The residuals follow a white noise distribution.
        Alternative Hypothesis (H1): The residuals do not follow a white noise distribution.

    Steps:
        1) Calculate the Variance of Residuals.
        2) Standardize the Residuals: convert residuals into standard normal form
        3) Compute the Chi-Square Statistic
        4) Degrees of Freedom: the degrees of freedom (df) for the chi-square distribution is equal to the number of
            residuals n minus the number of estimated parameters p in the model.
        5) Compare with the Critical Value: Determine the critical value from the chi-square distribution table at the
            desired significance level (e.g., 0.05) with n−p degrees of freedom.

    Decision Rule:
        * If the computed chi-square statistic is less than the critical value, you fail to reject the null
            hypothesis (the residuals are consistent with white noise).
        * If the computed chi-square statistic is greater than the critical value, you reject the null hypothesis
            (the residuals are not consistent with white noise).

    Args:
        residuals (numpy.ndarray): An array-like object of size N representing the residuals to be evaluated.
        p (int): number of estimated parameters affecting the residuals. For example if position(3) , clock(1),
            troposphere (1) and ionosphere (1) are estimated, p should be 6
        alpha(float): significance level for the chi-squared distribution (usually 0.05)
    Returns:
        str: a string report with the test results, to be written to a file
    """
    report = ""
    # Example residuals from Least Squares adjustment
    n = len(residuals)  # number of residuals

    # Calculate the variance of residuals
    variance_residuals = np.var(residuals)
    std_residuals = np.sqrt(variance_residuals)

    # Calculate the mean of residuals
    mean_residuals = np.mean(residuals)

    # Standardize the residuals
    standardized_residuals = (residuals - mean_residuals) / std_residuals

    # Degrees of freedom
    degrees_of_freedom = n - p

    # Compute the chi-square statistic
    chi_square_statistic = np.sum(standardized_residuals ** 2)

    # Determine the critical value from chi-square distribution at 0.05 significance level
    critical_value = stats.chi2.ppf(1 - alpha, degrees_of_freedom)

    # Decision rule
    if chi_square_statistic < critical_value:
        report += "Fail to reject the null hypothesis (chi_square_statistic < critical_value): " \
                  "residuals are consistent with white noise.\n"
    else:
        report += "Reject the null hypothesis (chi_square_statistic >= critical_value): " \
                  "residuals are not consistent with white noise.\n"
    report += f"\t* Chi-square statistic = {chi_square_statistic}\n"
    report += f"\t* Critical value ({alpha} significance level, {degrees_of_freedom} degrees of freedom)" \
              f" = {critical_value}\n"
    return report


def shapiro_test(residuals, alpha):
    """
    Performs the Shapiro-Wilk test for normality on a set of residuals.

    The Shapiro-Wilk test evaluates the null hypothesis that a given sample
    comes from a normally distributed population. It is particularly useful for
    small to moderate sample sizes and is commonly used to assess the
    normality of residuals in regression models.

    Args:
        residuals(numpy.ndarray) : an array-like with a sequence of residuals (or any numerical data) to be
            tested for normality.
        alpha(float): significance level for the chi-squared distribution (usually 0.05)

    Returns:
        str: A string report detailing the outcome of the Shapiro-Wilk test,
            including whether the null hypothesis of normality was rejected or not.
            The report contains:
            - The test statistic.
            - The p-value.
            - The decision based on the specified significance level.
    """
    report = ""

    statistic, p_value = stats.shapiro(residuals)

    # Check the p-value against the significance level
    if p_value > alpha:
        report += "Fail to reject the null hypothesis (p_value > alpha): " \
                  "The residuals follow a normal distribution.\n"
    else:
        report += "Reject the null hypothesis (p_value <= alpha): " \
                  "residuals are not consistent with a normal distribution.\n"

    report += f"\t* Statistic={statistic}\n\t* p_value={p_value}\n\t* significance_level={alpha}\n"
    return report
