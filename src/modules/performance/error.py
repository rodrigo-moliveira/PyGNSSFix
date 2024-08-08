import numpy as np
import scipy.stats as stats

from src import constants
from src.data_mng.csv.csv_data import CSVData
from src.models.frames import cartesian2geodetic, latlon2dcm_e_enu, latlon2dcm_e_ned


def error_vector(row, true_pos, rot_matrix=np.eye(3)):
    pos = np.array([row['x'], row['y'], row['z']])
    return rot_matrix @ (pos - true_pos)

def get_cov_local(row, rot_matrix=np.eye(3)):
    cov_ecef = np.array([[row['cov_xx'], row['cov_xy'], row['cov_xz']],
                         [row['cov_xy'], row['cov_yy'], row['cov_yz']],
                         [row['cov_xz'], row['cov_yz'], row['cov_zz']]])
    cov_local = rot_matrix @ cov_ecef @ rot_matrix.T
    return cov_local


def compute_latlon(position):
    pos_array = position.to_data_array()
    latlon = []
    for x in pos_array[:,0:3]:
        _lla = cartesian2geodetic(*x)
        latlon.append([_lla[0]*constants.RAD2DEG, _lla[1]*constants.RAD2DEG])
    return latlon

def compute_error_static(est_data: CSVData, true_data, true_pos, local="ENU"):
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
    cov_ecef = est_data.data.apply(lambda row: get_cov_local(row), axis=1)

    # Rotate covariance matrix from ECEF to local frame
    cov_local = est_data.data.apply(lambda row: get_cov_local(row, rot_e_local), axis=1)

    return error_matrix_ecef, error_matrix_local, list(cov_ecef), list(cov_local)


def compute_rms_error(error_matrix):
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

def chi_squared_test(residuals, p):
    """
    A chi-square test is a common method to check if the residuals from a Least Squares adjustment in GNSS positioning
    follow a white noise distribution.

    Formulate Hypotheses:
        Null Hypothesis (H0): The residuals follow a white noise distribution.
        Alternative Hypothesis (H1): The residuals do not follow a white noise distribution.

    Steps:
        Obtain the residuals (difference between observed and predicted measurements)
        Calculate the Variance of Residuals. Compute the variance of your residuals. This will be used to normalize the residuals.
        Standardize the Residuals: Convert residuals into standard normal form by dividing them by the variance
        Compute the Chi-Square Statistic
        Degrees of Freedom: The degrees of freedom (df) for the chi-square distribution is equal to the number of
        residuals n minus the number of estimated parameters p in the model.
        Compare with Critical Value: Determine the critical value from the chi-square distribution table at the desired
        significance level (e.g., 0.05) with n−p degrees of freedom.

        Decision Rule:
            * If the computed chi-square statistic is less than the critical value, you fail to reject the null
                hypothesis (the residuals are consistent with white noise).
            * If the computed chi-square statistic is greater than the critical value, you reject the null hypothesis
                (the residuals are not consistent with white noise).

    TODO: check two cases: one when the sigma is the variance of the residuals
        2) when the sigma is the covariance of the measurements.
    :return:
    """

    # Example residuals from Least Squares adjustment
    n = len(residuals)  # number of residuals
    p = 6  # number of estimated parameters, typically 4 for GNSS (x, y, z, clock bias)

    # Calculate the variance of residuals
    variance_residuals = np.var(residuals)
    std_residuals = np.sqrt(variance_residuals)

    # Calculate the mean of residuals
    mean_residuals = np.mean(residuals)

    # Standardize the residuals
    standardized_residuals = (residuals - mean_residuals) / std_residuals
    std_standard = np.std(standardized_residuals)
    print("\tmean", mean_residuals, "std", std_residuals, std_standard)
    # TODO check now that np.std() is 1

    # Degrees of freedom
    degrees_of_freedom = n - p

    # Compute the chi-square statistic
    chi_square_statistic = np.sum(standardized_residuals ** 2)

    # Determine the critical value from chi-square distribution at 0.05 significance level
    alpha = 0.05
    critical_value = stats.chi2.ppf(1 - alpha, degrees_of_freedom)

    # Decision rule
    if chi_square_statistic < critical_value:
        print("Fail to reject the null hypothesis: residuals are consistent with white noise.")
    else:
        print("Reject the null hypothesis: residuals are not consistent with white noise.")

    print(f"\tChi-square statistic: {chi_square_statistic}")
    print(f"\tCritical value (0.05 significance level): {critical_value}\n")
    #print(f"Degrees of freedom: {degrees_of_freedom}")


def shapiro_test(residuals):
    # Perform the Shapiro-Wilk test
    statistic, p_value = stats.shapiro(residuals)

    # Set the significance level (alpha)
    alpha = 0.05

    # Print the results
    print(f"Shapiro-Wilk Test Statistic: {statistic}")
    print(f"P-value: {p_value}")

    # Check the p-value against the significance level
    if p_value > alpha:
        print("The residuals follow a normal distribution (fail to reject H0)")
    else:
        print("The residuals do not follow a normal distribution (reject H0)")