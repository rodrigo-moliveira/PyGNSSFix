import numpy as np
from scipy import stats

def white_noise_test(fh, data, sensor):
    statistics, pvals = np.apply_along_axis(stats.shapiro, 0, np.array(data))

    for statistic, p_value, x_label in zip(statistics, pvals, ["X", "Y", "Z"]):

        # Set the significance level (alpha)
        alpha = 0.05

        # Print the results
        fh.write(f"[{sensor}] {x_label}-axis Shapiro-Wilk Test Statistic: {statistic}\n")
        fh.write(f"[{sensor}] {x_label}-axis P-value: {p_value}\n")

        # Check the p-value against the significance level
        if p_value > alpha:
            fh.write(f"[{sensor}] {x_label}-axis The data follows a normal distribution (fail to reject H0)\n")
        else:
            fh.write(f"[{sensor}] {x_label}-axis The data does not follow a normal distribution (reject H0)\n")
def Chi_Square_Residual_Test():
    # Example pseudorange observations and other required data
    P_r = np.array([...])  # observed pseudoranges
    rho_r = np.array([...])  # predicted pseudoranges
    c = 299792458  # speed of light in m/s
    dt_r = ...  # estimated receiver clock bias
    dT = np.array([...])  # satellite clock biases
    I_r = np.array([...])  # ionospheric delays
    T_r = np.array([...])  # tropospheric delays

    # Number of measurements and estimated parameters
    m = len(P_r)
    n = 4  # typically 4 for GNSS (x, y, z, clock bias)

    # Calculate residuals
    residuals = P_r - (rho_r + c * dt_r - c * dT + I_r + T_r)

    # Compute the variance and standard deviation of the residuals
    variance_residuals = np.var(residuals, ddof=1)
    std_residuals = np.sqrt(variance_residuals)

    # Standardize the residuals
    standardized_residuals = residuals / std_residuals

    # Compute the test statistic
    test_statistic = np.sum(standardized_residuals**2) / (m - n - 1)

    # Determine the critical value from chi-square distribution at alpha = 0.001
    alpha = 0.001
    degrees_of_freedom = m - n - 1
    critical_value = stats.chi2.ppf(1 - alpha, degrees_of_freedom)

    # Decision rule
    if test_statistic < critical_value:
        print("Fail to reject the null hypothesis: residuals are consistent with white noise.")
    else:
        print("Reject the null hypothesis: residuals are not consistent with white noise.")

    print(f"Test statistic: {test_statistic}")
    print(f"Critical value (alpha = 0.001): {critical_value}")
    print(f"Degrees of freedom: {degrees_of_freedom}")
