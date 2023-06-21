import numpy as np

from math import sin, cos

from PositioningSolver.src.data_types.basics.DataType import DataTypeFactory
from PositioningSolver.src.math_utils.Constants import Constant

C1 = DataTypeFactory("C1")
f1 = C1.freq  # L1 frequency for GPS


def ionosphereCorrection(user_lat, user_long, sv_el, sv_az, alfa, beta, GPS_time, frequency):
    """
    This function computes the ionosphere correction of pseudorange measurements for online single frequency users,
    using the a priori Klobuchar Ionospheric Model
    refs:
        * section 20.3.3.5.2.5 of **REF[3]**
        * https://gssc.esa.int/navipedia/index.php/Klobuchar_Ionospheric_Model

    GPS satellites broadcast the parameters of the Klobuchar ionospheric model for single frequency users.
    The Klobuchar model was designed to minimise user computational complexity and user computer storage as far
    as to keep a minimum number of coefficients to transmit on satellite-user link.
    The model is estimated to reduce about the 50% RMS ionospheric range error worldwide

    Args:
        user_lat (float) : user latitude [rad]
        user_long (float) : user longitude [rad]
        sv_el (float) : satellite elevation [rad]
        sv_az (float) : satellite azimuth [rad]
        alfa (list) : list of alfa parameters, length 4
        beta (list) : list of beta parameters, length 4
        GPS_time (src.data_types.basics.Epoch.Epoch) : GPS epoch to compute the iono
        frequency (src.data_types.basics.DataType.DataType) : user frequency (the iono is frequency dependent)


    Ionosphere delay is frequency dependent. The algorithm gives the delay with respect to the L1 frequency. For L2 / L5
    users, a correction is therefore needed.

    Return:
        float : ionosphere correction [m]

    """

    # Get angles in semicircles
    sv_el_semi = sv_el / Constant.PI
    user_lat_semi = user_lat / Constant.PI
    user_long_semi = user_long / Constant.PI

    # Calculate the earth-centred angle (elevation in semicircles)
    psi = 0.0137 / (sv_el_semi + 0.11) - 0.022

    # Compute and fix the latitude of the Ionospheric Pierce Point(IPP)
    lat_IPP = user_lat_semi + psi * cos(sv_az)
    if lat_IPP > 0.416:
        lat_IPP = 0.416
    elif lat_IPP < -0.416:
        lat_IPP = -0.416

    # Compute the longitude of the IPP
    long_IPP = user_long_semi + (psi * sin(sv_az)) / cos(lat_IPP * Constant.PI)

    # Find the geomagnetic latitude of the IPP
    lat_m = lat_IPP + 0.064 * cos((long_IPP - 1.617) * Constant.PI)

    # Find the local time at the IPP
    t = Constant.SECONDS_IN_DAY / 2 * long_IPP + GPS_time.seconds
    t = t % Constant.SECONDS_IN_DAY

    # Compute the amplitude of ionospheric delay.
    A_I = alfa[0] + alfa[1] * lat_m + alfa[2] * (lat_m ** 2) + alfa[3] * (lat_m ** 3)
    if A_I < 0:
        A_I = 0

    # Compute the period of ionospheric delay
    P_I = beta[0] + beta[1] * lat_m + beta[2] * (lat_m ** 2) + beta[3] * (lat_m ** 3)
    if P_I < 72000:
        P_I = 72000

    # Compute the phase of ionospheric delay
    X_I = 2 * Constant.PI * (t - 50400) / P_I

    # Compute the slant factor (elevation in semicircles).
    F = 1.0 + 16.0 * (0.53 - sv_el_semi) ** 3

    # Compute the ionospheric time delay [s]
    if abs(X_I) > 1.57:
        iono = 5E-9 * F
    else:
        iono = (5E-9 + A_I * (1 - (X_I ** 2) / 2 + (X_I ** 4) / 24)) * F

    # fix I for non L1 users
    if frequency != f1:
        iono = (f1.freq_value / frequency.freq_value) ** 2 * iono

    # get ionosphere in meters
    iono = iono * Constant.SPEED_OF_LIGHT
    return iono


# Constants for the tropospheric Saastamoinen model
LimLat = (15, 30, 45, 60, 75)

#                   P0(mbar) T0(K)    e0(mbar) beta(K/m)  lambda0
P_mean = np.array([[1013.25, 299.65,  26.31,   6.30e-3,   2.77],
                   [1017.25, 294.15,  21.79,   6.05e-3,   3.15],
                   [1015.75, 283.15,  11.66,   5.58e-3,   2.57],
                   [1011.75, 272.15,  6.78,    5.39e-3,   1.81],
                   [1013.00, 263.65,  4.11,    4.53e-3,   1.55]])

#                     DP(mbar)   DT(K)   De(mbar) Db(K/m)   dl
P_season = np.array([[0.00,      0.00,   0.00,    0.00e-3,  0.00],
                     [-3.75,     7.00,   8.85,    0.25e-3,  0.33],
                     [-2.25,     11.00,  7.24,    0.32e-3,  0.46],
                     [-1.75,     15.00,  5.36,    0.81e-3,  0.74],
                     [-0.50,     14.50,  3.39,    0.62e-3,  0.30]])


def troposphericCorrection(h, lat, DOY, el):
    """
    This function computes the tropospheric correction of pseudorange measurements for online users,
    using the a priori Saastamoinen Model

    This code has been adapted from:
        Open source PANG-NAV software (https://geodesy.noaa.gov/gps-toolbox/PANG-NAV.htm), tropo_correction0.m script,
        implementing the Saastamoinen Model


    Args:
        h (float) : user altitude (geodetic coordinate)     [m]
        lat (float) : user latitude (geodetic coordinate)   [rad]
        DOY (float) : Day of the year (from 1 to 365)       [1 - 365]
        el (float) : user elevation                         [rad]

    Return:
        float : tropospheric correction [m]


    """

    # convert lat to degrees
    lat = Constant.RAD2DEG * lat

    if lat < 0:
        D_star = 211
    else:
        D_star = 28

    if lat > LimLat[-1]:  # if lat > 75 [deg]
        P0 = P_mean[-1]
        DP = P_season[-1]

    elif lat < LimLat[0]:  # if lat < 15 [deg]
        P0 = P_mean[0]
        DP = P_season[0]

    else:
        _iis = [lat_i < lat for lat_i in LimLat]
        _iis.reverse()
        _i = 4 - _iis.index(True)
        m = (lat - LimLat[_i]) / (LimLat[_i + 1] - LimLat[_i])
        P0 = P_mean[_i] + (P_mean[_i + 1] - P_mean[_i]) * m
        DP = P_season[_i] + (P_season[_i + 1] - P_season[_i]) * m

    Par = P0 - DP * cos(2 * np.pi * (DOY - D_star) / 365.25)

    P, T, e, beta, Lambda = Par[:]

    fs = 1 - 0.00266 * cos(2 * lat) - 0.00000028 * h

    D_z_dry = (0.0022768 - 0.0000005) * P / fs

    D_z_wet = (0.002277 * (1255 / T + 0.05) * e) / fs

    M = 1.001 / np.sqrt(0.002001 + (sin(el)) ** 2)

    dD_tropo = (D_z_dry + D_z_wet) * M

    return dD_tropo
