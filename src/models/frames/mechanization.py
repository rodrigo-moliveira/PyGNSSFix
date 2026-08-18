# Utility functions for the mechanization equations of the INS
# TODO: update docstrings

import numpy as np

from src.constants import EARTH_ROTATION
from .frames import get_earth_radii
from src.utils.math_utils import vector2skew_symmetric

w_ie_e = np.array([0, 0, EARTH_ROTATION])  # in rad/s. This vector is also w_ie_i, although this is not needed


# Methods to compute angular velocities between frames e - n - b
def compute_w_en_n(v_eb_n, lld):
    v_n, v_e, v_z = v_eb_n[:]
    lat, long, d = lld[:]

    rm, rn = get_earth_radii(lat)

    rm_effective = rm - d
    rn_effective = rn - d

    w_en_n = np.zeros(3)

    w_en_n[0] = v_e / rn_effective  # wN
    w_en_n[1] = -v_n / rm_effective  # wE
    w_en_n[2] = -v_e * np.tan(lat) / rn_effective  # wD

    return w_en_n, rm_effective, rn_effective


def compute_w_ie_e():
    return w_ie_e


def compute_w_ie_n(lat):
    # lat in rad
    return np.array([np.cos(lat),
                     0,
                     -np.sin(lat)]) * EARTH_ROTATION


def compute_w_nb_b(euler_dot, euler):
    # euler angles (phi, theta, psi) = (roll, pitch, yaw)

    phi_dot, theta_dot, psi_dot = euler_dot[:]
    phi, theta, psi = euler[:]

    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    cos_theta = np.cos(theta)

    w_nb_b = np.zeros(3)

    w_nb_b[0] = phi_dot - psi_dot * np.sin(theta)
    w_nb_b[1] = theta_dot * cos_phi + psi_dot * cos_theta * sin_phi
    w_nb_b[2] = psi_dot * cos_phi * cos_theta - theta_dot * sin_phi

    return w_nb_b


###################################
# N-Frame Mechanization Equations #
###################################

def compute_w_ib_b(c_nb, w_ie_n, w_en_n, w_nb_b):
    # Gyro readout
    return c_nb @ w_ie_n + c_nb @ w_en_n + w_nb_b


def compute_f_ib_b(c_nb, c_en, w_en_n, w_ie_n, v_eb_n_dot, g_eb_e, v_eb_n):
    # accelerometer readout
    skew_en_n = vector2skew_symmetric(w_en_n)
    skew_ie_n = vector2skew_symmetric(w_ie_n)
    f_ib_n = v_eb_n_dot - c_en @ g_eb_e + (skew_en_n + 2 * skew_ie_n) @ v_eb_n
    f_ib_b = c_nb @ f_ib_n  # accelerometer readout
    return f_ib_b
