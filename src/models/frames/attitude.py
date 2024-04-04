"""Useful attitude functions and conversions
"""

import numpy as np

from .frames import cartesian2geodetic, latlon2dcm_e_ned
from src.utils.math_utils import vector2skew_symmetric, require_len_array, require_len_matrix


def euler2dcm(euler):
    """
    Convert Euler angles to direction cosine matrix (DCM). The conventional rotation sequence z-y-x is adopted.

    The DCM matrix is the following:
        [          cy*cz,          cy*sz,            -sy]
        [ sy*sx*cz-sz*cx, sy*sx*sz+cz*cx,          cy*sx]
        [ sy*cx*cz+sz*sx, sy*cx*sz-cz*sx,          cy*cx]
    where (x,y,z) are (roll,pitch,yaw), s is short for sin and c is short for cos

    Args:
        euler (np.ndarray or list): 3x1 Euler angles (roll, pitch, yaw), rad.
    Returns:
        np.ndarray: 3x3 DCM coordinate transformation matrix from n to b
    """
    require_len_array(euler)

    dcm = np.zeros((3, 3))
    c_angle = np.cos(euler)
    s_angle = np.sin(euler)

    # dcm = rot1(angles[0]) @ rot2(angles[1]) @ rot3(angles[2])
    dcm[0, 0] = c_angle[1] * c_angle[2]
    dcm[0, 1] = c_angle[1] * s_angle[2]
    dcm[0, 2] = -s_angle[1]
    dcm[1, 0] = s_angle[0] * s_angle[1] * c_angle[2] - c_angle[0] * s_angle[2]
    dcm[1, 1] = s_angle[0] * s_angle[1] * s_angle[2] + c_angle[0] * c_angle[2]
    dcm[1, 2] = c_angle[1] * s_angle[0]
    dcm[2, 0] = s_angle[1] * c_angle[0] * c_angle[2] + s_angle[2] * s_angle[0]
    dcm[2, 1] = s_angle[1] * c_angle[0] * s_angle[2] - c_angle[2] * s_angle[0]
    dcm[2, 2] = c_angle[1] * c_angle[0]

    return dcm


def dcm2euler(dcm):
    """
    Convert direction cosine matrix (DCM) to Euler angles
    The euler rotation sequence 3-2-1 is assumed
    Args:
        dcm (np.ndarray or list): 3x3 coordinate transformation matrix from n to b
    Returns:
        np.ndarray: 3x1 Euler angles, rad.
    """
    yaw = np.arctan2(dcm[0, 1], dcm[0, 0])  # arctan2(cy*sz, cy*cz) <=> arctan2(sz,cz)
    pitch = np.arcsin(-dcm[0, 2])
    roll = np.arctan2(dcm[1, 2], dcm[2, 2])

    return np.array([roll, pitch, yaw])


def exp_att2euler(pos, exp_att):
    """
    Convert exponential attitude to Euler Angles.

    More information about the exponential attitude conversions in:
    https://www.mecharithm.com/explicit-representation-orientation-exponential-coordinates/

    NOTE: The current implementation of this function is not the most correct because:
        input exp_att is the attitude from b to e
        output euler vector is the attitude from n to b
    There is an underlying change in the frame of the input and output attitude state.
    This function should just convert attitude from representation to the other, and not the frame.

    Args:
         pos(np.ndarray): 3x1 position vector associated with the current attitude state
         exp_att(np.ndarray): 3x1 attitude in exponential angle notation
    Returns:
        np.ndarray: 3x1 Euler angles, rad.
    """
    require_len_array(exp_att)
    require_len_array(pos)

    llh = cartesian2geodetic(*pos)
    theta = np.linalg.norm(exp_att)
    omega = exp_att / theta
    skew = vector2skew_symmetric(omega)

    # from exp att to dcm
    dcm_b_e = np.eye(3) + np.sin(theta) * skew + (1.0 - np.cos(theta)) * (skew @ skew)
    dcm_n_e = latlon2dcm_e_ned(llh[0], llh[1]).T  # from n to e
    dcm_n_b = dcm_b_e.T @ dcm_n_e  # from n to b

    # from dcm to euler
    return dcm2euler(dcm_n_b)


def euler2exp_att(pos, euler):
    """
    Convert Euler angles to exponential attitude

    More information about the exponential attitude conversions in:
    https://www.mecharithm.com/explicit-representation-orientation-exponential-coordinates/

    NOTE: The same applies to this function

    Args:
         pos(np.ndarray): 3x1 position vector associated with the current attitude state
         euler(np.ndarray): 3x1 Euler angles vector
    Returns:
        np.ndarray: 3x1 exponential attitude vector
    """
    require_len_array(euler)
    require_len_array(pos)

    llh = cartesian2geodetic(*pos)

    # from euler to dcm
    dcm_n_b = euler2dcm(euler)  # from n to b
    dcm_n_e = latlon2dcm_e_ned(llh[0], llh[1]).T  # from n to e
    dcm_b_e = dcm_n_e @ dcm_n_b.T  # from b to e

    # from dcm to exp att
    theta = np.arccos((np.trace(dcm_b_e) - 1.0) / 2.0)
    logRR = (dcm_b_e - dcm_b_e.T) * theta / (2.0 * np.sin(theta))
    # print(f"theta={theta} -> sin(theta)={sin(theta)} -> trace={np.trace(dcm_b_e)}")

    exp1 = logRR[2, 1]
    exp2 = logRR[0, 2]
    exp3 = logRR[1, 0]

    return np.array([exp1, exp2, exp3])


def euler2quat(pos, euler):
    require_len_array(euler)
    require_len_array(pos)

    llh = cartesian2geodetic(*pos)
    dcm_e_n = latlon2dcm_e_ned(llh[0], llh[1])  # from e to n
    dcm_n_b = euler2dcm(euler)  # from n to b
    dcm_e_b = dcm_n_b @ dcm_e_n

    return dcm2quat(dcm_e_b)


def quat_normalize(q):
    """
    Normalize a quaternion, scalar part is always non-negative
    Args:
        q(np.ndarray): input quaternion to be normalized
    Returns:
        (np.ndarray): normalized quaternion, scalar part is always non-negative
    """
    if q[0] < 0:
        q = -q
    q_norm = np.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3])
    qn = q / q_norm
    return qn


def dcm2quat(c):
    """
    Convert direction cosine matrix to quaternion
    Args:
        c(np.ndarray): direction cosine matrix
    Returns:
        np.ndarray: quaternion q, scalar first
    """
    require_len_matrix(c, 3, 3)

    tr = c[0, 0] + c[1, 1] + c[2, 2]
    tmp = np.array([0.0, 0.0, 0.0, 0.0])
    q = np.array([0.0, 0.0, 0.0, 0.0])
    if tr > 0.0:
        tmp[0] = 0.5 * np.sqrt(1.0 + tr)
        tmp[1] = 0.25 / tmp[0] * (c[1, 2] - c[2, 1])
        tmp[2] = 0.25 / tmp[0] * (c[2, 0] - c[0, 2])
        tmp[3] = 0.25 / tmp[0] * (c[0, 1] - c[1, 0])
    else:
        if (c[1, 1] > c[0, 0]) and (c[1, 1] > c[2, 2]):
            sqdip1 = np.sqrt(c[1, 1] - c[0, 0] - c[2, 2] + 1.0)
            tmp[2] = 0.5 * sqdip1
            if sqdip1 != 0.0:  # // if it equals 0, something is wrong
                sqdip1 = 0.5 / sqdip1
            tmp[0] = (c[2, 0] - c[0, 2]) * sqdip1
            tmp[1] = (c[0, 1] + c[1, 0]) * sqdip1
            tmp[3] = (c[1, 2] + c[2, 1]) * sqdip1
        elif c[2, 2] > c[0, 0]:
            sqdip1 = np.sqrt(c[2, 2] - c[0, 0] - c[1, 1] + 1.0)
            tmp[3] = 0.5 * sqdip1
            if sqdip1 != 0.0:  # // if it equals 0, something is wrong
                sqdip1 = 0.5 / sqdip1
            tmp[0] = (c[0, 1] - c[1, 0]) * sqdip1
            tmp[1] = (c[2, 0] + c[0, 2]) * sqdip1
            tmp[2] = (c[1, 2] + c[2, 1]) * sqdip1
        else:
            sqdip1 = np.sqrt(c[0, 0] - c[1, 1] - c[2, 2] + 1.0)
            tmp[1] = 0.5 * sqdip1
            if sqdip1 != 0.0:  # // if it equals 0, something is wrong
                sqdip1 = 0.5 / sqdip1
            tmp[0] = (c[1, 2] - c[2, 1]) * sqdip1
            tmp[2] = (c[0, 1] + c[1, 0]) * sqdip1
            tmp[3] = (c[2, 0] + c[0, 2]) * sqdip1

    # ensure q[0] is non-negative
    if tmp[0] < 0:
        q = -1.0 * tmp
    else:
        q = tmp

    return quat_normalize(q)  # quaternion normalization, *** no need if dcm is really a dcm
