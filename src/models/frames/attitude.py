import numpy as np

from .frames import cartesian2geodetic, latlon2dcm_e_ned
from src.utils.math_utils import vector2skew_symmetric, require_len_array, require_len_matrix


def euler2dcm(euler):
    require_len_array(euler)

    """
    Convert Euler angles to direction cosine matrix (DCM). The conventional rotation sequence z-y-x is adopted.
    
    Args:
        angles: 3x1 Euler angles (roll, pitch, yaw), rad.
    Returns:
        dcm: 3x3 coordinate transformation matrix from n to b
    """
    # the DCM matrix is the following
    #     [          cy*cz,          cy*sz,            -sy]
    #     [ sy*sx*cz-sz*cx, sy*sx*sz+cz*cx,          cy*sx]
    #     [ sy*cx*cz+sz*sx, sy*cx*sz-cz*sx,          cy*cx]
    # where (x,y,z) are (roll,pitch,yaw), s is short for sin and c is short for cos

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
        dcm: 3x3 coordinate transformation matrix from n to b
    Returns:
        angles: 3x1 Euler angles, rad.
    """
    # the DCM matrix is the following
    #     [          cy*cz,          cy*sz,            -sy]
    #     [ sy*sx*cz-sz*cx, sy*sx*sz+cz*cx,          cy*sx]
    #     [ sy*cx*cz+sz*sx, sy*cx*sz-cz*sx,          cy*cx]
    # where (x,y,z) are (roll,pitch,yaw), s is short for sin and c is short for cos

    yaw = np.arctan2(dcm[0, 1], dcm[0, 0])  # arctan2(cy*sz, cy*cz) <=> arctan2(sz,cz)
    pitch = np.arcsin(-dcm[0, 2])
    roll = np.arctan2(dcm[1, 2], dcm[2, 2])

    return np.array([roll, pitch, yaw])


def exp_att2euler(pos, exp_att):
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
    exponential attitude: https://www.mecharithm.com/explicit-representation-orientation-exponential-coordinates/
    Args:
         pos : ecef position
         euler : euler angles from nav to body frame
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
        q: quaternion
    Returns:
        qn: normalized quaternion, scalar part is always non-negative
    """
    if q[0] < 0:
        q = -q
    q_norm = np.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3])
    qn = q / q_norm
    return qn


def dcm2quat(c):
    require_len_matrix(c, 3, 3)

    """
    Convert direction cosine matrix to quaternion
    Args:
        c: direction cosine matrix
    Returns:
        q: quaternion, scalar first
    """
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
