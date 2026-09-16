""" Module with useful attitude functions and representation conversions """

import numpy as np

from src.utils.math_utils import require_len_array, require_len_matrix


def euler2dcm(euler):
    """
    Convert reference-to-body Euler angles to a direction cosine matrix (DCM).

    The Euler angles follow the 3-2-1 (yaw-pitch-roll) convention:
        euler = [roll, pitch, yaw] = [phi, theta, psi]

    The resulting DCM transforms vector components from the reference frame
    to the body frame:

        v_b = C_r^b @ v_r

    where:
        - r denotes the reference frame (e.g. navigation or ECEF)
        - b denotes the body frame
        - C_r^b is the output DCM from the reference frame (r) to the body frame (b)

    The DCM is given by:

        [          cy*cz,          cy*sz,            -sy]
        [ sy*sx*cz-sz*cx, sy*sx*sz+cz*cx,          cy*sx]
        [ sy*cx*cz+sz*sx, sy*cx*sz-cz*sx,          cy*cx]

    where:
        phi   = roll
        theta = pitch
        psi   = yaw

        s_x = sin(x)
        c_x = cos(x)

    Args:
        euler (numpy.ndarray): 3-element array containing [roll, pitch, yaw] in radians.
    Returns:
        numpy.ndarray: 3x3 DCM transforming vector components from the reference frame
            to the body frame (C_r^b).
    Raises:
        AttributeError: an exception is raised if the input argument does not have the correct shape
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
    Convert a reference-to-body DCM to 3-2-1 Euler angles.

    The DCM is assumed to transform vector components from the reference
    frame to the body frame:

        v_b = C_r^b @ v_r

    The Euler angles use the 3-2-1 (yaw-pitch-roll) convention:

        euler = [roll, pitch, yaw] = [phi, theta, psi]

    The angles are recovered from the DCM as:

        yaw   = atan2(C_r^b[0, 1], C_r^b[0, 0])
        pitch = asin(-C_r^b[0, 2])
        roll  = atan2(C_r^b[1, 2], C_r^b[2, 2])

    Args:
        dcm (numpy.ndarray): 3x3 DCM transforming vector components from the reference
            frame to the body frame (C_r^b).

    Returns:
        numpy.ndarray: 3x1 Euler angles, rad.
    Raises:
        AttributeError: an exception is raised if the input argument does not have the correct shape
    """
    require_len_matrix(dcm)

    yaw = np.arctan2(dcm[0, 1], dcm[0, 0])  # arctan2(cy*sz, cy*cz) <=> arctan2(sz,cz)

    # Protect against small numerical violations of |dcm[0, 2]| <= 1.
    pitch = np.arcsin(np.clip(-dcm[0, 2], -1.0, 1.0))

    roll = np.arctan2(dcm[1, 2], dcm[2, 2])

    return np.array([roll, pitch, yaw])


def quat_normalize(q):
    """
    Normalize a quaternion and enforce a non-negative scalar component.

    The quaternion uses the scalar-first convention:

        q = [q0, q1, q2, q3]

    Since a quaternion and its negative represent the same rotation, the
    sign is chosen such that q0 >= 0.

    Args:
        q (numpy.ndarray):
            4-element quaternion [q0, q1, q2, q3].

    Returns:
        numpy.ndarray:
            Unit-norm quaternion with q0 >= 0.

    Raises:
        AttributeError:
            If `q` does not have the required shape.
        ValueError:
            If `q` has zero magnitude.
    """
    require_len_array(q, length=4)

    q_norm = np.linalg.norm(q)

    if q_norm == 0.0:
        raise ValueError("Cannot normalize a zero quaternion.")

    qn = q / q_norm

    if qn[0] < 0.0:
        qn = -qn

    return qn

def quat_conjugate(q):
    """
    Compute the conjugate of a quaternion.

    The quaternion uses the scalar-first convention:

        q = [q0, q1, q2, q3]

    The conjugate is:

        q* = [q0, -q1, -q2, -q3]

    For a unit quaternion, the conjugate is also the inverse rotation.
    Therefore, if `q` represents a transformation from frame r to frame b,
    its conjugate represents the inverse transformation from frame b to
    frame r.

    Args:
        q (numpy.ndarray):
            4-element quaternion [q0, q1, q2, q3].

    Returns:
        numpy.ndarray:
            Conjugate quaternion [q0, -q1, -q2, -q3].

    Raises:
        AttributeError:
            If `q` does not have the required shape.
    """
    require_len_array(q, length=4)

    q_conj = q.copy()
    q_conj[1:] *= -1.0

    return q_conj

def quat_inverse(q):
    """
    Compute the inverse of a quaternion.

    The quaternion uses the scalar-first convention:

        q = [q0, q1, q2, q3]

    The inverse is defined as:

        q^-1 = q* / ||q||^2

    where q* is the quaternion conjugate. For a unit quaternion, the
    inverse is equal to the conjugate.

    Args:
        q (numpy.ndarray):
            4-element quaternion [q0, q1, q2, q3].

    Returns:
        numpy.ndarray:
            Inverse quaternion.

    Raises:
        AttributeError:
            If `q` does not have the required shape.
        ValueError:
            If `q` has zero magnitude.
    """
    require_len_array(q, length=4)

    q_norm_sq = np.dot(q, q)

    if q_norm_sq == 0.0:
        raise ValueError("Cannot invert a zero quaternion.")

    q_inv = q.copy()
    q_inv[1:] *= -1.0
    q_inv /= q_norm_sq

    return q_inv

def dcm2quat(c):
    """
    Convert a direction cosine matrix (DCM) to a scalar-first quaternion.

    The quaternion represents the same coordinate transformation as the
    input DCM. The transformation direction is therefore determined by the
    input DCM and is not imposed by this function.

    For example, if `c` transforms vector components from a reference frame
    to a body frame:

        v_b = C_b_r @ v_r

    then the returned quaternion represents the reference-to-body
    transformation. Conversely, if `c` is the body-to-reference DCM, the
    returned quaternion represents the body-to-reference transformation.

    The quaternion is represented using the scalar-first convention:

        q = [q0, q1, q2, q3] = [qw, qx, qy, qz]

    Since q and -q represent the same rotation, the sign is chosen such
    that q0 >= 0.

    Args:
        c (numpy.ndarray):
            3x3 direction cosine matrix representing the coordinate
            transformation to be converted.

    Returns:
        numpy.ndarray:
            4-element, unit-norm quaternion [q0, q1, q2, q3] representing
            the same coordinate transformation as `c`.

    Raises:
        AttributeError:
            If `c` does not have the required shape.
    """
    require_len_matrix(c, 3, 3)

    tr = c[0, 0] + c[1, 1] + c[2, 2]
    tmp = np.array([0.0, 0.0, 0.0, 0.0])
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

def quat_multiply(q1, q2):
    """
    Compute the Hamilton product of two quaternions.

    The quaternions use the scalar-first convention:

        q = [q0, q1, q2, q3] = [qw, qx, qy, qz]

    The Hamilton product is defined as:

        q = q1 ⊗ q2

    Quaternion multiplication is not commutative, i.e.:

        q1 ⊗ q2 != q2 ⊗ q1

    Args:
        q1 (numpy.ndarray):
            4-element quaternion [q0, q1, q2, q3].

        q2 (numpy.ndarray):
            4-element quaternion [q0, q1, q2, q3].

    Returns:
        numpy.ndarray:
            4-element quaternion representing q1 ⊗ q2.

    Raises:
        AttributeError:
            If either input does not have the required shape.
    """
    require_len_array(q1, length=4)
    require_len_array(q2, length=4)

    a0, a1, a2, a3 = q1
    b0, b1, b2, b3 = q2

    return np.array([
        a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
        a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
        a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
        a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0,
    ])

"""
def exp_att2euler(pos, exp_att):

    Convert exponential attitude to Euler Angles.

    More information about the exponential attitude conversions in:
    https://www.mecharithm.com/explicit-representation-orientation-exponential-coordinates/

    NOTE: The current implementation of this function is not the most correct because:
        input exp_att is the attitude from b to e
        output euler vector is the attitude from n to b
    There is an underlying change in the frame of the input and output attitude state.
    This function should just convert attitude from representation to the other, and not the frame.

    Args:
         pos(numpy.ndarray): 3x1 position vector associated with the current attitude state
         exp_att(numpy.ndarray): 3x1 attitude in exponential angle notation
    Returns:
        numpy.ndarray: 3x1 Euler angles, rad.
    Raises:
        AttributeError: an exception is raised if the input arguments do not have the correct shape

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
"""

"""
def euler2exp_att(pos, euler):

    Convert Euler angles to exponential attitude

    More information about the exponential attitude conversions in:
    https://www.mecharithm.com/explicit-representation-orientation-exponential-coordinates/

    NOTE: The same applies to this function

    Args:
         pos(numpy.ndarray): 3x1 position vector associated with the current attitude state
         euler(numpy.ndarray): 3x1 Euler angles vector
    Returns:
        numpy.ndarray: 3x1 exponential attitude vector
    Raises:
        AttributeError: an exception is raised if the input arguments do not have the correct shape

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
"""
"""
def euler2quat(pos, euler):

    Convert Euler angles to quaternion

    Args:
         pos(numpy.ndarray): 3x1 position vector associated with the current attitude state
         euler(numpy.ndarray): 3x1 Euler angles vector
    Returns:
        numpy.ndarray: 3x1 exponential attitude vector
    Raises:
        AttributeError: an exception is raised if the input arguments do not have the correct shape

    require_len_array(euler)
    require_len_array(pos)

    llh = cartesian2geodetic(*pos)
    dcm_e_n = latlon2dcm_e_ned(llh[0], llh[1])  # from e to n
    dcm_n_b = euler2dcm(euler)  # from n to b
    dcm_e_b = dcm_n_b @ dcm_e_n

    return dcm2quat(dcm_e_b)
"""