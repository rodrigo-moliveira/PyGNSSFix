""" Module with useful attitude functions and representation conversions """

import numpy as np
import scipy

from src.utils.math_utils import require_len_array, require_len_matrix, skew_symmetric2vector


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

    The quaternion uses the scalar-last convention:

        q = [q1, q2, q3, q4]

    Since a quaternion and its negative represent the same rotation, the
    sign is chosen such that q4 >= 0.

    Args:
        q (numpy.ndarray):
            4-element quaternion [q1, q2, q3, q4].

    Returns:
        numpy.ndarray:
            Unit-norm quaternion with q4 >= 0.

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

    if qn[3] < 0.0:
        qn = -qn

    return qn

def quat_conjugate(q):
    """
    Compute the conjugate of a quaternion.

    The quaternion uses the scalar-last convention:

        q = [q1, q2, q3, q4]

    The conjugate is:

        q* = [-q1, -q2, -q3, q4]

    For a unit quaternion, the conjugate is also the inverse rotation.
    Therefore, if `q` represents a transformation from frame r to frame b,
    its conjugate represents the inverse transformation from frame b to
    frame r.

    Args:
        q (numpy.ndarray):
            4-element quaternion [q1, q2, q3, q4].

    Returns:
        numpy.ndarray:
            Conjugate quaternion [-q1, -q2, -q3, q4].

    Raises:
        AttributeError:
            If `q` does not have the required shape.
    """
    require_len_array(q, length=4)

    q_conj = q.copy()
    q_conj[0:3] *= -1.0

    return q_conj

def quat_inverse(q):
    """
    Compute the inverse of a quaternion.

    The quaternion uses the scalar-last convention:

        q = [q1, q2, q3, q4]

    The inverse is defined as:

        q^-1 = q* / ||q||^2

    where q* is the quaternion conjugate. For a unit quaternion, the
    inverse is equal to the conjugate.

    Args:
        q (numpy.ndarray):
            4-element quaternion [q1, q2, q3, q4].

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
    q_inv[0:3] *= -1.0
    q_inv /= q_norm_sq

    return q_inv

def dcm2quat(c):
    """
    Convert a direction cosine matrix (DCM) to a scalar-last quaternion.

    The quaternion is represented as:
        q = [q1, q2, q3, q4] = [qx, qy, qz, qw]

    The returned quaternion represents the same coordinate transformation
    as the input DCM. The sign is chosen such that q4 >= 0.

    Args:
        c (numpy.ndarray):
            3x3 direction cosine matrix.

    Returns:
        numpy.ndarray:
            4-element unit quaternion [q1, q2, q3, q4].
    """
    require_len_matrix(c, 3, 3)

    tr = c[0, 0] + c[1, 1] + c[2, 2]
    q = np.zeros(4)

    if tr > 0.0:
        q[3] = 0.5 * np.sqrt(1.0 + tr)
        q[0] = 0.25 / q[3] * (c[1, 2] - c[2, 1])
        q[1] = 0.25 / q[3] * (c[2, 0] - c[0, 2])
        q[2] = 0.25 / q[3] * (c[0, 1] - c[1, 0])

    else:
        if (c[1, 1] > c[0, 0]) and (c[1, 1] > c[2, 2]):
            sqdip1 = np.sqrt(c[1, 1] - c[0, 0] - c[2, 2] + 1.0)
            q[1] = 0.5 * sqdip1

            if sqdip1 != 0.0:
                sqdip1 = 0.5 / sqdip1

            q[3] = (c[2, 0] - c[0, 2]) * sqdip1
            q[0] = (c[0, 1] + c[1, 0]) * sqdip1
            q[2] = (c[1, 2] + c[2, 1]) * sqdip1

        elif c[2, 2] > c[0, 0]:
            sqdip1 = np.sqrt(c[2, 2] - c[0, 0] - c[1, 1] + 1.0)
            q[2] = 0.5 * sqdip1

            if sqdip1 != 0.0:
                sqdip1 = 0.5 / sqdip1

            q[3] = (c[0, 1] - c[1, 0]) * sqdip1
            q[0] = (c[2, 0] + c[0, 2]) * sqdip1
            q[1] = (c[1, 2] + c[2, 1]) * sqdip1

        else:
            sqdip1 = np.sqrt(c[0, 0] - c[1, 1] - c[2, 2] + 1.0)
            q[0] = 0.5 * sqdip1

            if sqdip1 != 0.0:
                sqdip1 = 0.5 / sqdip1

            q[3] = (c[1, 2] - c[2, 1]) * sqdip1
            q[1] = (c[0, 1] + c[1, 0]) * sqdip1
            q[2] = (c[2, 0] + c[0, 2]) * sqdip1

    # Choose the representative with non-negative scalar part.
    if q[3] < 0.0:
        q = -q

    return quat_normalize(q)  # quaternion normalization, *** no need if dcm is really a dcm

def quat2dcm(q):
    """
    Convert a scalar-last quaternion to a direction cosine matrix (DCM).

    The quaternion uses the scalar-last convention:
        q = [q1, q2, q3, q4] = [qx, qy, qz, qw]

    The resulting DCM represents the same coordinate transformation as
    the input quaternion. The frame transformation is therefore determined
    by the interpretation of the input quaternion.

    For example:
        q_n^b -> C_n^b, with v_b = C_n^b @ v_n
        q_b^n -> C_b^n, with v_n = C_b^n @ v_b

    Args:
        q (numpy.ndarray): 4-element scalar-last quaternion.

    Returns:
        numpy.ndarray: 3x3 direction cosine matrix representing the same
            coordinate transformation as the input quaternion.

    Raises:
        AttributeError: If the input quaternion does not have the correct
            shape.
    """
    require_len_array(q, length=4)

    q1, q2, q3, q4 = q

    dcm = np.zeros((3, 3))

    dcm[0, 0] = q1*q1 - q2*q2 - q3*q3 + q4*q4
    dcm[1, 0] = 2.0 * (q1*q2 - q4*q3)
    dcm[2, 0] = 2.0 * (q1*q3 + q4*q2)

    dcm[0, 1] = 2.0 * (q1*q2 + q4*q3)
    dcm[1, 1] = q4*q4 - q1*q1 + q2*q2 - q3*q3
    dcm[2, 1] = 2.0 * (q2*q3 - q4*q1)

    dcm[0, 2] = 2.0 * (q1*q3 - q4*q2)
    dcm[1, 2] = 2.0 * (q2*q3 + q4*q1)
    dcm[2, 2] = q4*q4 - q1*q1 - q2*q2 + q3*q3

    return dcm

def quat_multiply(q1, q2):
    """
    Compute the Hamilton product of two quaternions.

    The quaternions must use the scalar-last convention:

        q = [q1, q2, q3, q4] = [qx, qy, qz, qw]

    The Hamilton product is defined as:

        q = q1 ⊗ q2

    Quaternion multiplication is not commutative, i.e.:

        q1 ⊗ q2 != q2 ⊗ q1

    The convention implemented follows the same right-to-left successive rotation as the DCMs, that is:
        DCM(q1 ⊗ q2) = DCM(q1) * DCM(q2)

    where DCM(.) is the quat2dcm function.

    Args:
        q1 (numpy.ndarray):
            4-element quaternion [q1, q2, q3, q4].

        q2 (numpy.ndarray):
            4-element quaternion [q1, q2, q3, q4].

    Returns:
        numpy.ndarray:
            4-element quaternion representing q1 ⊗ q2.

    Raises:
        AttributeError:
            If either input does not have the required shape.
    """
    require_len_array(q1, length=4)
    require_len_array(q2, length=4)

    a1, a2, a3, a4 = q1
    b1, b2, b3, b4 = q2

    return np.array([
         a4*b1 + a3*b2 - a2*b3 + a1*b4,
        -a3*b1 + a4*b2 + a1*b3 + a2*b4,
         a2*b1 - a1*b2 + a4*b3 + a3*b4,
        -a1*b1 - a2*b2 - a3*b3 + a4*b4,
    ])

def quat_rotate(q, x_A):
    """
    Transform a 3D vector from frame A to frame B using a quaternion.

    The quaternion uses the scalar-last convention:
        q = [q1, q2, q3, q4] = [qx, qy, qz, qw]

    The quaternion represents the coordinate transformation from frame A
    to frame B. The vector is embedded as a pure quaternion and transformed
    using the quaternion product defined by quat_multiply():

        x_B = q × x_A_quat × q_conj

    where:
        x_A_quat = [x_A[0], x_A[1], x_A[2], 0]

    The input quaternion must be unit length for the sandwich product to
    preserve the vector's magnitude.

    Args:
        q (numpy.ndarray):
            4-element unit quaternion representing the transformation
            from frame A to frame B, in scalar-last order.

        x_A (numpy.ndarray):
            3-element vector expressed in frame A.

    Returns:
        numpy.ndarray:
            3-element vector expressed in frame B.

    Raises:
        AttributeError:
            If either input does not have the required shape.
    """
    require_len_array(q, length=4)
    require_len_array(x_A, length=3)

    x_A_quat = np.array([x_A[0], x_A[1], x_A[2], 0.0])

    q_conj = quat_conjugate(q)

    x_B_quat = quat_multiply(
        quat_multiply(q, x_A_quat),
        q_conj
    )

    return x_B_quat[:3]

def rel_att_quat(quat, quat_next, dt):
    """
    Estimate the body angular velocity from two consecutive attitude
    quaternions using the relative rotation quaternion.

    The input quaternions represent the attitude transformation from
    the navigation frame (NED) to the body frame, using the convention

        v^b = C_n^b @ v^n

    where C_n^b is the direction cosine matrix (DCM) corresponding
    to the input quaternion.

    The relative quaternion is computed as

        q_rel = q_n^{b_{k+1}} ⊗ (q_n^{b_k})^{-1}

    where ⊗ denotes quaternion multiplication. Under the quaternion
    multiplication convention used by quat_multiply(), this relative
    quaternion represents the coordinate transformation from body
    frame b_k to body frame b_{k+1}, denoted q_{b_k}^{b_{k+1}}.

    This interpretation assumes that the navigation frame is held
    fixed between epochs t_k and t_{k+1}. Consequently, the estimated
    angular velocity is the average angular velocity of the body
    relative to the navigation frame, expressed in the body frame.

    If the navigation frame is a moving local-level NED frame, its
    rotation relative to inertial space is not included in this
    estimate. To obtain the inertial angular velocity measured by an
    ideal gyroscope, the navigation-frame angular velocity must also
    be accounted for.

    The relative quaternion is converted into a rotation vector using
    the axis-angle representation. For a unit quaternion,

        q_rel = [e * sin(theta / 2), cos(theta / 2)]

    where e is the unit rotation axis and theta is the rotation angle.

    The angle is recovered as

        theta = 2 * atan2(norm(qv), q4)

    where qv = q_rel[0:3] is the vector part and q4 = q_rel[3]
    is the scalar part (scalar-last quaternion convention).

    The rotation vector is then

        theta_vec = theta * qv / norm(qv)

    for a nonzero rotation. For very small rotations, the limiting
    approximation

        theta_vec ≈ 2 * qv

    is used.

    Finally, the interval-average angular velocity is estimated as

        w = theta_vec / dt

    Quaternion sign handling:
      1. The sign of quat_next is adjusted, if necessary, so that
         dot(quat, quat_next) >= 0. This avoids artificial sign jumps
         between consecutive quaternion samples.
      2. The relative quaternion is negated, if necessary, so that
         its scalar component is non-negative. This selects the
         principal rotation angle in [0, pi].

    Assumptions:
      - Both input quaternions are valid, normalized attitude
        quaternions.
      - quat_multiply(a, b) implements the quaternion composition
        convention consistent with the attitude-to-DCM mapping.
      - dt is positive and expressed in seconds.
      - The relative rotation over the interval is represented by
        its principal rotation, with angle in [0, pi].
      - The result is an interval-average angular velocity, not
        necessarily the instantaneous angular velocity at either
        endpoint.

    Args:
        quat (np.ndarray):
            Attitude quaternion at time t_k, mapping NED coordinates
            to body-frame coordinates. Scalar-last order [qx, qy, qz, qw].
        quat_next (np.ndarray) :
            Attitude quaternion at time t_{k+1}, using the same convention
            as quat.
        dt (float) :
            Time interval t_{k+1} - t_k, in seconds. Must be positive.

    Returns:
        np.ndarray:
            Estimated interval-average angular velocity of the body
            relative to the navigation frame, expressed in body
            coordinates, in rad/s.

    """

    # Ensure that consecutive quaternion samples use consistent signs and avoids artificial jumps due to the equivalence
    # q = -q
    d = np.dot(quat, quat_next)
    if d < 0.0:
        quat_next = -quat_next

    # q_rel is q_b,k^b,k+1, that is, rotation from body frame from t(k) to t(k+1)
    q_rel = quat_multiply(quat_next, quat_inverse(quat))

    # This ensures that we use the relative quaternion with a non-negative scalar component, corresponding to a
    # principal rotation angle in [0, pi]
    if q_rel[3] < 0:
        q_rel = -q_rel

    qv = q_rel[0:3]
    qv_norm = np.linalg.norm(qv)
    q4 = q_rel[3]

    theta = 2.0 * np.arctan2(qv_norm, q4)

    if np.allclose(theta, 0.0):
        theta_vec = 2.0 * qv
    else:
        theta_vec = theta * qv / qv_norm

    # Compute w_nb^b(t_k)
    w = theta_vec / dt

    return w

def rel_att_dcm(dcm, dcm_next, dt):
    """Estimate body angular velocity from consecutive attitude DCMs.

    The input DCMs represent the transformation from navigation-frame
    (NED) coordinates to body-frame coordinates, following the
    convention:

        v^b = C_n^b @ v^n

    where C_n^b is the direction cosine matrix (DCM) at the corresponding
    epoch.

    The relative DCM is computed as:

        C_rel = C_n^{b_{k+1}} @ (C_n^{b_k}).T

    This relative DCM represents the coordinate transformation from
    body frame b_k to body frame b_{k+1}, denoted C_{b_k}^{b_{k+1}},
    assuming the navigation frame is fixed between epochs t_k and
    t_{k+1}.

    The matrix logarithm is used to extract the relative rotation:

        L = logm(C_rel)

    For a valid rotation matrix, L is skew-symmetric in exact arithmetic.
    The rotation vector is recovered using the convention:

        theta_vec = vee(-L)

    where vee() converts a skew-symmetric matrix into its corresponding
    three-dimensional vector. The negative sign follows from the
    navigation-to-body DCM convention and the definition of the
    relative DCM above.

    The interval-average angular velocity is estimated as:

        w_nb^b = theta_vec / dt

    This gives the average angular velocity of the body relative to the
    navigation frame, expressed in body coordinates, over the interval
    [t_k, t_{k+1}].

    If the navigation frame is a moving local-level NED frame, its
    rotation relative to inertial space is not included in this estimate.
    To obtain the inertial angular velocity measured by an ideal
    gyroscope, the navigation-frame angular velocity must also be
    accounted for.

    Args:
        dcm (numpy.ndarray): Attitude DCM at time t_k, with shape (3, 3).
            Transforms vectors from navigation coordinates to body
            coordinates, i.e. v_b = dcm @ v_n.
        dcm_next (numpy.ndarray): Attitude DCM at time t_{k+1}, with
            shape (3, 3), using the same convention as dcm.
        dt (float): Positive time interval t_{k+1} - t_k, in seconds.

    Returns:
        numpy.ndarray: Estimated interval-average angular velocity
            w_nb^b, expressed in body coordinates, with shape (3,)
            and units of rad/s.

    Assumptions:
        - Both input DCMs are valid rotation matrices in SO(3).
        - Both DCMs use the same navigation-frame coordinate basis.
        - dt is positive and expressed in seconds.
        - The relative rotation is represented by the principal matrix
          logarithm, with rotation angle in [0, pi].
        - The result is an interval-average angular velocity, not
          necessarily the instantaneous angular velocity at either
          endpoint.

    """
    # dcm_rel is C_b,k^b,k+1, that is, rotation matrix from body frame from t(k) to t(k+1)
    dcm_rel = dcm_next @ dcm.T

    # The relation between rotation axis and DCM is
    # A(theta_vec) = exp(- skew(theta_vec)). We invert this
    # Note: a more numerical friendly implementation may be done in the future, if we get problems
    log = scipy.linalg.logm(dcm_rel)
    theta_vec = skew_symmetric2vector(-log)

    # Compute w_nb^b(t_k)
    w = theta_vec / dt

    return w


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