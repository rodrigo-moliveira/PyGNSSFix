""" Utility functions for the mechanization equations of the INS """


import numpy as np

from src.constants import EARTH_ROTATION
from . import quat_multiply, quat_inverse
from src.models.frames import frames
from src.utils.math_utils import vector2skew_symmetric, skew_symmetric2vector

w_ie_e = np.array([0, 0, EARTH_ROTATION])  # in rad/s. This vector is also w_ie_i, although this is not needed

# Methods to compute angular velocities
def compute_w_en_n(v_eb_n, lld):
    """Compute the transport rate w_en^n (rotation of NED frame with respect
    to EARTH frame in navigation frame coordinates).

    As the vehicle travels over the curved Earth, the local North-East-Down
    (NED) frame changes orientation relative to the Earth-fixed (ECEF) frame.
    This function computes that angular rate, commonly called the *transport
    rate*, from the vehicle's NED velocity and its geodetic position.

    The three components (expressed in N, E, D axes) are:

    .. math::

        \\omega_N = \\frac{v_E}{R_n^*}

        \\omega_E = -\\frac{v_N}{R_m^*}

        \\omega_D = -\\frac{v_E \\tan\\varphi}{R_n^*}

    where

    * \\(R_m^* = R_m - d\\) and \\(R_n^* = R_n - d\\) are the *effective*
      radii of curvature (see :func:`get_earth_radii`), accounting for the
      vehicle's altitude (note that *d* is positive in the down direction)
    * \\(\\varphi\\) is the geodetic latitude in radians.

    Args:
        v_eb_n (numpy.ndarray): Velocity of the body with respect to the ECEF frame,
            expressed in NED axes, shape ``(3,)`` in m/s.
            Order: ``[v_N, v_E, v_D]``.
        lld (numpy.ndarray) : Geodetic position, shape ``(3,)``:
            ``[lat (rad), lon (rad), down (m)]``.

    Returns:
        tuple : ``(w_en_n, rm_effective, rn_effective)``:

        * ``w_en_n`` – Transport rate w_en^n, shape ``(3,)`` in rad/s,
          expressed in NED axes: ``[w_N, w_E, w_D]``.
        * ``rm_effective`` – Effective meridian radius of curvature
          (Rm − d), in metres.
        * ``rn_effective`` – Effective prime-vertical radius of curvature
          (Rn − d), in metres.
    """
    v_n, v_e, _ = v_eb_n[:]
    lat, _, d = lld[:]

    rm, rn = frames.get_earth_radii(lat)

    rm_effective = rm - d
    rn_effective = rn - d

    w_en_n = np.zeros(3)

    w_en_n[0] = v_e / rn_effective  # wN
    w_en_n[1] = -v_n / rm_effective  # wE
    w_en_n[2] = -v_e * np.tan(lat) / rn_effective  # wD

    return w_en_n, rm_effective, rn_effective


def compute_w_ie_e():
    """Return the Earth's rotation rate vector expressed in the ECEF frame,
    that is the rotation of the e frame with respect to i frame, in e frame
    coordinates.

    Returns:
        numpy.ndarray: of shape ``(3,)`` in rad/s – the constant vector
        ``[0, 0, EARTH_ROTATION]``.
    """
    return w_ie_e


def compute_w_ie_n(lat):
    """Project the Earth rotation rate into the local NED frame.

    Rotates the constant ECEF Earth-rate vector
    ``[0, 0, w_ie]`` into the North-East-Down frame at the given geodetic
    latitude:

    Args:
        lat (float): Geodetic latitude in **radians**.

    Returns:
        numpy.ndarray: array of shape ``(3,)`` in rad/s – Earth rotation rate
        expressed in NED axes: ``[w_N, w_E, w_D]``.
    """
    return np.array([np.cos(lat),
                     0,
                     -np.sin(lat)]) * EARTH_ROTATION


def compute_w_nb_b_euler(euler_dot, euler):
    """Compute the body angular rate w_nb^b from Euler-angle rates,
    that is, implements the attitude kinematics equation in Euler-angle
    state space.

    Converts the time-derivative of the 3-2-1 (yaw-pitch-roll / Z-Y-X
    intrinsic) Euler angles into the corresponding body-frame angular
    velocity vector using the standard kinematic relationship:

        w_x = phi_dot - psi_dot * np.sin(theta)
        w_y = theta_dot * cos_phi + psi_dot * cos_theta * sin_phi
        w_z = psi_dot * cos_phi * cos_theta - theta_dot * sin_phi

    Args:
        euler_dot (numpy.ndarray): Time derivative of the Euler angles, shape ``(3,)`` in
            rad/s.  Order: ``[φ̇ (roll rate), θ̇ (pitch rate), ψ̇ (yaw rate)]``.
        euler (numpy.ndarray): Current Euler angles, shape ``(3,)`` in radians.
            Order: ``[φ (roll), θ (pitch), ψ (yaw)]``.

    Returns:
        numpy.ndarray : array of shape ``(3,)`` in rad/s – angular velocity of
        the body frame with respect to the navigation frame, expressed in
        body axes: ``[w_x^b, w_y^b, w_z^b]``.

    Note:
        This formulation is singular at θ = ±90° (gimbal lock), just as
        the Euler-angle parameterisation itself.
    """

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

def compute_w_nb_b_dcm(dcm_dot, dcm_n_b):
    """Compute the body angular rate w_nb^b from the DCM and its derivative,
    that is, implements the attitude kinematics equation in DCM
    state space.

    Uses the kinematic DCM equation

        d(C_n^b)/dt = C_n^b * Skew_bn^b

    and solves for the angular-rate vector by inverting the skew operation:

        Skew_nb^b = -d(C_n^b)/dt * (C_n^b)^T

    Args:
        dcm_dot (numpy.ndarray): Time derivative of the direction-cosine matrix C_n^b
            (navigation-to-body), shape ``(3, 3)``.
        dcm_n_b (numpy.ndarray): Current direction-cosine matrix C_n^b, shape ``(3, 3)``.

    Returns:
        (numpy.ndarray) : array of shape ``(3,)`` in rad/s – angular velocity of
        the body frame with respect to the navigation frame, expressed in
        body axes.

    Note:
        The result is extracted from the skew-symmetric matrix via
        :func:`skew_symmetric2vector`, which reads the three independent
        off-diagonal elements.
    """
    skew_nb_b = -dcm_dot @ dcm_n_b.T
    w_nb_b = skew_symmetric2vector(skew_nb_b)

    return w_nb_b

def compute_w_nb_b_quat(quat_dot, q_n_b):
    """Compute the body angular rate w_nb^b from the quaternion and its derivative,
    that is, implements the attitude kinematics equation in quaternion
    state space.

    Uses the kinematic quaternion equation

        q_dot = 1/2 * w_nb^b ⊗ q_n^b

    and isolates the rate vector:

        w_nb^b = 2 * q_dot ⊗ q_n^b^{-1}

    The scalar (fourth) component of the resulting quaternion product is
    zero by construction and is discarded.

    Args:
        quat_dot (numpy.ndarray): Time derivative of the unit quaternion q_n^b, shape
            ``(4,)``.  Convention: scalar last, that is ``[x, y, z, w]``
        q_n_b (numpy.ndarray): Current unit quaternion q_n^b, shape ``(4,)``

    Returns:
        numpy.ndarray : array of shape ``(3,)`` in rad/s – angular velocity of
        the body frame with respect to the navigation frame, expressed in
        body axes.

    """
    w = 2.0 * quat_multiply(quat_dot, quat_inverse(q_n_b))
    w_nb_b = w[0:3]
    return w_nb_b

###################################
# N-Frame Mechanization Equations #
###################################

def compute_w_ib_b(c_nb, w_ie_n, w_en_n, w_nb_b):
    """Compute the ideal (true) gyroscope readout via N-frame mechanization.

    The gyroscope measures the total angular rate of the body frame (b)
    with respect to the inertial frame (i), expressed in body axes:

        w_ib^b = C_nb^b (w_ie^n + w_en^n) + w_nb^b

    where the chain of frames is **I → E → N → B**:

    * ``w_ie^n`` – Earth rotation rate expressed in the NED frame
      (see :func:`compute_w_ie_n`).
    * ``w_en^n`` – Transport rate due to vehicle motion over the curved
      Earth surface (see :func:`compute_w_en_n`).
    * ``w_nb^b`` – Body angular rate with respect to the NED frame
      (see :func:`compute_w_nb_b_euler` / ``_dcm`` / ``_quat``).
    * ``C_n^b`` – Direction-cosine matrix from NED to body axes.

    Args:
        c_nb (numpy.ndarray): DCM from navigation (NED) to body frame, shape ``(3, 3)``.
        w_ie_n (numpy.ndarray): Earth rotation rate in NED axes, shape ``(3,)`` in rad/s.
        w_en_n (numpy.ndarray): Transport (NED-frame) rate due to vehicle motion,
            shape ``(3,)`` in rad/s.
        w_nb_b (numpy.ndarray): Body angular rate with respect to the NED frame,
            expressed in body axes, shape ``(3,)`` in rad/s.

    Returns:
        numpy.ndarray: array of shape ``(3,)`` in rad/s – the ideal gyroscope
        readout (specific angular rate of b wrt i, expressed in b).
    """
    # w_ib^b = C_nb^b (w_ie^n + w_en^n) + w_nb^b
    return c_nb @ w_ie_n + c_nb @ w_en_n + w_nb_b


def compute_f_ib_b(c_nb, c_en, w_en_n, w_ie_n, v_eb_n_dot, g_eb_e, v_eb_n):
    """Compute the ideal (true) accelerometer readout via N-frame mechanization.
    -> explain that g is gravity with respect to E and not I
    The accelerometer measures the specific force (non-gravitational
    acceleration) of the body with respect to the inertial frame, expressed
    in body axes.  In the NED frame the equation of motion gives:

        f_ib_b = c_nb * (v_eb_n_dot - c_en @ g_eb_e + (skew_en_n + 2 * skew_ie_n) @ v_eb_n)

    and the body-frame readout is obtained by projecting into the body
    axes:

    where

    * ``d(v_eb^n)/dt`` – time derivative of the body velocity in the NED frame
      (computed from the reference trajectory - differentiation).
    * ``g_eb^e`` – local gravity vector, that is, the gravity of the body with
        respect to the ECEF frame (Coriolis effect must be removed from inertial gravity)
    * ``w_en^n x v`` – transport-rate Coriolis term (curvature of the
      local-level frame).
    * ``2 w_ie^n x v`` – Earth-rotation Coriolis term.

    Args:
        c_nb: DCM from navigation (NED) to body frame, shape ``(3, 3)``.
        c_en: DCM from ECEF to NED frame, shape ``(3, 3)``.
        w_en_n: Transport rate in NED axes, shape ``(3,)`` in rad/s.
        w_ie_n: Earth rotation rate in NED axes, shape ``(3,)`` in rad/s.
        v_eb_n_dot: Time derivative of the ECEF velocity expressed in NED,
            shape ``(3,)`` in m/s².
        g_eb_e: local gravity vector in ECEF axes,
            shape ``(3,)`` in m/s².
        v_eb_n: Velocity of the body with respect to ECEF, expressed in
            NED axes, shape ``(3,)`` in m/s.

    Returns:
        np.ndarray : array of shape ``(3,)`` in m/s² – the ideal accelerometer
        readout (specific force of b wrt i, expressed in b).
    """
    # accelerometer readout
    skew_en_n = vector2skew_symmetric(w_en_n)
    skew_ie_n = vector2skew_symmetric(w_ie_n)
    f_ib_n = v_eb_n_dot - c_en @ g_eb_e + (skew_en_n + 2 * skew_ie_n) @ v_eb_n
    f_ib_b = c_nb @ f_ib_n  # accelerometer readout
    return f_ib_b

def compute_imu_readouts(lld, w_nb_b, v_eb_n, v_eb_n_dot, c_nb):
    """Compute ideal (error-free) IMU readouts via N-frame mechanization.

        Given the vehicle's geodetic position, attitude, and velocity state,
        this function evaluates the full N-frame mechanization equations to
        produce the gyroscopic angular rate and the specific-force (accelerometer)
        readout that a perfect IMU would observe.

        Pipeline
        --------
        1. **Frame geometry** – Build the ECEF→NED rotation matrix ``C_en``
           from the current latitude/longitude and compute the local gravity
           vector at the vehicle's ECEF position.
        2. **Angular rates** – Evaluate the Earth-rotation rate in NED
           (w_ie^n), the transport rate from vehicle motion (w_en^n), and
           accept the body rate w_nb^b as input.
        3. **Mechanization** – Combine the rates and velocities through
           :func:`compute_w_ib_b` and :func:`compute_f_ib_b` to obtain the
           body-frame readouts.

        Args:
            lld (np.ndarray): Geodetic position, shape ``(3,)``:
                ``[lat (rad), lon (rad), down (m)]``.
                **Convention:** altitude is *positive downwards* (i.e. depth
                below the ellipsoid).  For an aircraft or satellite above the
                ellipsoid, pass ``alt`` as a **negative** value.
            w_nb_b (np.ndarray): Body angular rate with respect to the navigation frame,
                expressed in body axes, shape ``(3,)`` in rad/s.
            v_eb_n (np.ndarray): Velocity of the body with respect to ECEF, expressed in
                NED axes, shape ``(3,)`` in m/s.  Order:
                ``[v_N, v_E, v_D]``.
            v_eb_n_dot (np.ndarray): Time derivative of *v_eb_n*, shape ``(3,)`` in m/s².
            c_nb (np.ndarray): DCM from NED to body frame, shape ``(3, 3)``.

        Returns:
            tuple [np.ndarray,np.ndarray] : A tuple ``(w_ib_b, f_ib_b)``:
            * ``w_ib_b`` – Ideal gyroscope readout (specific angular rate of
              B wrt I, in body axes), shape ``(3,)`` in rad/s.
            * ``f_ib_b`` – Ideal accelerometer readout (specific force of
              B wrt I, in body axes), shape ``(3,)`` in m/s².
        """
    # Get rotation vectors between frames i-e, e-n and n-b, and local gravity vector
    c_en = frames.latlon2dcm_e_ned(lld[0], lld[1])  # matrix from e to n
    w_ie_n = compute_w_ie_n(lld[0])
    w_en_n, _, _ = compute_w_en_n(v_eb_n, lld)
    r_eb_e = np.array(
        frames.geodetic2cartesian(lld[0], lld[1], -lld[2]))  # remember that altitude is positive downwards!
    g_eb_e = frames.grav_acceleration(r_eb_e, mode="earth")

    # Apply the mechanization equations (in the n-frame)
    w_ib_b = compute_w_ib_b(c_nb, w_ie_n, w_en_n, w_nb_b)
    f_ib_b = compute_f_ib_b(c_nb, c_en, w_en_n, w_ie_n, v_eb_n_dot, g_eb_e, v_eb_n)

    return w_ib_b, f_ib_b