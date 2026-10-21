"""Module that handles convertion between state-space representations of Position, Velocity and Attitude (INS algorithms)"""
import numpy as np
from collections.abc import Callable

from src.models.frames import geodetic2cartesian, cartesian2geodetic, euler2dcm, dcm2quat, quat2dcm, dcm2euler, latlon2dcm_e_ned


def base_convertion(
    data_in: np.ndarray,
    operation: Callable[..., np.ndarray],
    *args,
) -> np.ndarray:
    """Apply a conversion operation to one or more 3D, 4D or 9D vectors.

    Converts the input data to a NumPy array, validates its dimensions,
    and passes it to the specified operation along with any additional
    positional arguments. Single-vector inputs are returned as a single
    vector, while batch inputs retain their two-dimensional shape.

    Args:
        data_in: Input vector with shape (3,), (4,) or (9,), or a batch of
            vectors with shape (N, 3) or (N, 4) or (N, 9).
        operation: Callable that performs the conversion. It must accept
            the input array and any additional positional arguments.
        *args: Additional positional arguments passed to the operation.

    Returns:
        Converted vector for a single-vector input, or a two-dimensional
        array for a batch input. The output dimensions depend on the
        operation.

    Raises:
        ValueError: If the input does not have one or two dimensions,
            or if its vector dimension is not 3 or 4 or 9.
    """
    data = np.asarray(data_in)

    if data.ndim == 1:
        if data.shape[0] not in (3, 4, 9):
            raise ValueError(
                f"Expected a vector with dimension 3 or 4 or 9, got shape {data.shape}"
            )
        was_single = True
        data = data[None, :]

    elif data.ndim == 2:
        if data.shape[1] not in (3, 4, 9):
            raise ValueError(
                f"Expected each vector to have dimension 3 or 4 or 9, got shape {data.shape}"
            )
        was_single = False

    else:
        raise ValueError(
            f"Expected input with shape (3,), (4,) (9,), (N, 3), (N, 4) or (N, 9), "
            f"got {data.shape}"
        )

    data_out = operation(data, *args)

    return data_out[0] if was_single else data_out


def sc_lld2ecef(data_in):
    """Convert geodetic coordinates to Earth-Centered, Earth-Fixed coordinates.

    Converts each input vector from geodetic latitude, longitude, and
    altitude to Earth-Centered, Earth-Fixed (ECEF) Cartesian coordinates.
    Altitude is assumed to be positive downwards and is negated before
    conversion.

    Args:
        data_in: Array of geodetic coordinates with shape (N, 3), where
            each row contains latitude, longitude, and altitude, in that
            order. Units must match those expected by
            ``geodetic2cartesian``.

    Returns:
        NumPy array of ECEF Cartesian coordinates with shape (N, 3).
    """
    data_out = np.zeros(data_in.shape)
    for t in range(data_in.shape[0]):
        # altitude is positive downwards, hence the minus sign in h
        data_out[t, :] = geodetic2cartesian(data_in[t, 0], data_in[t, 1], -data_in[t, 2])
    return data_out


def sc_ecef2lld(data_in):
    """Convert Earth-Centered, Earth-Fixed coordinates to geodetic coordinates

    Converts each input vector from Earth-Centered, Earth-Fixed (ECEF)
    Cartesian coordinates to geodetic latitude, longitude, and
    altitude.
    Altitude is assumed to be positive downwards and is negated after
    conversion.

    Args:
        data_in: Array of geodetic coordinates with shape (N, 3), where
            each row contains x, y and z coordinates. Units must match
            those expected by ``cartesian2geodetic``.

    Returns:
        NumPy array of geodetic coordinates with shape (N, 3).
    """
    data_out = np.zeros(data_in.shape)
    for t in range(data_in.shape[0]):
        data_out[t, :] = cartesian2geodetic(data_in[t, 0], data_in[t, 1], data_in[t, 2])
        data_out[t, 2] *= -1  # altitude is positive downwards
    return data_out


def sc_lla2lld(data_in):
    """Convert altitude sign convention from positive-up to positive-down.

    Copies the input geodetic coordinates and reverses the sign of the
    altitude component.

    Args:
        data_in: Array of geodetic coordinates with shape (N, 3), where
            each row contains latitude, longitude, and altitude, in that
            order. The input altitude is positive upwards.

    Returns:
        Copy of the input array with altitude positive downwards.
        The input array remains unchanged.
    """
    data_out = data_in.copy()
    data_out[:, 2] *= -1
    return data_out

def sc_euler2quat(data_in):
    """Convert Euler angles to quaternions.

    Converts each input Euler-angle vector to a quaternion using the
    Euler-angle-to-DCM and DCM-to-quaternion transformations.

    Args:
        data_in: Array of Euler angles with shape (N, 3), where each row
            contains roll, pitch, and yaw, respectively. The angles define
            the rotation from the reference frame to the body frame,
            following the 3-2-1 (yaw-pitch-roll) rotation sequence.

    Returns:
        Array of quaternions with shape (N, 4), using the scalar-last
        convention [qx, qy, qz, qw].
    """
    data_out = np.zeros((len(data_in), 4))
    for t in range(data_in.shape[0]):
        euler = data_in[t, :]
        dcm = euler2dcm(euler)
        quat = dcm2quat(dcm)
        data_out[t, :] = quat
    return data_out


def sc_quat2euler(data_in):
    """Convert quaternions to Euler angles.

    Converts each input quaternion to Euler angles using the
    quaternion-to-DCM and DCM-to-Euler-angle transformations.

    Args:
        data_in: Array of quaternions with shape (N, 4), using the
            scalar-last convention [qx, qy, qz, qw].

    Returns:
        Array of Euler angles with shape (N, 3), where each row contains
        roll, pitch, and yaw, respectively. The angles describe the
        rotation from the reference frame to the body frame, following
        the 3-2-1 (yaw-pitch-roll) rotation sequence.
    """
    data_out = np.zeros((len(data_in), 3))
    for t in range(data_in.shape[0]):
        quat = data_in[t, :]
        dcm = quat2dcm(quat)
        euler = dcm2euler(dcm)
        data_out[t, :] = euler
    return data_out

def sc_vel_n2e(data_in, pos_lld):
    """Transform velocity vectors from NED to ECEF coordinates.

    Transforms each velocity vector from the local North-East-Down (NED)
    frame to the Earth-Centered, Earth-Fixed (ECEF) frame using the
    geodetic latitude and longitude at the corresponding position.

    Args:
        data_in: Array of velocity vectors in the NED frame, with shape
            (N, 3). Each row contains north, east, and down velocity
            components.
        pos_lld: Array of geodetic positions with shape (N, 3). Each row
            contains latitude, longitude, and altitude, in that order.
            Latitude and longitude units must match those expected by
            ``latlon2dcm_e_ned``.

    Returns:
        Array of velocity vectors in the ECEF frame, with shape (N, 3).
        Each row contains the ECEF x, y, and z velocity components.
    """
    data_out = np.zeros((len(data_in), 3))
    for t in range(data_in.shape[0]):
        v_eb_n = data_in[t, :]
        lld = pos_lld[t, :]
        dcm_e_ned = latlon2dcm_e_ned(lld[0], lld[1]) # DCM from the ECEF frame to the NED frame
        v_eb_e = dcm_e_ned.T @ v_eb_n # v_eb^e = (dcm_e^n).T @ v_eb^n
        data_out[t, :] = v_eb_e

    return data_out

def sc_vel_e2n(data_in, pos_lld):
    """Transform velocity vectors from ECEF to NED coordinates.

    Transforms each velocity vector from the Earth-Centered,
    Earth-Fixed (ECEF) frame to the local North-East-Down (NED) frame
    using the geodetic latitude and longitude at the corresponding
    position.

    Args:
        data_in: Array of velocity vectors in the ECEF frame, with shape
            (N, 3). Each row contains the ECEF x, y, and z velocity
            components.
        pos_lld: Array of geodetic positions with shape (N, 3). Each row
            contains latitude, longitude, and altitude, in that order.
            Latitude and longitude units must match those expected by
            ``latlon2dcm_e_ned``.

    Returns:
        Array of velocity vectors in the NED frame, with shape (N, 3).
        Each row contains north, east, and down velocity components.
    """
    data_out = np.zeros((len(data_in), 3))
    for t in range(data_in.shape[0]):
        v_eb_e = data_in[t, :]
        lld = pos_lld[t, :]
        dcm_e_ned = latlon2dcm_e_ned(lld[0], lld[1]) # DCM from the ECEF frame to the NED frame
        v_eb_n = dcm_e_ned @ v_eb_e # v_eb^n = dcm_e^n @ v_eb^e
        data_out[t, :] = v_eb_n

    return data_out

def sc_euler2dcm(data_in):
    """Convert Euler angles to DCM.

    Converts each input Euler angles to a DCM matrix using the
    Euler-to-DCM-angle transformations.

    Args:
        data_in: Array of Euler angles with shape (N, 3), where each row contains
        roll, pitch, and yaw, respectively. The angles describe the
        rotation from the reference frame to the body frame, following
        the 3-2-1 (yaw-pitch-roll) rotation sequence.

    Returns:
        Array of DCM attitude with shape (N, 9), where each row contains
        the flattened 3x3 DCM matrix. It can be converted back to a
        3x3 numpy array using the ``reshape(3, 3)`` method.
    """
    data_out = np.zeros((len(data_in), 9))
    for t in range(data_in.shape[0]):
        euler = data_in[t, :]
        dcm = euler2dcm(euler)
        data_out[t, :] = dcm.flatten()
    return data_out