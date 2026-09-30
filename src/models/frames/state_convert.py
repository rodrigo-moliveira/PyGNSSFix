import numpy as np

# TODO: add docstrings

from src.models.frames import geodetic2cartesian, cartesian2geodetic, euler2dcm, dcm2quat, quat2dcm, dcm2euler


def base_convertion(data_in: np.ndarray, operation: callable) -> np.ndarray:
    data = np.asarray(data_in)

    if data.ndim == 1:
        if data.shape[0] != 3:
            raise ValueError(
                f"Expected a vector with dimension 3, got shape {data.shape}"
            )
        was_single = True
        data = data[None, :]       # (3,) -> (1, 3)

    elif data.ndim == 2:
        if data.shape[1] != 3:
            raise ValueError(
                f"Expected each vector to have dimension 3, got shape {data.shape}"
            )
        was_single = False

    else:
        raise ValueError(
            f"Expected input with shape (3,) or (N, 3), got {data.shape}"
        )

    data_out = operation(data_in)

    return data_out[0] if was_single else data_out


def new_lld2ecef(data_in):
    data_out = np.zeros(data_in.shape)
    for t in range(data_in.shape[0]):
        # altitude is positive downwards, hence the minus sign in h
        data_out[t, :] = geodetic2cartesian(data_in[t, 0], data_in[t, 1], -data_in[t, 2])
    return data_out


def new_ecef2lld(data_in):
    data_out = np.zeros(data_in.shape)
    for t in range(data_in.shape[0]):
        data_out[t, :] = cartesian2geodetic(data_in[t, 0], data_in[t, 1], data_in[t, 2])
        data_out[t, 2] *= -1  # altitude is positive downwards
    return data_out


def new_lla2lld(data_in):
    data_out = data_in.copy()
    data_out[:, 2] *= -1
    return data_out

def new_euler2quat(data_in):
    data_out = np.zeros((len(data_in), 4))
    for t in range(data_in.shape[0]):
        euler = data_in[t, :]
        dcm = euler2dcm(euler)
        quat = dcm2quat(dcm)
        data_out[t, :] = quat
    return data_out


def new_quat2euler(data_in):
    data_out = np.zeros((len(data_in), 3))
    for t in range(data_in.shape[0]):
        quat = data_in[t, :]
        dcm = quat2dcm(quat)
        euler = dcm2euler(dcm)
        data_out[t, :] = euler
    return data_out
