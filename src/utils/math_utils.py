import numpy as np
from numpy import sin, cos


def rot1(theta):
    """
    Args:
        theta (float): Angle in radians
    Return:
        numpy.ndarray : 3x3 rotation matrix of angle theta around the X-axis
    """
    return np.array(
        [
            [1, 0, 0],
            [0, cos(theta), sin(theta)],
            [0, -sin(theta), cos(theta)],
        ]
    )


def rot2(theta):
    """
    Args:
        theta (float): Angle in radians
    Return:
        numpy.ndarray : 3x3 rotation matrix of angle theta around the Y-axis
    """
    return np.array(
        [
            [cos(theta), 0, -sin(theta)],
            [0, 1, 0],
            [sin(theta), 0, cos(theta)],
        ]
    )


def rot3(theta):
    """
    Args:
        theta (float): Angle in radians
    Return:
        numpy.ndarray : 3x3 rotation matrix of angle theta around the Z-axis
    """
    return np.array(
        [
            [cos(theta), sin(theta), 0],
            [-sin(theta), cos(theta), 0],
            [0, 0, 1],
        ]
    )


def vector2skew_symmetric(a):
    """
    Compute skew symmetric matrix of angular vector a
    Args:
        a: (3,) array. Angular vector, in [rad]
    Returns:
        skew: (3,3) Skew symmetric matrix
    """

    if a.size != 3:
        raise TypeError(f"Vector {a} must have shape (3,). Instead, it has shape {a.shape}")

    skew = np.array([[0.0, -a[2], a[1]],
                     [a[2], 0.0, -a[0]],
                     [-a[1], a[0], 0.0]])
    return skew


def skew_symmetric2vector(skew):
    """
    Compute the angular vector associated to the provided skew symmetric matrix
    Args:
        skew: (3,3) Skew symmetric matrix
    Returns:
        a: (3,) array. Angular vector, in [rad]
    """

    if skew.shape != (3, 3):
        raise TypeError(f"Vector {skew} must have shape (3,3). Instead, it has shape {skew.shape}")

    # note: we assume that the matrix is skew symmetric (skew^T = -skew)
    return np.array([-skew[1, 2], skew[0, 2], -skew[0, 1]])
