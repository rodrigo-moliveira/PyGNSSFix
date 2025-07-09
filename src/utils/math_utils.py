""" Module with useful mathematical functions to support the library """

__all__ = ["rot1", "rot2", "rot3", "vector2skew_symmetric", "skew_symmetric2vector", "require_len_array",
           "require_len_matrix"]

import numpy as np
from numpy import sin, cos
from src.errors import ArraySizeError


def rot1(theta):
    """
    DCM of rotation of `theta` around first axis (Rotation 1).

    Args:
        theta (float): Angle in radians
    Returns:
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
    DCM of rotation of `theta` around second axis (Rotation 2).

    Args:
        theta (float): Angle in radians
    Returns:
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
    DCM of rotation of `theta` around third axis (Rotation 3).

    Args:
        theta (float): Angle in radians
    Returns:
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
    Compute skew symmetric matrix of angular vector a.

    Args:
        a(numpy.ndarray): a (3,) array. Angular vector, in [rad]
    Returns:
        numpy.ndarray: Skew symmetric matrix, (3,3)
    Raises:
        AttributeError: an exception is raised if the input argument does not have requires shape (3,)
    """

    if a.size != 3:
        raise AttributeError(f"Vector {a} must have shape (3,). Instead, it has shape {a.shape}")

    skew = np.array([[0.0, -a[2], a[1]],
                     [a[2], 0.0, -a[0]],
                     [-a[1], a[0], 0.0]])
    return skew


def skew_symmetric2vector(skew):
    """
    Compute the angular vector associated to the provided skew symmetric matrix
    Args:
        skew(numpy.ndarray): Skew symmetric matrix, (3,3)
    Returns:
        a(numpy.ndarray): a (3,) array. Angular vector, in [rad]
    Raises:
        AttributeError: an exception is raised if the input argument does not have requires shape (3,)
    """

    if skew.shape != (3, 3):
        raise AttributeError(f"Vector {skew} must have shape (3,3). Instead, it has shape {skew.shape}")

    # note: we assume that the matrix is skew symmetric (skew^T = -skew)
    return np.array([-skew[1, 2], skew[0, 2], -skew[0, 1]])


def require_len_array(arr: np.ndarray, length=3) -> None:
    """
    Checks the length of a numpy flat array.

    Args:
        arr(numpy.ndarray): input array to be checked
        length(int): required length of the array. Default value is 3
    Raises:
        AttributeError: if the input array is not a numpy object, the exception is raised
        ArraySizeError: if the input array does not have the correct length, the exception is raised
    """
    if not (type(arr) == np.ndarray):
        raise AttributeError(f"Not a valid numpy array object (numpy.ndarray). It is of type {type(arr)}")

    if len(arr) != length or arr.size != length:
        raise ArraySizeError(f"Invalid shape of array {arr}: should be an array of size {length}")


def require_len_matrix(arr: np.ndarray, nrows=3, ncols=3):
    """
    Checks the length of a numpy array (matrix).

    Args:
        arr(numpy.ndarray): input matrix to be checked
        nrows(int): required number of rows. Default value is 3
        ncols(int): required number of columns. Default value is 3
    Raises:
        AttributeError: if the input array is not a numpy object, the exception is raised
        ArraySizeError: if the input array (matrix) does not have the correct shape, the exception is raised
    """
    if not (type(arr) == np.ndarray):
        raise AttributeError(f"Not a valid numpy array object (numpy.ndarray). It is of type {type(arr)}")

    rows, cols = arr.shape

    if rows != nrows or cols != ncols:
        raise ArraySizeError(f"Invalid shape of matrix {arr}: should be of dimension {nrows}x{ncols} but it is "
                             f"a {rows}x{cols} matrix")


def delete_state(x_in, P_in, index):
    """
    Delete the state at the given index from state vector `x_in` and covariance matrix `P_in`.
    Removes both the corresponding row and column on the covariance matrix.

    Args:
        x_in(numpy.ndarray): state vector of size n
        P_in(numpy.ndarray): covariance matrix of size (nxn)
        index(int): index of the state to be removed

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: returns the new state vector and covariance matrix without the deleted
            state (size n-1).
    """
    x_out = np.delete(x_in, index)
    P_out = np.delete(P_in, index, axis=0)  # Delete row
    P_out = np.delete(P_out, index, axis=1)  # Delete column
    return x_out, P_out


def add_state(x_in, P_in, index, new_state, new_var):
    """
    Add a new state value at the given index in state vector `x_in` and covariance matrix `P_in`.
    The new state is assumed uncorrelated with existing ones.

    Args:
        x_in(numpy.ndarray): state vector, shape (n,)
        P_in(numpy.ndarray): covariance matrix, shape (n, n)
        index(int): position to insert the new state
        new_state(float): value of the new state (scalar)
        new_var(float) : variance of the new state (scalar)

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: tuple with:
            * x_out : extended state vector, shape (n+1,)
            * P_out : extended covariance matrix, shape (n+1, n+1)
    """
    # Insert new value in state vector
    x_out = np.insert(x_in, index, new_state)

    # Create new covariance matrix with zeros
    n = P_in.shape[0]
    P_out = np.zeros((n + 1, n + 1))

    # Copy existing blocks
    P_out[:index, :index] = P_in[:index, :index]
    P_out[:index, index+1:] = P_in[:index, index:]
    P_out[index+1:, :index] = P_in[index:, :index]
    P_out[index+1:, index+1:] = P_in[index:, index:]

    # Set the variance of the new state
    P_out[index, index] = new_var

    return x_out, P_out
