import numpy as np
from src.errors import ArraySizeError


def require_len_array(arr: np.ndarray, length=3):
    if not (type(arr) == np.ndarray):
        raise TypeError(f"not a valid numpy array object (numpy.ndarray). It is of type {type(arr)}")

    if len(arr) != length or arr.size != length:
        raise ArraySizeError(f"invalid shape of array {arr}: should be an array of size {length}")


def require_len_matrix(arr: np.ndarray, nrows, ncols):
    if not (type(arr) == np.ndarray):
        raise TypeError(f"not a valid numpy array object (numpy.ndarray). It is of type {type(arr)}")

    rows, cols = arr.shape

    if rows != nrows or cols != ncols:
        raise ArraySizeError(f"invalid shape of matrix {arr}: should be of dimension {nrows}x{ncols} but it is "
                             f"a {rows}x{cols} matrix")
