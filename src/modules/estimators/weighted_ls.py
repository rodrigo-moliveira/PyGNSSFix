""" Module with the implementation of the Weighted Least-Squares solver """

import numpy as np


class WeightedLeastSquares:
    """
    Solve the Weighted Least-Squares (WLS) equation of the linear system y = G @ x, where
        * y - observation vector
        * x - state vector
        * G - state matrix (relates state vector with the observation vector)

    Given the vector of observations y, the estimated state x^{hat} is given by:
        x^{hat} = (G.T @ W @ G)^-1 @ G.T @ W @ y

    where W is the weight matrix

    The provided vectors must have the following dimensions
        dim(x) = nx1
        dim(y) = mx1
        dim(G) = mxn
        dim(W) = mxm
    If an inconsistency in these shapes is detected, an AttributeError exception is raised
    """

    def __init__(self, y: np.ndarray, G: np.ndarray, W: np.ndarray = None):
        """
        Constructor of the WeightedLeastSquares object.

        Args:
            y (numpy.ndarray) : the observation vector (mx1)
            G (numpy.ndarray) : the state matrix (mxn)
            W (numpy.ndarray or None) : the weight matrix. If not provided, then the solution resorts to a
                classical Least Squares, without weighting the observations (mxm)

        Raises:
            AttributeError : if the dimensions of the provided arrays do not match, or if the input vectors are not
                            numpy ndarray objects
        """
        # Convert y to a 2D column vector if it is 1D
        if y.ndim == 1:
            y = y[:, np.newaxis]

        # Check the shapes of y and G
        my, m1 = y.shape
        mg, ng = G.shape

        if m1 != 1:
            raise AttributeError(f"Observation vector y should be (m, 1) but got {y.shape}")
        if mg != my:
            raise AttributeError(
                f"Observation vector y and state matrix G have inconsistent shapes (number of rows), y has shape "
                f"{y.shape} and G {G.shape}.")

        # If W is None, set it to an identity matrix of appropriate size
        if W is None:
            W = np.eye(my)
        else:
            mw, mw1 = W.shape
            if mw != mw1 or mw != my:
                raise AttributeError(
                    f"Weight matrix W should be (m, m) and match the number of rows of y (m={my}), "
                    f"but got shape {W.shape}")

        # all sizes are verified
        self._m = mg
        self._n = ng

        self._x = None
        self._y = y
        self._G = G
        self._W = W
        self._S = None  # S = (G.T @ W @ G)^-1

    def solve(self):
        """
        Solves the Weighted Least Squares normal equation.

        Raises:
            LinAlgError: this function calls the `numpy.linalg.inv`.
                This exception is raised if the matrix inversion fails.
        """
        self._S = np.linalg.inv(self._G.T @ self._W @ self._G)
        self._x = self._S @ self._G.T @ self._W @ self._y

    def get_solution(self):
        """
        Returns:
            numpy.ndarray: returns the solution array `x`
        """
        return self._x.flatten()

    def get_cov(self):
        """
        Returns:
            numpy.ndarray: returns the covariance matrix `S = (G.T @ W @ G)^-1`
        """
        return self._S
