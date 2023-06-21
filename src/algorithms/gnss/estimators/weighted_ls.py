import numpy as np


class WeightedLeastSquares:
    """
    Solve Weighted Least Squares equation of the linear system y = G @ x, where
        * y - observation vector
        * x - state vector
        * G - state matrix (relates state vector with the observation vector)

    Given the vector of observations y, the estimated state x^{hat} is given by:
        * x^{hat} = (G.T @ W @ G)^-1 @ G.T @ W @ y

    where:
        * W is the weight matrix


    in the case of iterated least squares, the system of equations takes the form dy = G @ dx, where
        * dy are the prefit residuals (observation - predicted_observation)
        * dx is the incremental state (to be added to the initial guess)


    The provided vectors must have the following dimensions
        dim(x) = nx1
        dim(y) = mx1
        dim(G) = mxn
        dim(W) = mxm
    If an inconsistency in these shapes is detected, a ValueError exception is raised
    """

    def __init__(self, y: np.ndarray, G: np.ndarray, W: np.ndarray = None):
        """

        Args:
            y (numpy.ndarray) : the observation vector
            G (numpy.ndarray) : the state matrix
            W (numpy.ndarray, optional) : the weight matrix. If not provided, then the solution resorts to a
                                        classical Least Squares, without weighting the observations

        Raises:
            AttributeError : if the dimensions of the provided arrays do not match, or if the input vectors are not
                            numpy ndarray objects
        """
        self._check_args(y, G, W)

        self._x = None
        self._y = y
        self._G = G
        self._W = W if W is not None else np.eye(self._m)
        self._S = None  # S = (G.T @ W @ G)^-1
        self._postfit_res = None  # postfit residuals = y - G @ x^{hat}

    def _check_args(self, y, G, W):
        if not isinstance(y, np.ndarray):
            raise AttributeError(f"Parameter 'y' must be of type {np.ndarray}")
        if not isinstance(G, np.ndarray):
            raise AttributeError(f"Parameter 'G' must be of type {np.ndarray}")
        if W is not None and not isinstance(W, np.ndarray):
            raise AttributeError(f"Parameter 'W' must be of type {np.ndarray}")

        # unpack shapes
        shape_y = list(y.shape)  # y vector should have shape m x 1
        shape_G = list(G.shape)  # G matrix should have shape m x n

        if len(shape_y) != 1:
            raise AttributeError(f"Parameter 'y' must be unidimensional, that is, of shape mx1. "
                                 f"Provided shape is {y.shape}")
        _m_y = shape_y[0]
        if len(shape_G) != 2:
            raise AttributeError(f"Parameter 'G' must be a 2D matrix, that is, of shape mxn. "
                                 f"Provided shape is {G.shape}")
        _m_G = shape_G[0]
        _n_G = shape_G[1]

        if _m_y != _m_G:
            raise AttributeError(f"Parameters 'y' and 'G' should have consistent shapes (mx1 and mxn, respectively)."
                                 f"'y' is of shape {y.shape} and 'G' is of shape {G.shape}")
        self._m = _m_y  # m (observation dimension)
        self._n = _n_G  # n (state dimension)

        if W is not None:
            _raise = False
            if len(list(W.shape)) != 2:
                _raise = True
            else:
                _m_W, _n_W = W.shape  # W matrix should have shape m x m
                if not (_m_W == self._m and _n_W == self._m):
                    _raise = True

            if _raise:
                raise AttributeError(f"Parameter 'W' must be a square matrix with dimension mxm, consistent with 'y' "
                                     f"of shape mx1. "
                                     f"'W' is of shape {W.shape} whereas y is of shape {y.shape}.")

    def solve(self):

        self._S = np.linalg.inv(self._G.T @ self._W @ self._G)
        self._x = self._S @ self._G.T @ self._W @ self._y

    def get_solution(self):
        return self._x
