""" Module with the implementation of the Extended Kalman Filter solver """

import numpy as np


class EKF:
    """
    A simplified implementation of the Extended Kalman Filter (EKF) for discrete-time systems.

    Unlike the typical EKF formulation, this version does not take in non-linear state transition (f)
    and observation (h) functions. Instead, it assumes that the user provides the *state transition Jacobian*
    and *measurement Jacobian* directly, along with the respective covariance matrices.

    This class is intended for cases where the user has already linearized the system, or is working with
    systems where the state transition and observation models are already linear or approximately linear.

    Methods:
    --------
    predict(x_in, P_in, time_step, F, Q_c):
        Performs the prediction step of the EKF given the state transition Jacobian F and continuous-time
        process noise covariance Q_c.

    update(obs_vector, P_in, x_in, H, R):
        Performs the update step of the EKF using the measurement Jacobian H and measurement noise covariance R.
    """

    @staticmethod
    def predict(x_in, P_in, time_step, F, Q, continuous=False):
        """
        EKF prediction step using the provided linearized state transition matrix.

        Args:
            x_in (np.ndarray): Prior state estimate (n x 1).
            P_in (np.ndarray): Prior error covariance matrix (n x n).
            time_step (float): Time step between predictions.
            F (np.ndarray): State transition Jacobian matrix (n x n).
            Q (np.ndarray): Process noise covariance matrix (n x n).
            continuous(bool): if True then Q is continuous-time process noise, and is discretized
                if False, then Q is already the discrete-time covariance

        Returns:
            tuple[np.ndarray, np.ndarray]: tuple with predicted state estimate and predicted error covariance matrix.
        """
        if continuous:
            # discretization of Q
            Q_d = F @ Q @ F.T * time_step
        else:
            # Q_d is constructed manually for each state based on its noise model
            # e.g., Q_d[i, i] = sigma_i**2 * delta_t   (for random walk)
            # or    Q_d[i, i] = sigma_i**2 * (1 - alpha**2)  (for Gauss-Markov)
            Q_d = Q

        # state and covariance prediction
        x_out = F @ x_in
        P_out = F @ P_in @ F.T + Q_d

        return x_out, P_out

    @staticmethod
    def update(obs_vector, P_in, x_in, H, R):
        """
        EKF update step using the provided linearized measurement model.

        Args:
            obs_vector (np.ndarray): Innovation vector (measured - predicted observation), shape (m x 1).
            P_in (np.ndarray): Prior error covariance matrix (n x n).
            x_in (np.ndarray): Prior state estimate (n x 1).
            H (np.ndarray): Measurement Jacobian matrix (m x n).
            R (np.ndarray): Measurement noise covariance matrix (m x m).

        Returns:
            tuple [np.ndarray, np.ndarray]: tuple with updated state estimate and predicted error covariance matrix.
        """
        # kalman gain
        invS = np.linalg.inv(H @ P_in @ H.T + R)
        K = P_in @ H.T @ invS

        # update step
        P_out = P_in - K @ H @ P_in
        x_out = (x_in + K @ obs_vector)

        return x_out, P_out
