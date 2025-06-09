""" Module with the implementation of the Extended Kalman Filter solver """

import numpy as np


# TODO: add docstring
# NOTE: this EKF performs the predict step based on a linearized state transition function F
# (suitable for GNSS static/kinematic applications). The algorithm must be improved if a proper state function
# and jacobian are required
# the same goes for the observation function and its jacobian
#
class EKF:

    @staticmethod
    def predict(x_in, P_in, time_step, F, Q_c):
        # discretization of Q
        Q_d = F @ Q_c @ F.T * time_step

        # state and covariance prediction
        x_out = F @ x_in
        P_out = F @ P_in @ F.T + Q_d

        return x_out, P_out

    @staticmethod
    def update(obs_vector, P_in, x_in, H, R):
        # kalman gain
        invS = np.linalg.inv(H @ P_in @ H.T + R)
        K = P_in @ H.T @ invS

        # update step
        P_out = P_in - K @ H @ P_in
        x_out = (x_in + K @ obs_vector)

        return x_out, P_out
