import numpy as np
from src.models.frames import mechanization, attitude, frames
from src.utils.finite_diff import finite_difference
from src.utils.math_utils import rot1, rot2, rot3


# TODO: update docstrings

class SensorEmulator:
    def __init__(self, trace_path, data_manager):
        self.trace_path = trace_path
        self.data_manager = data_manager

    def __str__(self):
        return "InsAlgorithm(Sensor Emulation Algorithm)"

    def compute(self):
        # fetch PVAT inputs from the data manager
        ref_pos_lld = self.data_manager.get_data("ref_pos")
        ref_vel_eb_n = self.data_manager.get_data("ref_vel")
        ref_euler_att = self.data_manager.get_data("ref_att")
        time = np.array(ref_pos_lld.to_time_array())

        true_gyro, true_accel = self._compute_true_imu_readouts(time, ref_pos_lld, ref_vel_eb_n, ref_euler_att)


        # compute error readouts
        self._corrupt_noise_imu(true_gyro, true_accel, time)
        # self._corrupt_noise_gps(ref_pos_lld, ref_vel_eb_n)

    def _corrupt_noise_imu(self, true_gyro, true_accel, time):
        sampling_time = time[1, 0] - time[0, 0]

        readout_gyro = self._add_errors_imu("gyroscope", true_gyro, sampling_time)
        readout_accel = self._add_errors_imu("accelerometer", true_accel, sampling_time)

        self.data_manager.gyro.init_data(time, readout_gyro)
        self.data_manager.accel.init_data(time, readout_accel)

    def _add_misalignment(self, readouts, misalignment):
        for i in range(len(readouts)):
            M = rot1(misalignment[i, 0]) @ rot2(misalignment[i, 1]) @ rot3(misalignment[i, 2])
            readouts[i, :] = M @ readouts[i, :]

    def _add_errors_imu(self, sensor, true_read, sampling_time):
        user_models = self.data_manager.imu_sensor[sensor]
        n_epochs = len(true_read)

        misalignment = user_models["misalignment"].gen(n_epochs, sampling_time)
        scale_factor = user_models["scale_factor"].gen(n_epochs, sampling_time)
        bias_constant = user_models["bias_constant"].gen(n_epochs, sampling_time)
        bias_drift = user_models["bias_drift"].gen(n_epochs, sampling_time)
        noise = user_models["observation_noise"].gen(n_epochs, sampling_time)

        # attach misalignment contribution
        self._add_misalignment(true_read, misalignment)

        # attach all other contributions
        readout = (1 + scale_factor) * true_read + bias_drift + bias_constant + noise

        return readout

    def _corrupt_noise_gps(self, ref_pos_lld, ref_vel_ned):
        # compute gps_pos_ecef, gps_vel_ecef (ECEF format)
        # compute gps_pos_lld, gps_vel_ned (NED format)

        # compute noise vectors
        pos_noise = self.gps.get_stochastic_process(len(ref_pos_lld), "position").compute()
        vel_noise = self.gps.get_stochastic_process(len(ref_vel_ned), "velocity").compute()

        # compute reference ecef position
        ref_pos_ecef = lld2ecef(ref_pos_lld)
        gps_pos_ecef = np.zeros(np.shape(ref_pos_ecef))

        # add noise to NED velocity
        gps_vel_ned = ref_vel_ned + vel_noise
        gps_vel_ecef = np.zeros(np.shape(gps_vel_ned))

        # epoch loop
        for i in range(len(ref_pos_lld)):
            # get rotation matrix C_e^n (ECEF to NED)
            c_e_n = matrix_ecef2ned(ref_pos_lld[i, 0], ref_pos_lld[i, 1])

            # add noise to gps ecef position
            gps_pos_ecef[i, :] = ref_pos_ecef[i, :] + c_e_n.T @ pos_noise[i, :]

            # convert ned velocity to ecef
            gps_vel_ecef[i, :] = c_e_n.T @ gps_vel_ned[i, :]

        # transform gps ecef position to LLD form
        gps_pos_lld = ecef2lld(gps_pos_ecef)

        # ecef output
        gps_ecef = np.concatenate((gps_pos_ecef, gps_vel_ecef), axis=1)
        self.results.append(gps_ecef)

        # ned output
        gps_ned = np.concatenate((gps_pos_lld, gps_vel_ned), axis=1)
        self.results.append(gps_ned)

    def _compute_true_imu_readouts(self, time, pos_lld, vel_eb_n, euler_att):
        # n-frame mechanization

        # initialize outputs
        w_ib_b = np.zeros((len(pos_lld.data), 3))  # gyro is w_ib_b
        f_ib_b = np.zeros((len(pos_lld.data), 3))  # accel is f_ib_b

        # iterate in time
        for i, pos_row in pos_lld.data.iloc[1:].iterrows():
            # 1 - unpack data for this epoch
            t = time[i,0]
            step = t - time[i-1,0]
            lld = np.array(pos_row[1:4])
            v_eb_n = np.array(vel_eb_n.data.iloc[i,1:4])
            v_eb_n_prev = np.array(vel_eb_n.data.iloc[i-1, 1:4])
            euler = np.array(euler_att.data.iloc[i,1:4])
            euler_prev = np.array(euler_att.data.iloc[i-1, 1:4])

            # 2 - compute rotation matrices
            c_nb = attitude.euler2dcm(euler)  # matrix from n to b
            c_en = frames.latlon2dcm_e_ned(lld[0], lld[1])  # matrix from e to n

            # 3 - apply finite differences to euler angles
            euler_dot = finite_difference(euler_prev, euler, step)

            # 4 - apply finite differences to velocity vector (in n frame!)
            v_eb_n_dot = finite_difference(v_eb_n_prev, v_eb_n, step)

            # 5 - get rotation vectors between frames i-e, e-n and n-b, and local gravity vector
            w_nb_b = mechanization.compute_w_nb_b(euler_dot, euler)
            w_ie_n = mechanization.compute_w_ie_n(lld[0])
            w_en_n, _, _ = mechanization.compute_w_en_n(v_eb_n, lld)
            r_eb_e = np.array(frames.geodetic2cartesian(lld[0], lld[1], -lld[2])) # remember that altitude is positive downwards!
            g_eb_e = frames.grav_acceleration(r_eb_e)

            # 6 - finally, apply the mechanization equations (in the n-frame)
            w_ib_b[i] = mechanization.compute_w_ib_b(c_nb, w_ie_n, w_en_n, w_nb_b)
            f_ib_b[i] = mechanization.compute_f_ib_b(c_nb, c_en, w_en_n, w_ie_n, v_eb_n_dot, g_eb_e, v_eb_n)

        # append to results
        self.data_manager.ref_gyro.init_data(time, w_ib_b)
        self.data_manager.ref_accel.init_data(time, f_ib_b)

        return w_ib_b, f_ib_b