import numpy as np
import pandas as pd
import os

from src.common_log import get_logger, INS_ALG_LOG
from src.errors import FiniteDifferenceError, SensorEmulationError
from src.models.frames import mechanization, base_convertion, sc_euler2dcm, rel_att_quat, rel_att_dcm
from src.utils.finite_diff import NumericalDifferentiator
from src.io.config import config_dict, EnumAttitudeDiffState, EnumAttitudeDiffMethod
from src.utils.math_utils import rot1, rot2, rot3


# TODO: update docstrings

class SensorEmulator:
    def __init__(self, trace_path, data_manager):
        self.trace_dir = None
        self.data_manager = data_manager
        self.log = get_logger(INS_ALG_LOG)

        if config_dict.get("algorithm", "trace_files"):
            self.trace_dir = f"{trace_path}\\sensor_emulator"
            try:
                os.makedirs(self.trace_dir)
            except:
                raise IOError(f"Cannot create dir: {self.trace_dir}")

        # attitude configuration
        self.att_diff_state = EnumAttitudeDiffState.init_model(config_dict.get("sensor_emulator", "attitude", "diff_state"))
        self.att_diff_method = EnumAttitudeDiffMethod.init_model(config_dict.get("sensor_emulator", "attitude", "method"))
        if self.att_diff_method == EnumAttitudeDiffMethod.FIN_DEF:
            self.att_differentiator = NumericalDifferentiator(config_dict.get("sensor_emulator", "attitude", "finite_diff_method"))
        else:
            self.att_differentiator = None
            if self.att_diff_state == EnumAttitudeDiffState.RPY:
                raise SensorEmulationError(f"When Relative Attitude is defined, the attitude state cannot be set to RPY. "
                                           f"Please revise the attitude state to either QUAT or DCM.")

        # velocity configuration
        self.vel_differentiator = NumericalDifferentiator(config_dict.get("sensor_emulator", "velocity", "finite_diff_method"))

        # Printing to log
        self.log.info(f"Sensor Emulator initialized with the following parameters:")
        self.log.info(f"Attitude Difference State: {self.att_diff_state}")
        self.log.info(f"Attitude Differentiation Method: {self.att_diff_method}")
        if self.att_diff_method == EnumAttitudeDiffMethod.FIN_DEF:
            self.log.info(f"Attitude Finite Difference Method: {config_dict.get("sensor_emulator", "attitude", "finite_diff_method")}")
        self.log.info(f"Velocity Finite Difference Method: {config_dict.get("sensor_emulator", "velocity", "finite_diff_method")}")


    def __str__(self):
        return "InsAlgorithm(Sensor Emulation Algorithm)"

    def compute(self):
        # compute errorless readouts
        true_gyro, true_accel, time = self._compute_true_imu_readouts()

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

    def _compute_true_imu_readouts(self):
        # n-frame mechanization

        if self.trace_dir is not None:
            fh_omega = open(f"{self.trace_dir}\\w_nb_b.txt", "w")
            fh_force = open(f"{self.trace_dir}\\v_eb_n_dot.txt", "w")

            fh_omega.write("epoch,w_ib_b X [rad/s],w_ib_b Y [rad/s],w_ib_b Z [rad/s]\n")
            fh_force.write("epoch,f_ib_b X [m/s^2],f_ib_b Y [m/s^2],f_ib_b Z [m/s^2]\n")
        else:
            fh_omega = fh_force = None

        # fetch PVAT inputs from the data manager
        ref_pos_lld = self.data_manager.get_data("ref_pos")
        ref_vel_eb_n = self.data_manager.get_data("ref_vel_eb_n")
        ref_euler_att = self.data_manager.get_data("ref_att_euler")
        ref_quat_att = self.data_manager.get_data("ref_att_quat")
        time = np.array(ref_pos_lld.to_time_array())
        # TODO: consider if ref_dcm should be a state variable
        ref_dcm_att = base_convertion(ref_euler_att.to_data_array(), sc_euler2dcm)
        ref_dcm_att = pd.DataFrame(np.concatenate((time, ref_dcm_att), axis=1))

        # initialize outputs
        w_ib_b = np.ones((len(ref_pos_lld.data), 3))  * np.nan # gyro is w_ib_b
        f_ib_b = np.ones((len(ref_pos_lld.data), 3)) * np.nan  # accel is f_ib_b

        # iterate in time
        for i, pos_row in ref_pos_lld.data.iloc[:].iterrows():
            # 1 - unpack data for this epoch
            t = time[i, 0]
            lld = np.array(pos_row[1:4])
            v_eb_n = np.array(ref_vel_eb_n.data.iloc[i, 1:4])
            euler = np.array(ref_euler_att.data.iloc[i, 1:4])
            quat = np.array(ref_quat_att.data.iloc[i, 1:])

            # 2 - compute rotation matrices
            # c_nb = attitude.euler2dcm(euler)  # matrix from n to b
            c_nb = np.array(ref_dcm_att.iloc[i, 1:]).reshape(3, 3)

            # 3. compute attitude difference (w_nb_b)
            w_nb_b = np.zeros(3)
            if self.att_differentiator is not None:
                # finite difference method
                try:
                    # TODO: check for each case, how to improve the correctness of the finite difference (for example, call quat_normalize)
                    if self.att_diff_state == EnumAttitudeDiffState.RPY:
                        euler_dot = self.att_differentiator.compute(ref_euler_att.data, 0, ref_euler_att.data_cols, t)
                        w_nb_b = mechanization.compute_w_nb_b_euler(euler_dot, euler)
                    elif self.att_diff_state == EnumAttitudeDiffState.QUAT:
                        quat_dot = self.att_differentiator.compute(ref_quat_att.data, 0, ref_quat_att.data_cols, t)
                        w_nb_b = mechanization.compute_w_nb_b_quat(quat_dot, quat)
                    elif self.att_diff_state == EnumAttitudeDiffState.DCM:
                        dcm_dot = self.att_differentiator.compute(ref_dcm_att, 0, list(range(1,10)), t)
                        w_nb_b = mechanization.compute_w_nb_b_dcm(np.array(dcm_dot).reshape(3,3), c_nb)
                except FiniteDifferenceError as e:
                    self.log.warning(f"FiniteDifferenceError caught when computing attitude finite difference due to: {e}")
                    continue
            else:
                try:
                    # relative attitude
                    dt = time[i+1, 0] - t
                    if self.att_diff_state == EnumAttitudeDiffState.QUAT:
                        quat_next = np.array(ref_quat_att.data.iloc[i+1, 1:])
                        w_nb_b = rel_att_quat(quat, quat_next, dt)
                    elif self.att_diff_state == EnumAttitudeDiffState.DCM:
                        dcm_next = np.array(ref_dcm_att.iloc[i + 1, 1:]).reshape(3, 3)
                        w_nb_b = rel_att_dcm(c_nb, dcm_next, dt)
                except Exception as e:
                    self.log.warning(f"Exception caught when computing relative attitude difference due to: {e}")
                    continue

            # 4. compute velocity finite difference (v_eb_n_dot)
            try:
                v_eb_n_dot = self.vel_differentiator.compute(ref_vel_eb_n.data, 0, ref_vel_eb_n.data_cols, t)
            except FiniteDifferenceError as e:
                self.log.warning(f"FiniteDifferenceError caught when computing velocity finite difference due to: {e}")
                continue

            # 5. Apply the mechanization equations (in the n-frame)
            _w_ib_b, _f_ib_b = mechanization.compute_imu_readouts(lld, w_nb_b, v_eb_n, v_eb_n_dot, c_nb)
            w_ib_b[i] = _w_ib_b
            f_ib_b[i] = _f_ib_b

            # 6 - save to trace file
            if fh_omega is not None:
                fh_omega.write(f"{t},{w_nb_b[0]},{w_nb_b[1]},{w_nb_b[2]}\n")
            if fh_force is not None:
                fh_force.write(f"{t},{v_eb_n_dot[0]},{v_eb_n_dot[1]},{v_eb_n_dot[2]}\n")

        # append to results
        self.data_manager.true_gyro.init_data(time, w_ib_b)
        self.data_manager.true_accel.init_data(time, f_ib_b)


        if fh_omega is not None and fh_force is not None:
            fh_omega.close()
            fh_force.close()

        return w_ib_b, f_ib_b, time