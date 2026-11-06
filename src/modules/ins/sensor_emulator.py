""" Module with the Sensor Emulation Algorithm"""
import numpy as np
import pandas as pd
import os

from src.common_log import get_logger, INS_ALG_LOG
from src.errors import FiniteDifferenceError, SensorEmulationError
from src.models.frames import mechanization, base_convertion, sc_euler2dcm, rel_att_quat, rel_att_dcm
from src.utils.finite_diff import NumericalDifferentiator
from src.io.config import config_dict, EnumAttitudeDiffState, EnumAttitudeDiffMethod
from src.utils.math_utils import rot1, rot2, rot3


class SensorEmulator:
    """Emulates sensor measurements (gyroscope and accelerometer) for algorithm validation.

    This class replays or emulates sensor traces to feed INS/navigation algorithms
    with controlled, realistic measurement data.  It handles attitude-differencing
    configuration (state representation, method, finite-difference scheme) and
    velocity differentiation, and optionally writes intermediate trace files to disk.

    Sensor noise model
    ------------------
    Each gyro / accelerometer channel is modelled as::

        y_meas = (1 + scale_factor) * M @ y_true + bias_drift + bias_constant + noise

    where

    * **scale_factor** – per-axis scale-factor error (dimensionless).
    * **M** – misalignment matrix (small-angle tilt between sensor and body axes).
    * **y_true** – ideal (bias-free, noise-free) sensor reading.
    * **bias_drift** – slowly time-varying bias (e.g. first-order Markov / random walk).
    * **bias_constant** – fixed bias offset.
    * **noise** – zero-mean white measurement noise (typically white Gaussian).

    Attributes:
        trace_dir: Absolute path to the directory where sensor trace files are
            written, or ``None`` if trace-file output is disabled.
        data_manager: Reference to the central data manager for trace I/O.
        log: Logger instance scoped to the INS algorithm log channel.
        att_diff_state: Enum describing the attitude representation used for
            differencing (e.g. RPY, QUAT, DCM).
        att_diff_method: Enum describing how the attitude difference is computed
            (e.g. FIN_DEF for finite-difference).
        att_differentiator: A :class:`NumericalDifferentiator` instance for attitude,
            or ``None`` when a relative attitude difference is used.
        vel_differentiator: A :class:`NumericalDifferentiator` instance for velocity.
    """
    def __init__(self, trace_path, data_manager):
        """Initialise the sensor emulator from the global configuration.

        Reads the ``sensor_emulator`` section of *config_dict* to determine the
        attitude-differencing state, method, and finite-difference schemes, then
        sets up the required differentiators.  If the *algorithm* section requests
        trace-file output, a dedicated sub-directory is created under *trace_path*.

        Args:
            trace_path: Base directory under which sensor trace files will be
                stored
            data_manager: Instance responsible for reading/writing trace data.

        Raises:
            IOError: If the trace directory cannot be created.
            SensorEmulationError: If the attitude difference method is set to
                ``FIN_DEF`` but the attitude state is ``RPY``, which is an
                unsupported combination.

        """
        self.trace_dir = None
        self.data_manager = data_manager
        self.log = get_logger(INS_ALG_LOG)

        if config_dict.get("algorithm", "trace_files"):
            self.trace_dir = f"{trace_path}\\sensor_emulator"
            try:
                os.makedirs(self.trace_dir)
            except:
                raise IOError(f"Cannot create dir: {self.trace_dir}")

        # configuration from the user
        self.att_diff_state = EnumAttitudeDiffState.init_model(config_dict.get("sensor_emulator", "attitude", "diff_state"))
        self.att_diff_method = EnumAttitudeDiffMethod.init_model(config_dict.get("sensor_emulator", "attitude", "method"))
        if self.att_diff_method == EnumAttitudeDiffMethod.FIN_DEF:
            self.att_differentiator = NumericalDifferentiator(config_dict.get("sensor_emulator", "attitude", "finite_diff_method"))
        else:
            self.att_differentiator = None
            if self.att_diff_state == EnumAttitudeDiffState.RPY:
                raise SensorEmulationError(f"When Relative Attitude is defined, the attitude state cannot be set to RPY. "
                                           f"Please revise the attitude state to either QUAT or DCM.")
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
        """Generate the full set of (noisy) IMU sensor readouts for the current run.

        This is the main entry point of the :class:`SensorEmulator`.  It orchestrates
        a two-stage pipeline:

        1. **True readouts** – Compute bias-free, noise-free gyroscope and
           accelerometer outputs from the reference trajectory.
        2. **Error injection** – Corrupt the true readouts with the configured
           misalignment, scale-factor, bias, and noise models (see the class-level
           *Sensor noise model* section).

        The resulting corrupted data is written into
        ``self.data_manager.gyro`` and ``self.data_manager.accel`` for downstream
        algorithm consumption.
        """
        # compute errorless readouts
        true_gyro, true_accel, time = self._compute_true_imu_readouts()

        # compute error readouts
        self._corrupt_noise_imu(true_gyro, true_accel, time)
        # self._corrupt_noise_gps(ref_pos_lld, ref_vel_eb_n)

    def _corrupt_noise_imu(self, true_gyro, true_accel, time):
        """Apply the full IMU error model and store corrupted readouts.

        Delegates to :meth:`_add_errors_imu` for each sensor (gyroscope and
        accelerometer) and writes the resulting arrays into the data manager so
        that downstream algorithms read the noisy values.

        Args:
            true_gyro: Ideal gyroscope readings
            true_accel: Ideal accelerometer readings
            time: Time-stamp array
        """
        sampling_time = time[1, 0] - time[0, 0]

        readout_gyro = self._add_errors_imu("gyroscope", true_gyro, sampling_time)
        readout_accel = self._add_errors_imu("accelerometer", true_accel, sampling_time)

        self.data_manager.gyro.init_data(time, readout_gyro)
        self.data_manager.accel.init_data(time, readout_accel)

    def _add_misalignment(self, readouts, misalignment):
        """Apply per-epoch misalignment rotations in place.

        For each epoch, a 3-2-1 (Z-Y-X) intrinsic rotation matrix is constructed
        from the three Euler angles supplied in *misalignment* and multiplied
        onto the corresponding readout row.

        Args:
            readouts: Sensor reading matrix, shape ``(n_epochs, 3)``.  Modified
                in place.
            misalignment: Per-epoch Euler angles, shape ``(n_epochs, 3)``, where
                the columns correspond to rotation about axis 1, axis 2, and
                axis 3 (radians).
        """
        for i in range(len(readouts)):
            M = rot1(misalignment[i, 0]) @ rot2(misalignment[i, 1]) @ rot3(misalignment[i, 2])
            readouts[i, :] = M @ readouts[i, :]

    def _add_errors_imu(self, sensor, true_read, sampling_time):
        """Generate a corrupted sensor readout by stacking all error terms.

        Retrieves the user-configured error models for *sensor* from
        ``self.data_manager.imu_sensor``, generates samples for each model, and
        combines them according to::

            readout = (1 + scale_factor) * (M @ true_read)
                      + bias_drift + bias_constant + noise

        where ``M`` is the per-epoch misalignment rotation.

        Args:
            sensor: Sensor identifier – ``"gyroscope"`` or ``"accelerometer"``.
            true_read: Ideal sensor readings, shape ``(n_epochs, 3)``.  The
                misalignment step modifies this array in place; callers that need
                the original values should pass a copy.
            sampling_time: Time interval between consecutive samples in seconds,
                used by each error model for temporal correlation (e.g. Allan
                variance, random walk).

        Returns:
            Corrupted readout matrix, shape ``(n_epochs, 3)``, with the same
            physical units as *true_read*.
        """
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

    """def _corrupt_noise_gps(self, ref_pos_lld, ref_vel_ned):
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
        self.results.append(gps_ned)"""

    def _compute_true_imu_readouts(self):
        """Compute bias-free, noise-free IMU readouts via N-frame mechanization.

        Reads the reference trajectory in the navigation (N) frame and
        extracts the corresponding ideal gyroscope and accelerometer outputs
        that an IMU mounted on the vehicle would observe.

        The mechanization equations are formulated in the N frame as described
        in [1]_.

            [1] Springer Handbook of Global Navigation Satellite Systems,
            Peter J.G. Teunissen, Oliver Montenbruck, Springer Cham, 2017.
            Chapter 28

        Returns:
            A tuple of three arrays:

            * ``true_gyro`` – Ideal gyroscope body-rate readings,
              shape ``(n_epochs, 3)``, in rad/s.
            * ``true_accel`` – Ideal accelerometer specific-force readings,
              shape ``(n_epochs, 3)``, in m/s².
            * ``time`` – Time-stamp matrix, shape ``(n_epochs, 2)``, where
              column 0 is the epoch time (s) and column 1 is the sampling
              offset.
        """

        if self.trace_dir is not None:
            fh_omega = open(f"{self.trace_dir}\\w_nb_b.txt", "w")
            fh_force = open(f"{self.trace_dir}\\v_eb_n_dot.txt", "w")

            fh_omega.write("epoch,w_ib^b X [rad/s],w_ib^b Y [rad/s],w_ib^b Z [rad/s]\n")
            fh_force.write("epoch,f_ib^b X [m/s^2],f_ib^b Y [m/s^2],f_ib^b Z [m/s^2]\n")
        else:
            fh_omega = fh_force = None

        # fetch PVAT inputs from the data manager
        ref_pos_lld = self.data_manager.get_data("ref_pos")
        ref_vel_eb_n = self.data_manager.get_data("ref_vel_eb_n")
        ref_euler_att = self.data_manager.get_data("ref_att_euler")
        ref_quat_att = self.data_manager.get_data("ref_att_quat")
        time = np.array(ref_pos_lld.to_time_array())
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
                    # NOTE: we may need to check for each case, how to improve the correctness of the finite difference (for example, call quat_normalize)
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

        # Append to results to data manager
        self.data_manager.true_gyro.init_data(time, w_ib_b)
        self.data_manager.true_accel.init_data(time, f_ib_b)

        # Close trace files
        if fh_omega is not None and fh_force is not None:
            fh_omega.close()
            fh_force.close()

        return w_ib_b, f_ib_b, time