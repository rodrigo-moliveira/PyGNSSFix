""" Data Manager for INS Algorithms
"""
import os

from src import WORKSPACE_PATH
from src.data_mng.csv.csv_data import CSVData
from src.data_mng.unit_conversions import convert_unit
from src.io.config import config_dict, EnumPositionForm, EnumAlgorithmINS, EnumVelocityFrame, EnumAttitudeForm
from src.io.config.import_imu import read_imu_file
from src.models.frames import *
from src.data_mng import Container
from src.common_log import IO_LOG, get_logger

__all__ = ["InsDataManager"]

from src.models.sensor.imu import IMU


class InsDataManager(Container):
    """Manage reference, navigation, and sensor data for an INS/GNSS run.

    Stores the data required for an Inertial Navigation System (INS) and
    Global Navigation Satellite System (GNSS) processing run. The class
    derives from :py:class:`Container` and organizes reference position,
    velocity, and attitude data; reference and measured IMU observations;
    GNSS data; and IMU sensor configuration.

    Most attributes are instances of :py:class:`CSVData`, which provide
    structured storage and metadata for time-tagged data. The ``gps`` and
    ``imu_sensor`` attributes hold their respective data or configuration
    objects.

    Attributes:
        ref_pos (CSVData): Reference geodetic position in latitude,
            longitude, and down coordinates (LLD). Columns represent
            latitude [rad], longitude [rad], and down [m].
        ref_pos_ecef (CSVData): Reference position in Earth-Centered,
            Earth-Fixed (ECEF) Cartesian coordinates.
        ref_vel_eb_n (CSVData): Reference velocity of the body with respect
            to the Earth, expressed in the local North-East-Down (NED) frame.
        ref_vel_eb_e (CSVData): Reference velocity of the body with respect
            to the Earth, expressed in the ECEF frame.
        ref_att_euler (CSVData): Reference attitude represented by Euler
            angles, using the navigation-to-body rotation convention, RPY_n^b.
        ref_att_quat (CSVData): Reference attitude represented by
            quaternions in scalar-last order [qx, qy, qz, qw], q_n^b.
        ref_gyro (CSVData): Reference angular-rate data from the gyroscope, w_ib^b.
        ref_accel (CSVData): Reference acceleration data from the
            accelerometer, f_ib^b.
        gyro (CSVData): Gyroscope measurements or simulated observations.
        accel (CSVData): Accelerometer measurements or simulated observations.
        gps: GNSS data object with GNSS data-source configuration.
        imu_sensor: IMU sensor object with sensor configuration.
    """
    __slots__ = [
        "ref_pos",
        "ref_pos_ecef",
        "ref_vel_eb_n",
        "ref_vel_eb_e",
        "ref_att_euler",
        "ref_att_quat",
        "ref_gyro",
        "ref_accel",
        "gyro",
        "accel",
        "gps",
        "imu_sensor"]

    def __init__(self):
        """Initialize the data manager and its data containers.

        Calls the parent :py:class:`Container` constructor and initializes
        the attributes used to store reference navigation data, IMU data,
        GNSS data, and sensor configuration.
        """
        super().__init__()

        #######################
        #    PVA Variables    #
        #######################

        self.ref_pos = CSVData(name="ref_pos",
            description="Reference position in the navigation (NED) frame, in LLD form",
            units=['rad', 'rad', 'm'],
            legend=['pos_lat', 'pos_lon', 'pos_down'],
            title="Position (LLD)",
            time_cols=[0],
            data_cols=[1, 2, 3])

        self.ref_pos_ecef = CSVData(name="ref_pos_ecef",
            description="Reference position in the ECEF frame",
            units=['m', 'm', 'm'],
            legend=['pos_x', 'pos_y', 'pos_z'],
            title="Position (ECEF)",
            time_cols=[0],
            data_cols=[1, 2, 3])

        self.ref_vel_eb_n = CSVData(name="ref_vel_eb^n",
            description="Reference velocity in the navigation frame (NED)",
            units=['m/s', 'm/s', 'm/s'],
            legend=['vel_n', 'vel_e', 'vel_d'],
            title="Velocity v_eb^n (NED)",
            time_cols=[0],
            data_cols=[1, 2, 3])

        self.ref_vel_eb_e = CSVData(name="ref_vel_eb^e",
            description="Reference velocity in the Earth Fixed Frame (ECEF)",
            units=['m/s', 'm/s', 'm/s'],
            legend=['vel_x', 'vel_y', 'vel_z'],
            title="Velocity v_eb^e (ECEF)",
            time_cols=[0],
            data_cols=[1, 2, 3])

        self.ref_att_euler = CSVData(name="ref_att_euler",
            description="Reference attitude (Euler angles, ZYX convention) from navigation n to body frame b",
            units=['rad', 'rad', 'rad'],
            legend=['ref_roll', 'ref_pitch', 'ref_yaw'],
            title="Attitude (Euler Angles)",
            time_cols=[0],
            data_cols=[1, 2, 3])

        self.ref_att_quat = CSVData(name="ref_att_quat",
             description="Reference attitude quaternion (scalar first, Q0 convention) from navigation n to body frame b, q_n^b",
             units=['', '', '', ''],
             legend=['q0', 'q1', 'q2', 'q3'],
             title="Attitude quaternion q_n^b",
             time_cols=[0],
             data_cols=[1, 2, 3, 4])

        #######################
        # Sensor Measurements #
        #######################

        # True (errorless) gyro readout measurements
        self.ref_gyro = CSVData(name='ref_gyro',
            description='true angular velocity in the body frame (w_ib_b)',
            units=['rad/s', 'rad/s', 'rad/s'],
            legend=['gyro_x', 'gyro_y', 'gyro_z'],
            title="Gyroscope True w_ib_b",
            time_cols=[0],
            data_cols=[1, 2, 3])

        # True (errorless) accelerometer readout measurements
        self.ref_accel = CSVData(name='ref_accel',
            description='true acceleration in the body frame (f_ib_b)',
            units=['m/s^2', 'm/s^2', 'm/s^2'],
            legend=['accel_x', 'accel_y', 'accel_z'],
            title="Accelerometer True f_ib_b",
            time_cols=[0],
            data_cols=[1, 2, 3])

        # Real (with error) gyro readout measurements
        self.gyro = CSVData(name='gyro',
            description='gyro measurements w_ib_b',
            units=['rad/s', 'rad/s', 'rad/s'],
            legend=['gyro_x', 'gyro_y', 'gyro_z'],
            title="Gyroscope Readouts w_ib_b",
            time_cols=[0],
            data_cols=[1, 2, 3])

        # Real (with error) accelerometer readout measurements
        self.accel = CSVData(name='accel',
             description='accelerometer measurements f_ib_b',
             units=['m/s^2', 'm/s^2', 'm/s^2'],
             legend=['accel_x', 'accel_y', 'accel_z'],
             title="Accelerometer Readouts f_ib_b",
             time_cols=[0],
             data_cols=[1, 2, 3])

        # Sensors
        self.imu_sensor = IMU()
        self.gps = None

    def __str__(self):
        return f'{type(self).__name__}( DataManager for INS algorithms )'

    def __repr__(self):
        return str(self)

    def add_data(self, data_name, data):
        """
        Adds the provided data to this DataManager instance.

        Args:
            data_name(str): the string name of the data to fill, matching one of the class attributes
            data(Any): data object to be stored

        Raises:
            AttributeError: an `AttributeError` is raised if the `data_name` string does not match one of the
                class attributes
        """
        if data_name in self.__slots__:
            setattr(self, data_name, data)
        else:
            raise AttributeError(f"Unsupported data: {data_name}, not in {self.__slots__}")

    def get_data(self, data_name):
        """
        Args:
            data_name(str): the data to be fetched
        Returns:
            Any: returns the queried data for the provided `data_name` attribute
        """
        # single data
        if isinstance(data_name, str):
            if data_name in self.__slots__:
                return getattr(self, data_name)
            else:
                raise AttributeError(f'{data_name} is not available.')
        return None

    def _read_input_pva(self):
        """ Reads input Position, Velocity and Attitude data and fills the corresponding attributes """
        log = get_logger(IO_LOG)

        # Import Position Data
        file_ref_pos = config_dict.get("inputs", "reference_position", "file")
        pos_form = EnumPositionForm.init_model(config_dict.get("inputs", "reference_position", "form"))
        pos_units = config_dict.get("inputs", "reference_position", "units")

        log.info(f"reading position data stored in file {file_ref_pos} (pos_form={pos_form})...")

        # transform LLA or ECEF forms to LLD
        if pos_form == EnumPositionForm.POS_LLA or pos_form == EnumPositionForm.POS_LLD:
            # Input Position provided in LLA Form
            self.ref_pos.read_data(f"{WORKSPACE_PATH}/{file_ref_pos}")
            time = np.array(self.ref_pos.to_time_array())
            self.ref_pos_ecef.init_data(time, np.zeros((len(time), 3)))

            if pos_form == EnumPositionForm.POS_LLA:
                pos_lla = convert_unit(self.ref_pos.to_data_array(), pos_units, self.ref_pos.units)
                pos_lld = base_convertion(pos_lla, lla2lld)
            else:
                pos_lld = convert_unit(self.ref_pos.to_data_array(), pos_units, self.ref_pos.units)
            pos_ecef = base_convertion(pos_lld, lld2ecef)

        elif pos_form == EnumPositionForm.POS_ECEF:
            # Input Position provided in ECEF Form
            self.ref_pos_ecef.read_data(f"{WORKSPACE_PATH}/{file_ref_pos}")
            time = np.array(self.ref_pos_ecef.to_time_array())
            self.ref_pos.init_data(time, np.zeros((len(time), 3)))

            pos_ecef = convert_unit(self.ref_pos_ecef.to_data_array(), pos_units, self.ref_pos_ecef.units)
            pos_lld = base_convertion(pos_ecef, ecef2lld)

        else:
            raise AttributeError(f"Position form is not valid. Should be one of "
                                 f"{EnumPositionForm.show_options()}.")
        self.ref_pos.update_data_array(pos_lld)
        self.ref_pos_ecef.update_data_array(pos_ecef)

        # Import Attitude Data
        file_ref_att = config_dict.get("inputs", "reference_attitude", "file")
        att_form = EnumAttitudeForm.init_model(config_dict.get("inputs", "reference_attitude", "form"))
        att_units = config_dict.get("inputs", "reference_attitude", "units")

        log.info(f"reading attitude data stored in file {file_ref_att} (pos_form={att_form})...")

        if att_form == EnumAttitudeForm.RPY:
            self.ref_att_euler.read_data(f"{WORKSPACE_PATH}/{file_ref_att}")
            time = np.array(self.ref_att_euler.to_time_array())
            self.ref_att_quat.init_data(time, np.zeros((len(time), 4)))
            att_euler = convert_unit(self.ref_att_euler.to_data_array(), att_units, self.ref_att_euler.units)
            att_quat = base_convertion(att_euler, euler2quat)
        elif att_form == EnumAttitudeForm.QUAT:
            self.ref_att_quat.read_data(f"{WORKSPACE_PATH}/{file_ref_att}")
            time = np.array(self.ref_att_quat.to_time_array())
            self.ref_att_euler.init_data(time, np.zeros((len(time), 3)))
            att_quat = np.array(self.ref_att_quat.to_data_array())
            att_euler = base_convertion(att_quat, quat2euler)
        else:
            raise AttributeError(
                f"Attitude form is currently not available. Only form available is RPY.")

        self.ref_att_euler.update_data_array(att_euler)
        self.ref_att_quat.update_data_array(att_quat)

        # Import Velocity Data
        file_ref_vel = config_dict.get("inputs", "reference_velocity", "file")
        vel_frame = EnumVelocityFrame.init_model(config_dict.get("inputs", "reference_velocity", "frame"))
        vel_units = config_dict.get("inputs", "reference_velocity", "units")
        log.info(f"reading data '{self.ref_vel_eb_n.name}' stored in file {file_ref_vel}...")

        if vel_frame == EnumVelocityFrame.FRAME_N:
            self.ref_vel_eb_n.read_data(f"{WORKSPACE_PATH}/{file_ref_vel}")
            time = np.array(self.ref_vel_eb_n.to_time_array())
            self.ref_vel_eb_e.init_data(time, np.zeros((len(time), 3)))

            vel_eb_n = convert_unit(self.ref_vel_eb_n.to_data_array(), vel_units, self.ref_vel_eb_n.units)
            vel_eb_e = base_convertion(vel_eb_n, vel_n2e, np.array(self.ref_pos.to_data_array()))

        elif vel_frame == EnumVelocityFrame.FRAME_E:
            self.ref_vel_eb_e.read_data(f"{WORKSPACE_PATH}/{file_ref_vel}")
            time = np.array(self.ref_vel_eb_e.to_time_array())
            self.ref_vel_eb_n.init_data(time, np.zeros((len(time), 3)))

            vel_eb_e = convert_unit(self.ref_vel_eb_e.to_data_array(), vel_units, self.ref_vel_eb_e.units)
            vel_eb_n = base_convertion(vel_eb_e, vel_e2n, np.array(self.ref_pos.to_data_array()))
        else:
            raise AttributeError(
                f"Velocity frame is currently not available. Only navigation frame N and E is available.")
        self.ref_vel_eb_n.update_data_array(vel_eb_n)
        self.ref_vel_eb_e.update_data_array(vel_eb_e)


    def read_inputs(self, ins_alg: EnumAlgorithmINS, trace_dir):
        """
        Function to read the input data. The data files are read from the user configurations
        according to the chosen algorithm.
        There are mandatory inputs for each INS Algorithm.

        Args:
            ins_alg(EnumAlgorithmINS): INS Algorithm Enumeration
            trace_dir(str): path to write trace files
        Raises:
            IOError: an exception is raised if the trace directory is not created successfully
        """
        log = get_logger(IO_LOG)

        if ins_alg == EnumAlgorithmINS.SENSOR_EMULATOR:
            log.info(f"Reading Sensor Emulator Inputs.")

            # Input Reference PVA data
            self._read_input_pva()

            # Sensor (IMU/GNSS) configuration
            imu_path = config_dict.get("inputs", "imu_sensor_file")
            read_imu_file(f"{WORKSPACE_PATH}/{imu_path}", self.imu_sensor, log=log)
            # Add GNSS file here

        else:
            raise IOError(f"Unknown Model {ins_alg}")

        if config_dict.get("inputs", "trace_files"):
            self._trace_files(trace_dir, ins_alg)

    def _trace_files(self, trace_dir, ins_alg: EnumAlgorithmINS):
        inputs_dir = f"{trace_dir}\\inputs"
        try:
            os.makedirs(inputs_dir)
        except:
            raise IOError(f"Cannot create dir: {inputs_dir}")

        # trace data files
        if ins_alg == EnumAlgorithmINS.SENSOR_EMULATOR:
            with open(f"{inputs_dir}\\InputPositionLLD.txt", "w", newline="") as file:
                file.write(self.ref_pos.to_file())

            with open(f"{inputs_dir}\\InputPositionECEF.txt", "w", newline="") as file:
                file.write(self.ref_pos_ecef.to_file())

            with open(f"{inputs_dir}\\InputVelocityNED.txt", "w", newline="") as file:
                file.write(self.ref_vel_eb_n.to_file())

            with open(f"{inputs_dir}\\InputVelocityECEF.txt", "w", newline="") as file:
                file.write(self.ref_vel_eb_e.to_file())

            with open(f"{inputs_dir}\\InputAttitudeEuler.txt", "w", newline="") as file:
                file.write(self.ref_att_euler.to_file())

            with open(f"{inputs_dir}\\InputAttitudeQuat.txt", "w", newline="") as file:
                file.write(self.ref_att_quat.to_file())

    def save_data(self, directory):
        log = get_logger(IO_LOG)

        for attribute in self.__slots__:
            data = getattr(self, attribute, None)

            if isinstance(data, CSVData):
                if not data.is_empty():
                    if "ref" in data.name:
                        # Skip (input) reference data (these should go to trace files)
                        continue
                    log.info(f"Saving data {data.name} to {directory}")
                    with open(f"{directory}\\{data.name}.txt", "w", newline="") as file:
                        file.write(data.to_file())
