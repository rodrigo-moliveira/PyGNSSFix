""" Data Manager for GNSS Algorithms
"""
import os

import logging

from src import WORKSPACE_PATH
from src.data_mng.csv.csv_data import CSVData
from src.data_mng.unit_conversions import convert_unit
from src.io.states import OUTPUT_FILENAME_MAP, get_file_header, export_to_file
from src.io.config import config_dict, EnumAlgorithmPNT, EnumPositionForm, EnumAlgorithmINS, EnumVelocityFrame, EnumAttitudeForm
from src.io.rinex_parser import RinexNavReader, RinexObsReader, AntexReader
from src.models.frames import lla2lld
from src.models.noise.noise_manager import GNSSNoiseManager
from src.data_mng import Container
from src.common_log import IO_LOG, get_logger

__all__ = ["InsDataManager"]


class InsDataManager(Container):
    """
    TODO: to update
    Data Manager class that stores all necessary data for a GNSS run. Derives from the :py:class:`Container`
    class.

    Attributes:


    """
    __slots__ = [
        "ref_pos",
        "ref_pos_ecef",
        "ref_vel",
        "ref_att",
        "ref_gyro",
        "ref_accel",
        "pos",
        "vel",
        "att",
        "gyro",
        "accel",
        "gps"]

    def __init__(self):
        """ Default constructor with no arguments """
        super().__init__()

        ########################
        #    PVAT Variables    #
        ########################

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

        self.ref_vel = CSVData(name="ref_vel",
            description="Reference velocity in the navigation frame (NED)",
            units=['m/s', 'm/s', 'm/s'],
            legend=['vel_n', 'vel_e', 'vel_d'],
            title="Velocity v_eb^n (NED)",
            time_cols=[0],
            data_cols=[1, 2, 3])

        self.ref_att = CSVData(name="ref_att",
            description="Reference attitude (Euler angles, ZYX convention)",
            units=['rad', 'rad', 'rad'],
            legend=['ref_roll', 'ref_pitch', 'ref_yaw'],
            title="Attitude (Euler Angles)",
            time_cols=[0],
            data_cols=[1, 2, 3])

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
        log = get_logger(IO_LOG)

        # Import Position Data
        try:
            file_ref_pos = config_dict.get("inputs", "reference_position", "file")
            pos_form = EnumPositionForm.init_model(config_dict.get("inputs", "reference_position", "form"))
            pos_units = config_dict.get("inputs", "reference_position", "units")

            log.info(f"reading data '{self.ref_pos.name}' stored in file {file_ref_pos}...")
            self.ref_pos.read_data(f"{WORKSPACE_PATH}/{file_ref_pos}")

            # transform LLA or ECEF forms to LLD
            if pos_form == EnumPositionForm.POS_LLA:
                position = convert_unit(self.ref_pos.to_data_array(), pos_units, self.ref_pos.units)
                position = lla2lld(position)
            # TODO: add here the ref_pos_ecef!!!
            # elif pos_form == EnumPositionForm.POS_ECEF:
            #    position = convert_unit(self.ref_pos, pos_units, ("m", "m", "m"))
            #    position = ecef2lld(position)
            elif pos_form == EnumPositionForm.POS_LLD:  # position is already in LLD form
                position = convert_unit(self.ref_pos.to_data_array(), pos_units, self.ref_pos.units)
            else:
                raise AttributeError(f"Position form is not valid. Should be one of "
                                     f"{EnumPositionForm.show_options()}.")
            self.ref_pos.update_data_array(position)
        except Exception as e:
            log.warning(f"Did not successfully import data '{self.ref_pos.name}' due to: {repr(e)}")

        # Import Velocity Data
        try:
            file_ref_vel = config_dict.get("inputs", "reference_velocity", "file")
            vel_frame = EnumVelocityFrame.init_model(config_dict.get("inputs", "reference_velocity", "frame"))
            vel_units = config_dict.get("inputs", "reference_velocity", "units")
            log.info(f"reading data '{self.ref_vel.name}' stored in file {file_ref_vel}...")
            self.ref_vel.read_data(f"{WORKSPACE_PATH}/{file_ref_vel}")

            if vel_frame == EnumVelocityFrame.FRAME_N:
                velocity = convert_unit(self.ref_vel.to_data_array(), vel_units, self.ref_vel.units)
            # TODO: add here imports for FRAME_B and FRAME_E
            else:
                raise AttributeError(
                    f"Velocity frame is currently not available. Only navigation frame N is available.")

            self.ref_vel.update_data_array(velocity)

        except Exception as e:
            log.warning(f"Did not successfully import data '{self.ref_vel.name}' due to: {repr(e)}")

        # Import Attitude Data
        try:
            file_ref_att = config_dict.get("inputs", "reference_attitude", "file")
            att_form = EnumAttitudeForm.init_model(config_dict.get("inputs", "reference_attitude", "form"))
            att_units = config_dict.get("inputs", "reference_attitude", "units")
            log.info(f"reading data '{self.ref_att.name}' stored in file {file_ref_att}...")
            self.ref_att.read_data(f"{WORKSPACE_PATH}/{file_ref_att}")

            if att_form == EnumAttitudeForm.RPY:
                attitude = convert_unit(self.ref_att.to_data_array(), att_units, self.ref_att.units)
            # TODO: add here imports for QUATERNION
            else:
                raise AttributeError(
                    f"Attitude form is currently not available. Only form available is RPY.")

            self.ref_att.update_data_array(attitude)

        except Exception as e:
            log.warning(f"Did not successfully import data '{self.ref_vel.name}' due to: {repr(e)}")

    def read_inputs(self, ins_alg: EnumAlgorithmINS, trace_dir):
        """
        Function to read the input data. The data files are read from the user configurations
        according to the chosen algorithm.
        There are mandatory inputs for each GNSS Algorithm (SPS, PR-PPP).

        Args:
            gnss_alg(EnumAlgorithmPNT): GNSS Algorithm Enumeration
            trace_dir(str): path to write trace files
        Raises:
            IOError: an exception is raised if the trace directory is not created successfully
        """
        log = get_logger(IO_LOG)

        # read observation data
        # TODO: depending on algorithm selected, the inputs to be read are different
        if ins_alg == EnumAlgorithmINS.SENSOR_EMULATOR:
            log.info(f"Reading Sensor Emulator Inputs.")
            self._read_input_pva()

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
            with open(f"{inputs_dir}\\InputPosition.txt", "w") as file:
                file.write(str(self.get_data("ref_pos")))

            with open(f"{inputs_dir}\\InputVelocity.txt", "w") as file:
                file.write(str(self.get_data("ref_vel")))

            with open(f"{inputs_dir}\\InputAttitude.txt", "w") as file:
                file.write(str(self.get_data("ref_att")))


    @staticmethod
    def check_input_list(folder_name: str, input_list: list, log, warning=False):
        if len(input_list) == 0:
            if warning:
                log.warning(f"Optional input files {folder_name} not provided. "
                            f"Please check if they are necessary for some of the selected models.")
            else:
                raise IOError(f"Mandatory input files {folder_name} not provided. Please check the configurations.")
        if len(input_list) == 1 and not input_list[0]:
            if warning:
                log.warning(f"Optional input files {folder_name} not provided. "
                            f"Please check if they are necessary for some of the selected models.")
            else:
                raise IOError(f"Mandatory input files {folder_name} not provided. Please check the configurations.")
