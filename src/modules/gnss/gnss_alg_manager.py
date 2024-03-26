import os
import time
import traceback
import numpy as np

from src import RUNS_PATH
from src.data_mng.gnss import GnssDataManager
from src.io.config import config_dict, EnumPositioningMode
from src.models.frames import cartesian2geodetic, latlon2dcm_e_enu
from src.errors import ConfigError
from src.common_log import MAIN_LOG, get_logger

from .solver.gnss_solver import GnssSolver
from .preprocessor import PreprocessorManager


class GnssAlgorithmManager:

    def __init__(self):
        # create output folder
        data_dir = config_dict.get("output", "output_path")
        self.data_dir = self._check_data_dir(data_dir)

        # create data members
        self.data_manager = GnssDataManager()

    def run(self):
        MAIN_LOGGER = get_logger(MAIN_LOG)
        MAIN_LOGGER.info("Starting GNSS Algorithm Manager")
        model = config_dict.get('model', 'mode')

        if model not in (EnumPositioningMode.SPS, EnumPositioningMode.SPS_IF):
            raise ConfigError(f"Selected Model {model} not valid. Available options are "
                              f"SPS, SPS_IF")
        MAIN_LOGGER.info(f"Running GNSS algorithm {model}")

        # Input Reader Module
        try:
            MAIN_LOGGER.info(f"Starting Input Reader Module...")
            self.data_manager.read_inputs(f"{self.data_dir}\\trace")
        except Exception as e:
            MAIN_LOGGER.error(f"Stopping execution of program due to error in execution of Input Reader Module: {e}")
            print(traceback.format_exc())
            exit(-1)

        # Main Algorithm Module
        try:
            MAIN_LOGGER.info(f"Starting Main Algorithm Module...")
            self.compute(self.data_manager, f"{self.data_dir}\\trace")
        except Exception as e:
            MAIN_LOGGER.error(f"Stopping execution of program due to error in execution of Main Algorithm Module: {e}")
            print(traceback.format_exc())
            exit(-1)

        # Output Writer Module
        try:
            MAIN_LOGGER.info(f"Starting Output Writer Module...")
            self.data_manager.save_data(f"{self.data_dir}\\output")
        except Exception as e:
            MAIN_LOGGER.error(f"Stopping execution of program due to error in execution of Output Writer Module: "
                              f"{str(e)}")
            print(traceback.format_exc())
            exit(-1)

        MAIN_LOGGER.info(f"Successfully executed GNSS algorithm {model}")

    def _check_data_dir(self, data_dir):
        """
        check if data_dir is a valid dir. If not, use the default dir.
        check if the data_dir exists. If not, create it.
        Args:
            data_dir: all generated files are saved in data_dir
        Returns:
            data_dir: valid data dir.
        """
        # check data dir
        # data_dir is not specified, automatically create one
        if data_dir is None or data_dir == '' or data_dir == "default":
            data_dir = str(RUNS_PATH)
            if data_dir[-1] != '//':
                data_dir = data_dir + '//'
            data_dir = data_dir + time.strftime('%Y-%m-%dT%HH%MM%SS', time.localtime()) + '//'
            data_dir = os.path.abspath(data_dir)

        # try to create data dir
        if not os.path.exists(data_dir):
            try:
                data_dir = os.path.abspath(data_dir)
                os.makedirs(data_dir)
                os.makedirs(f"{data_dir}\\trace")
                os.makedirs(f"{data_dir}\\output")
            except:
                raise IOError(f"Cannot create dir: {data_dir}")
        return data_dir

    @classmethod
    def compute_dop(cls, log, data_manager):
        log.info("Computing Dilution of Precision (DOP) metrics in ECEF and local (ENU) frames")
        sol = data_manager.get_data("nav_solution")

        for state in sol:
            dop_matrix = state.get_additional_info("dop_matrix")  # DOP matrix

            # DOPs
            geometry_dop = np.sqrt(dop_matrix[0, 0] + dop_matrix[1, 1] + dop_matrix[2, 2] + dop_matrix[3, 3])
            position_dop = np.sqrt(dop_matrix[0, 0] + dop_matrix[1, 1] + dop_matrix[2, 2])
            time_dop = np.sqrt(dop_matrix[3, 3])
            x_ecef = np.sqrt(dop_matrix[0, 0])
            y_ecef = np.sqrt(dop_matrix[1, 1])
            z_ecef = np.sqrt(dop_matrix[2, 2])
            dop_ecef = {"geometry": geometry_dop,
                        "position": position_dop,
                        "time": time_dop,
                        "x_ecef": x_ecef,
                        "y_ecef": y_ecef,
                        "z_ecef": z_ecef}
            state.add_additional_info("dop_ecef", dop_ecef)

            # computing DOPs in ENU frame
            lat, long, _ = cartesian2geodetic(*state.position)
            R = latlon2dcm_e_enu(lat, long)  # rotation matrix from ECEF to ENU
            dop_enu = R @ dop_matrix[0:3, 0:3] @ R.T

            dop_east = np.sqrt(dop_enu[0, 0])
            dop_north = np.sqrt(dop_enu[1, 1])
            dop_up = np.sqrt(dop_enu[2, 2])
            dop_horizontal = np.sqrt(dop_enu[0, 0] + dop_enu[1, 1])
            dop_enu = {"east": dop_east,
                       "north": dop_north,
                       "up": dop_up,
                       "horizontal": dop_horizontal}
            state.add_additional_info("dop_local", dop_enu)

            log.info(f"DOPs for epoch {str(state.epoch)}: geometry = {geometry_dop}, horizontal = {dop_horizontal}, "
                     f"vertical = {dop_up}")

    def compute(self, data_manager, trace_path):
        log = get_logger("MAIN_LOG")

        # get the input raw obs data
        nav_data = data_manager.get_data("nav_data")

        # perform pre-processing here
        log.info(f"Starting Preprocessor Module")
        preprocessor = PreprocessorManager(trace_path, data_manager)
        preprocessor.compute()  # this is the gnss_models data to actually process

        # run estimation algorithm
        log.info(f"Running estimation algorithm...")
        solver = GnssSolver(data_manager.get_clean_obs_data(), data_manager.get_raw_obs_data(), nav_data)
        solver.solve()

        data_manager.add_data("nav_solution", solver.solution)

        # compute DOPs in ECEF and local (ENU) frame
        self.compute_dop(log, data_manager)
