""" GNSS Algorithm Manager Module """

import traceback
import numpy as np

from src.io.config import config_dict, EnumAlgorithmPNT
from src.models.frames import cartesian2geodetic, latlon2dcm_e_enu
from src.errors import ConfigError
from src.common_log import MAIN_LOG, get_logger

from .solver.gnss_solver import GnssSolver
from .preprocessor import PreprocessorManager
from ...data_mng.gnss.gnss_data_mng import GnssDataManager


class GnssAlgorithmManager:
    """
    GNSS Algorithm Manager.

    This class executes the GNSS PNT solver algorithm. It is structured in the following tasks:
        * Constructs the :py:class:`GnssDataManager` instance
        * Reads all the input data (:py:meth:`GnssDataManager.read_inputs`)
        * Launches the main PNT algorithm (function `compute`), which performs the following sub-tasks:
            * launches the :py:class:`PreprocessorManager` algorithm
            * launches the :py:class:`GnssSolver` algorithm
            * compute DOPs (see :py:meth:`compute_dop`)

        Attributes:
            data_dir(str): absolute path where the output data is stored
            data_manager(GnssDataManager): instance of the GNSS Data Manager
            main_log(logging.Logger): logger instance
    """

    def __init__(self, data_dir):
        # create data members
        self.data_manager = GnssDataManager()
        self.main_log = None
        self.data_dir = data_dir

    def run(self):
        """ Main function that executes the GNSS Algorithm """
        self.main_log = get_logger(MAIN_LOG)
        self.main_log.info("Starting GNSS Algorithm Manager")
        gnss_alg = config_dict.get('gnss_alg')

        if gnss_alg not in (EnumAlgorithmPNT.SPS, EnumAlgorithmPNT.PR_PPP, EnumAlgorithmPNT.CP_PPP):
            raise ConfigError(f"Selected Model {gnss_alg} not valid. Available options are "
                              f"SPS, PR-PPP, CP-PPP.")
        self.main_log.info(f"Running GNSS algorithm {gnss_alg}")

        config_dict.get_obs_std()

        # Input Reader Module
        try:
            self.main_log.info(f"Starting Input Reader Module...")
            self.data_manager.read_inputs(gnss_alg, f"{self.data_dir}\\trace")
        except Exception as e:
            self.main_log.error(f"Stopping execution of program due to error in execution of Input Reader Module: {e}")
            print(traceback.format_exc())
            exit(-1)

        # Main Algorithm Module
        try:
            self.main_log.info(f"Starting Main Algorithm Module...")
            self._compute(self.data_manager, f"{self.data_dir}\\trace")
        except Exception as e:
            self.main_log.error(f"Stopping execution of program due to error in execution of Main Algorithm "
                                f"Module: {e}")
            print(traceback.format_exc())
            exit(-1)

        # Output Writer Module
        try:
            self.main_log.info(f"Starting Output Writer Module...")
            self.data_manager.save_data(f"{self.data_dir}\\output")
        except Exception as e:
            self.main_log.error(f"Stopping execution of program due to error in execution of Output Writer Module: "
                                f"{str(e)}")
            print(traceback.format_exc())
            exit(-1)

        self.main_log.info(f"Successfully executed GNSS algorithm {gnss_alg}")

    def _compute_dop(self):
        """ Internal function for the computation of the DOPs """
        self.main_log.info("Computing Dilution of Precision (DOP) metrics in ECEF and local (ENU) frames")
        sol = self.data_manager.get_data("nav_solution")

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

            self.main_log.info(f"DOPs for epoch {str(state.epoch)}: geometry = {geometry_dop}, "
                               f"horizontal = {dop_horizontal}, vertical = {dop_up}")

    def _compute(self, data_manager, trace_path):
        """ Internal function for the computation of the PNT solver task """
        # perform pre-processing here
        self.main_log.info(f"Starting Preprocessor Module")
        preprocessor = PreprocessorManager(trace_path, data_manager)
        preprocessor.compute()

        # run estimation algorithm
        self.main_log.info(f"Running estimation algorithm...")
        solver = GnssSolver(data_manager, trace_path)
        solver.solve()

        data_manager.add_data("nav_solution", solver.solution)

        # compute DOPs in ECEF and local (ENU) frame
        self._compute_dop()
