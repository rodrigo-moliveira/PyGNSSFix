# GNSS Single Point Solution Algorithm
import numpy as np

from .algorithm import Algorithm
from .gnss.solver.gnss_solver import GnssSolver
from .gnss.preprocessor.preprocessor_manager import PreprocessorManager
from ..common_log import get_logger
from ..models.frames import cartesian2geodetic, latlon2dcm_e_enu


class GnssSinglePointSolution(Algorithm):
    def __init__(self):
        super().__init__()
        self.name = "GNSS Single Point Solution Algorithm"

    def __str__(self):
        return f"{self.name}"

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
        raw_obs_data = data_manager.get_data("obs_data")
        nav_data = data_manager.get_data("nav_data")

        # perform pre-processing here
        log.info(f"Starting Preprocessor Module")
        preprocessor = PreprocessorManager(trace_path, data_manager)
        preprocessor.compute()  # this is the gnss_obs data to actually process
        exit()
        # run estimation algorithm
        log.info(f"Running estimation algorithm...")
        solver = GnssSolver(obs_data, nav_data)
        solver.solve()

        data_manager.add_data("nav_solution", solver.solution)

        # compute DOPs in ECEF and local (ENU) frame
        self.compute_dop(log, data_manager)
