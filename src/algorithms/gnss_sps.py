# GNSS Single Point Solution Algorithm
import numpy as np

from .algorithm import Algorithm
from .gnss.solver.gnss_solver import GnssSolver
from .gnss.preprocessor.preprocessor_manager import PreprocessorManager
from ..common_log import get_logger


class GnssSinglePointSolution(Algorithm):
    def __init__(self):
        super().__init__()
        self.name = "GNSS Single Point Solution Algorithm"

    def __str__(self):
        return f"{self.name}"

    def compute_dop(self, log, data_manager):
        log.info("Computing Dilution of Precision (DOP) metrics in ECEF and NED frames")
        sol = data_manager.get_data("nav_solution")

        for state in sol:
            system_matrix = state.get_solver_info("geometry_matrix")  # this is G matrix in LS
            dop_matrix = np.linalg.inv(system_matrix.T @ system_matrix)

            # geometry, position and time DOPs
            geometry_dop = np.sqrt(dop_matrix[0, 0] + dop_matrix[1, 1] + dop_matrix[2, 2] + dop_matrix[3, 3])
            position_dop = np.sqrt(dop_matrix[0, 0] + dop_matrix[1, 1] + dop_matrix[2, 2])
            time_dop = np.sqrt(dop_matrix[3, 3])

            # x, y, z DOPs
            x_ecef = np.sqrt(dop_matrix[0, 0])
            y_ecef = np.sqrt(dop_matrix[1, 1])
            z_ecef = np.sqrt(dop_matrix[2, 2])
            info = {"geometry": geometry_dop,
                    "position": position_dop,
                    "time": time_dop,
                    "x_ecef": x_ecef,
                    "y_ecef": y_ecef,
                    "z_ecef": z_ecef}
            state.add_solver_info("dop_ecef", info)

        """exit()
        for epoch, DOPs in self.items():
            # get receiver position
            receiver = receiver_pos.get_data_for_epoch(epoch).copy()
            receiver.form = "geodetic"
            lat, long, h = receiver

            # get ECEF DOP matrix
            DOP_matrix = DOPs.matrix



            # get DOPs with respect to ENU coordinates
            R = rot1((Constant.PI / 2 - lat)) @ rot3((Constant.PI / 2 + long))  # rotation matrix from ECEF to ENU
            DOP_ENU = R @ DOP_matrix[0:3, 0:3] @ R.T

            # east, north, up DOPs
            DOPs.east = np.sqrt(DOP_ENU[0, 0])
            DOPs.north = np.sqrt(DOP_ENU[1, 1])
            DOPs.up = np.sqrt(DOP_ENU[2, 2])

            # horizontal DOP
            DOPs.horizontal = np.sqrt(DOP_ENU[0, 0] + DOP_ENU[1, 1])
        """

    def compute(self, data_manager, trace_path):
        log = get_logger("MAIN_LOG")

        # get the input raw obs data
        raw_obs_data = data_manager.get_data("obs_data")
        nav_data = data_manager.get_data("nav_data")

        # perform pre-processing here
        log.info(f"Starting Preprocessor Module")
        preprocessor = PreprocessorManager(trace_path, raw_obs_data, nav_data)
        obs_data = preprocessor.compute()  # this is the observation data to actually process

        # run estimation algorithm
        log.info(f"Running estimation algorithm...")
        solver = GnssSolver(obs_data, nav_data)
        solver.solve()

        data_manager.add_data("nav_solution", solver.solution)

        # compute DOPs in ECEF and local (NED) frame
        self.compute_dop(log, data_manager)
