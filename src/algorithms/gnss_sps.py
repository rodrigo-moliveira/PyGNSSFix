# GNSS Single Point Solution Algorithm

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
