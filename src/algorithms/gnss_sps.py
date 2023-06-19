# GNSS Single Point Solution Algorithm

from .algorithm import Algorithm
from .gnss.preprocessor.preprocessor_manager import PreprocessorManager
from ..common_log import get_logger
from src.io.config import config_dict

class GnssSinglePointSolution(Algorithm):
    def __init__(self):
        super().__init__()
        self.inputs = {"observation_data": None, "navigation_data": None}
        self.outputs = {"nav_solution": None}
        self.name = "GNSS Single Point Solution Algorithm"

    def __str__(self):
        return f"Algorithm({self.name})"

    def compute(self, data_manager, trace_path):
        log = get_logger("GNSS_ALG")

        # get the input raw obs data
        raw_obs_data = data_manager.get_data("obs_data")

        # perform pre-processing here
        log.info(f"Starting Preprocessor Module")
        preprocessor = PreprocessorManager(trace_path, raw_obs_data)
        obs_data = preprocessor.compute()  # this is the observation data to actually process

        # run algorithm
        log.info(f"Running estimation algorithm...")
        # ...

        # Perform Least-Squares over input pseudorange observables..

        for time in range(1, 100):
            pass

        self.outputs["nav_solution"] = 1
