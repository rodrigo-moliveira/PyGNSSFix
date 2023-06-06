# GNSS Single Point Solution Algorithm

from .algorithm import Algorithm


class GnssSinglePointSolution(Algorithm):
    def __init__(self):
        super().__init__()
        self.inputs = {"observation_data": None, "navigation_data": None}
        self.outputs = {"nav_solution": None}
        self.name = "GNSS Single Point Solution Algorithm"

    def __str__(self):
        return f"Algorithm({self.name})"

    def compute(self, data_manager):
        # compute results and append them
        # get nav and obs data
        # perform pre-processing here..
        # run preprocessor
        main_log.info(f"Running preprocessor...")

        # run algorithm
        main_log.info(f"Running estimation algorithm...")
        # ...

        # Perform Least-Squares over input pseudorange observables..

        for time in range(1, 100):
            pass

        self.outputs["nav_solution"] = 1
