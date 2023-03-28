# GNSS Single Point Solution Algorithm

from .algorithm import Algorithm


class GnssSinglePointSolution(Algorithm):
    def __init__(self):
        super().__init__()
        self.inputs = {"observation_data": None, "navigation_data": None}
        self.outputs = {"output_pos": None}
        self.name = "GNSS Single Point Solution Algorithm"

    def __str__(self):
        return f"Algorithm({self.name})"

    def compute(self):
        # compute results and append them

        # ...

        # Perform Least-Squares over input pseudorange observables..

        self.get_results()["output_pos"] = 1
