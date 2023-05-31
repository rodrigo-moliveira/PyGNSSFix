# base filter class
from PositioningSolver.src.gnss.data_types.Observation import Observation
from PositioningSolver.src.gnss.data_types.Satellite import Satellite
from PositioningSolver.src.data_types.basics.Epoch import Epoch


# Filter classes inherited from this must implement the ``is_applicable`` method
class Filter:
    def __init__(self):
        pass

    def is_applicable(self, sat: Satellite, epoch: Epoch, observation: Observation):
        # return True to remove observables, return False to keep observable
        return True

    def apply(self, sat: Satellite, epoch: Epoch, observation: Observation, v_removable: list):
        v_removable.append(observation)
