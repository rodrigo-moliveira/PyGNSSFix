# base filter class
from src.data_types.gnss.observation import Observation
from src.data_types.gnss.satellite import Satellite
from src.data_types.date.date import Epoch


# Filter classes inherited from this must implement the ``is_applicable`` method
class Filter:
    def __init__(self):
        pass

    def is_applicable(self, sat: Satellite, epoch: Epoch, observation: Observation, **kwargs):
        # return True to remove observables, return False to keep observable
        return True

    def apply(self, sat: Satellite, epoch: Epoch, observation: Observation, v_removable: list):
        v_removable.append(observation)
