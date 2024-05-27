""" Constellation Filter Module """
from src.modules.gnss.preprocessor.filter import Filter


# This is currently unused
class ConstellationFilter(Filter):

    def __init__(self, constellations):
        super().__init__()
        self.constellations = constellations

    def apply(self, sat, epoch, observation, v_removable):

        if sat.sat_system not in self.constellations:
            v_removable.append(observation)

