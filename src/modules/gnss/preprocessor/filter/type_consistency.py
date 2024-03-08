from src.data_types.gnss.data_type import DataType
from . import Filter


class TypeConsistencyFilter(Filter):

    def __init__(self, types):
        super().__init__()
        self.types = types

        # build required datatypes (pseudorange and carrier phase)
        # NOTE: Doppler and SNR are not mandatory
        self.mandatory = {}
        for const, observations in self.types.items():
            self.mandatory[const] = []
            for obs in observations:
                if DataType.is_code(obs) or DataType.is_carrier(obs):
                    self.mandatory[const].append(obs)

    def is_applicable(self, sat, epoch, observation, **kwargs):
        # return False to keep this observable

        # check if this gnss_models is part of the required observations
        ret_val = observation.datatype not in self.types[sat.sat_system]

        if not ret_val:

            # check if all mandatory observations are available for this epoch
            for obs_required in self.mandatory[sat.sat_system]:
                ret_val = True
                for obs in kwargs["obs_list"]:
                    if obs.datatype == obs_required:
                        ret_val = False
                        break
        return ret_val
