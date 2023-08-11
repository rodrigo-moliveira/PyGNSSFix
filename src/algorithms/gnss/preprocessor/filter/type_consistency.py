from . import Filter


class TypeConsistencyFilter(Filter):

    def __init__(self, types):
        super().__init__()
        self.types = types

    def is_applicable(self, sat, epoch, observation, **kwargs):
        # return False to keep this observable

        # check if this observation is part of the required observations
        ret_val = observation.datatype not in self.types[sat.sat_system]

        if not ret_val:
            # check if all required observations (self.types) are available for this epoch
            for obs_required in self.types[sat.sat_system]:
                ret_val = True
                for obs in kwargs["obs_list"]:
                    if obs.datatype == obs_required:
                        ret_val = False
                        break
        return ret_val
