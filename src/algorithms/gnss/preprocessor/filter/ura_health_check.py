from src.data_types.gnss.data_type import DataType
from . import Filter


class SatFilterHealthURA(Filter):

    def __init__(self, navigation_data, ura_check, ura_threshold, health_check):
        super().__init__()
        self.navigation_data = navigation_data
        self.ura_check = ura_check
        self.ura_threshold = ura_threshold
        self.health_check = health_check

    def apply(self, sat, epoch, observation, v_removable):
        # return False to keep this observable

        # only apply this filter to Signal Observables
        if not DataType.is_signal(observation.datatype):
            return

        if observation.value < self.snr_threshold:
            # remove all observables associated with this frequency
            v_observables = self.observation_data.get_observables_at_epoch(epoch, sat)

            # iterate over all observables and flag those with the same frequency as 'observation'
            for this_obs in v_observables:
                if this_obs.datatype.freq == observation.datatype.freq:
                    if this_obs not in v_removable:
                        v_removable.append(this_obs)
