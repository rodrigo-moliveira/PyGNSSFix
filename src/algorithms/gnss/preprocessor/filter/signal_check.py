from PositioningSolver.src.data_types.basics.DataType import DataType
from . import Filter
from PositioningSolver.src.config import config


class SignalCheckFilter(Filter):

    def __init__(self, observation_data):
        super().__init__()
        self.observation_data = observation_data
        self.snr_threshold = config["gps_solver"]["signal_strength_filter"]["select"]

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
