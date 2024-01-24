from datetime import timedelta

from . import Functor

from src.data_types.gnss.observation_data import ObservationData
from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.observation import Observation
from src.errors import NonExistentObservable

# The SmoothFunctor class implements the Hatch filter, which smooths pseudorange observables using time differences
# of carrier-phase observables, thus reducing raw pseudorange noise.


class SmoothFunctor(Functor):

    def __init__(self, time_constant: float, obs_rate: float):
        super().__init__()
        self.rate = obs_rate
        self.time_constant = time_constant

        # local map to save the SPRs
        self._SPRs = ObservationData()

    @staticmethod
    def _get_carrier(v_obs, freq):
        for obs in v_obs:
            if DataType.is_carrier(obs.datatype) or DataType.is_iono_free_carrier(obs.datatype):
                if obs.datatype.freq_number == freq:
                    return obs
        return None

    def get_alfa(self, t, t_prev, t_0):
        """
        alfa = (t_i - t_{i-1}) / (t_i - t_0) if t_i - t_0 < Tc
        alfa = (t_i - t_{i-1}) / Tc if t_i - t_0 >= Tc
        """
        delta_initial = (t - t_0).total_seconds()
        delta = (t - t_prev).total_seconds()

        if delta_initial < self.time_constant:
            return delta / delta_initial
        else:
            return delta / self.time_constant

    def __call__(self, obs_data_in, epoch, sat):
        # smooth function:
        # SPR(t_i) = alfa PR(t_i) + (1 - alfa) * (SPR(t_{i-1}) + L(t_i) - L(t_{i-1})
        v_obs_out = []
        epoch_prev = epoch + timedelta(seconds=-self.rate)

        # get observations at this epoch and previous epoch (Input Raw ObservationData)
        try:
            raw_obs = obs_data_in.get_observables_at_epoch(epoch, sat)
        except NonExistentObservable:
            # in case we don't have data for the current epoch, we cannot apply smoothing
            return []  # safe return

        # iterate over all observations and only process codes
        for obs in raw_obs:
            if DataType.is_code(obs.datatype) or DataType.is_iono_free_code(obs.datatype):
                pseudorange = obs  # PR(t_i)
                carrier = self._get_carrier(raw_obs, pseudorange.datatype.freq_number)  # L(t_i)

                # compute smooth observable
                smooth_obs = self.smooth_function(obs_data_in, sat, epoch, epoch_prev, pseudorange, carrier)

                # append it to output list
                v_obs_out.append(smooth_obs)

        return v_obs_out

    def smooth_function(self, obs_data_in, sat, epoch, epoch_prev, pseudorange, carrier):

        # get SPR datatype (C1 + L1 = SPR1, C2 + L2 = SPR2, C12 + L12 = SPR12, ...)
        smooth_type = DataType.get_smooth_datatype(pseudorange.datatype)

        try:
            # SPR(t_{i-1})
            spr_prev = self._SPRs.get_observable_at_epoch(sat, epoch_prev, smooth_type)
            carrier_prev = obs_data_in.get_observable_at_epoch(sat, epoch_prev, carrier.datatype)  # L(t_{i-1})

            # get alfa weight # TODO: review this
            alfa = self.get_alfa(epoch, epoch_prev, obs_data_in.get_first_arc_epoch(sat, epoch, self.rate))

            # smooth computation
            carrier_diff = carrier.value - carrier_prev.value
            value = alfa * pseudorange.value + (1 - alfa) * (spr_prev.value + carrier_diff)
            smooth_obs = Observation(smooth_type, value)

        except:
            # an exception was raised in one of the above functions.
            # possible reasons:
            #   * first epoch (no previous-epoch data)

            # set smooth data equal to raw pseudorange (unfiltered)
            smooth_obs = Observation(smooth_type, pseudorange.value)

        # save this SPR in the local functor memory
        self._SPRs.set_observation(epoch, sat, smooth_obs)

        return smooth_obs
