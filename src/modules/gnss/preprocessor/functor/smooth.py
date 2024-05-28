"""Smoother Module with the implementation of the Hatch Filter"""

from datetime import timedelta

from . import Functor

from src.data_mng.gnss.observation_data import ObservationData
from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.observation import Observation
from src.errors import TimeSeriesError


class SmoothFunctor(Functor):
    """
    The SmoothFunctor class implements the Hatch filter, which smooths pseudorange observables using time differences
    of carrier-phase observables, thus reducing raw pseudorange noise.

    The smooth pseudorange observation is computed for epoch t_i as:
        SPR(t_i) = alfa PR(t_i) + (1 - alfa) * (SPR(t_{i-1}) + L(t_i) - L(t_{i-1})

    where `alfa` is the driver of the smoother and is defined as:
        alfa = (t_i - t_{i-1}) / (t_i - t_0) if t_i - t_0 < Tc
        alfa = (t_i - t_{i-1}) / Tc if t_i - t_0 >= Tc

    Tc is the user-defined `time_constant`, and t_0 is the initial epoch

    NOTE:
        Currently, no cycle slip detection is implemented, and thus problems will arise if cycle slips appear

    References:
        [1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
            Springer Cham, 2017 - Chapter 20.4

    """

    def __init__(self, time_constant: float, obs_rate: float):
        """Constructor of the SmoothFunctor algorithm

        Parameters:
            time_constant(float): user-defined time constant for the Hatch Filter, in seconds
            obs_rate(float): data rate of the observation data (period of observations), in seconds
        """
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
        delta_initial = (t - t_0).total_seconds()
        delta = (t - t_prev).total_seconds()

        if delta_initial < self.time_constant:
            return delta / delta_initial
        else:
            return delta / self.time_constant

    def __call__(self, obs_data_in, epoch, sat):
        v_obs_out = []
        epoch_prev = epoch + timedelta(seconds=-self.rate)

        # get observations at this epoch
        try:
            raw_obs = obs_data_in.get_observables_at_epoch(epoch, sat)
        except TimeSeriesError:
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
            spr_prev = self._SPRs.get_observable_at_epoch(sat, epoch_prev, smooth_type)  # SPR(t_{i-1})
            carrier_prev = obs_data_in.get_observable_at_epoch(sat, epoch_prev, carrier.datatype)  # L(t_{i-1})

            # get alfa weight
            alfa = self.get_alfa(epoch, epoch_prev, obs_data_in.get_first_arc_epoch(sat, epoch, self.rate))

            # smooth computation
            carrier_diff = carrier.value - carrier_prev.value
            value = alfa * pseudorange.value + (1 - alfa) * (spr_prev.value + carrier_diff)
            smooth_obs = Observation(smooth_type, value)

        except (TimeSeriesError, KeyError):
            # an exception was raised in one of the above functions (no data for the previous epoch).
            # ACTION: set smooth data equal to raw pseudorange (unfiltered)
            smooth_obs = Observation(smooth_type, pseudorange.value)

        # save this SPR in the local functor memory
        self._SPRs.set_observation(epoch, sat, smooth_obs)

        return smooth_obs
