from . import Functor
from PositioningSolver.src.data_types.basics.DataType import DataType
from PositioningSolver.src.gnss.data_types.Observation import Observation


class IonoFreeFunctor(Functor):

    def __init__(self, observations):
        super().__init__()

        self.observations = list(observations)
        try:
            self.base_freq_index = int(self.observations[0][0])
            self.second_freq_index = int(self.observations[1][0])
        except Exception as e:
            raise AttributeError(f"Only one frequency is available, {self.observations}. User needs to select two "
                                 f"frequencies for the computation of Iono Free observables. "
                                 f"\nTraceback Reason: {e}")

    def get_iono_free_datatype(self, v_obs_in):
        # get code and carrier for first frequency (C1, L1)
        C1 = C2 = L1 = L2 = None

        for obs in v_obs_in:
            if obs.datatype.freq_number == self.base_freq_index:
                if DataType.is_code(obs.datatype):
                    C1 = obs
                elif DataType.is_carrier(obs.datatype):
                    L1 = obs

        # get code and carrier for second frequency (C2, L2)
        for obs in v_obs_in:
            if obs.datatype.freq_number == self.second_freq_index:
                if DataType.is_code(obs.datatype):
                    C2 = obs
                elif DataType.is_carrier(obs.datatype):
                    L2 = obs

        return C1, C2, L1, L2

    @staticmethod
    def compute_iono_free(type1: Observation, type2: Observation):
        # get frequency values
        f1 = type1.datatype.freq.freq_value
        f2 = type2.datatype.freq.freq_value

        gama1 = f1 * f1 / (f1 * f1 - f2 * f2)
        gama2 = f2 * f2 / (f1 * f1 - f2 * f2)

        iono_free = gama1 * type1.value - gama2 * type2.value
        return iono_free

    def __call__(self, obs_data_in, epoch, sat):
        v_obs_in = obs_data_in.get_observables_at_epoch(epoch, sat)
        v_obs_out = []

        C1, C2, L1, L2 = self.get_iono_free_datatype(v_obs_in)

        # get iono free code
        if C1 is not None and C2 is not None:
            C12 = DataType.get_iono_free_datatype(C1.datatype, C2.datatype)

            # get iono-free value
            iono_free = self.compute_iono_free(C1, C2)
            v_obs_out.append(Observation(C12, iono_free))

        # get iono free carrier
        if L1 is not None and L2 is not None:
            L12 = DataType.get_iono_free_datatype(L1.datatype, L2.datatype)

            # get iono-free value
            iono_free = self.compute_iono_free(L1, L2)
            v_obs_out.append(Observation(L12, iono_free))

        return v_obs_out
