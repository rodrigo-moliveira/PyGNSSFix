from src.io.config import config_dict
from . import Functor
from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.observation import Observation


class IonoFreeFunctor(Functor):

    def __init__(self, constellation, observations):
        super().__init__()

        self.constellation = constellation
        self.observations = list(observations)
        try:
            self.base_freq_index = int(self.observations[0][0])
            self.second_freq_index = int(self.observations[1][0])
        except Exception as e:
            raise AttributeError(f"Problem getting base and second frequencies in Iono Free Computation for "
                                 f"constellation {self.constellation}, observations provided are {self.observations}. "
                                 f"Traceback Reason: {e}")

        self.compute_iono_std(constellation)

    def compute_iono_std(self, constellation):
        std_dict = config_dict.get_obs_std()[constellation]
        datatypes = list(std_dict.keys())
        if len(datatypes) != 2:
            raise AttributeError(f"Error computing observation stds for constellation {constellation}: "
                                 f"There should be 2 datatypes defined in std list {config_dict.get_obs_std()}")
        datatype_1 = datatypes[0]
        pr_std_1 = std_dict[datatype_1]
        f1 = datatype_1.freq.freq_value

        datatype_2 = datatypes[1]
        pr_std_2 = std_dict[datatype_2]
        f2 = datatype_2.freq.freq_value

        gama1 = f1 * f1 / (f1 * f1 - f2 * f2)
        gama2 = f2 * f2 / (f1 * f1 - f2 * f2)

        pr_std_if = gama1 * pr_std_1 - gama2 * pr_std_2
        datatype_if = DataType.get_iono_free_datatype(datatype_1, datatype_2, constellation)
        config_dict.set_obs_std(constellation, datatype_if, pr_std_if)

        # TODO: add log message "computing observation std for observations...."
        print(f"Computing iono-free observation std for constellation {constellation}. {datatype_1}: {pr_std_1} - "
              f"{datatype_2}: {pr_std_2} -> {datatype_if}: {pr_std_if}")

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
            C12 = DataType.get_iono_free_datatype(C1.datatype, C2.datatype, C1.datatype.constellation)

            # get iono-free value
            iono_free = self.compute_iono_free(C1, C2)
            v_obs_out.append(Observation(C12, iono_free))

        if L1 is not None and L2 is not None:
            L12 = DataType.get_iono_free_datatype(L1.datatype, L2.datatype, L1.datatype.constellation)

            # get iono-free value
            iono_free = self.compute_iono_free(L1, L2)
            v_obs_out.append(Observation(L12, iono_free))

        return v_obs_out
