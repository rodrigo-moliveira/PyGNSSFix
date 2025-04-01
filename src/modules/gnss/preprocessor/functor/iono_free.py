""" Module with the implementation of the Ionospheric Free Functor """

from src.io.config import config_dict
from src.common_log import MODEL_LOG, get_logger
from . import Functor
from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.observation import Observation


class IonoFreeFunctor(Functor):
    """
    The `IonoFreeFunctor` computes 1st order ionospheric free observations.

    The pseudorange and carrier phase iono free combinations are computed with the following equations:

        PR12 = (f1^2 * PR1 - f2^2 * PR2) / (f1^2 - f2^2)
        CP12 = (f1^2 * CP1 - f2^2 * CP2) / (f1^2 - f2^2)

    where PRi and CPi are the pseudorange and carrier phase observables for frequency fi.

    More information can be found in:
    https://gssc.esa.int/navipedia/index.php/Ionosphere-free_Combination_for_Dual_Frequency_Receivers
    """
    def __init__(self, constellation: str, base_freq: DataType, second_freq: DataType):
        super().__init__()
        self.constellation = constellation
        self.base_freq = base_freq
        self.second_freq = second_freq
        self.gama1 = 0
        self.gama2 = 0
        self.pr1 = None
        self.pr2 = None
        self.pr_if = None
        self.cp1 = None
        self.cp2 = None
        self.cp_if = None
        self._compute_iono_std(constellation)

    def _compute_iono_std(self, constellation):
        """ Definition of auxiliary data """
        std_dict = config_dict.get_obs_std()[constellation]
        datatypes = list(std_dict.keys())

        for datatype in datatypes:
            if datatype.freq == self.base_freq and DataType.is_code(datatype):
                self.pr1 = datatype
                self.cp1 = DataType.get_carrier_from_code(self.pr1)
            if datatype.freq == self.second_freq and DataType.is_code(datatype):
                self.pr2 = datatype
                self.cp2 = DataType.get_carrier_from_code(self.pr2)

        if self.pr1 is not None and self.pr2 is not None:
            pr_std_1 = std_dict[self.pr1]
            f1 = self.pr1.freq.freq_value

            pr_std_2 = std_dict[self.pr2]
            f2 = self.pr2.freq.freq_value

            self.gama1 = f1 * f1 / (f1 * f1 - f2 * f2)
            self.gama2 = f2 * f2 / (f1 * f1 - f2 * f2)

            pr_std_if = abs(self.gama1 * pr_std_1 - self.gama2 * pr_std_2)
            self.pr_if = DataType.get_iono_free_datatype(self.pr1, self.pr2, constellation)
            self.cp_if = DataType.get_iono_free_datatype(self.cp1, self.cp2, constellation)
            config_dict.update_obs_std(constellation, self.pr_if, pr_std_if)

            log = get_logger(MODEL_LOG)
            log.info(f"Computing iono-free observation std for constellation {constellation}. {self.pr1}: {pr_std_1} - "
                     f"{self.pr2}: {pr_std_2} -> {self.pr_if}: {pr_std_if}")
        else:
            raise AttributeError(f"Error computing observation stds for constellation {constellation}: "
                                 f"Mismatch between observation std list {std_dict} and base and second frequencies "
                                 f"{self.base_freq} and {self.second_freq}")

    def _get_iono_free_observations(self, v_obs_in) -> list[DataType]:
        """ Gets the Iono Free DataType """
        # get code and carrier for first frequency (C1, L1)
        C1 = C2 = L1 = L2 = None

        for obs in v_obs_in:
            if obs.datatype == self.pr1:
                C1 = obs
            if obs.datatype == self.pr2:
                C2 = obs
            if obs.datatype == self.cp1:
                L1 = obs
            if obs.datatype == self.cp2:
                L2 = obs

        return C1, C2, L1, L2

    def _compute_iono_free(self, obs1: Observation, obs2: Observation) -> float:
        """ Computes the iono free combination (value) """
        iono_free = self.gama1 * obs1.value - self.gama2 * obs2.value
        return iono_free

    def __call__(self, obs_data_in, epoch, sat):
        v_obs_in = obs_data_in.get_observables_at_epoch(epoch, sat)
        v_obs_out = []

        C1, C2, L1, L2 = self._get_iono_free_observations(v_obs_in)

        # get iono free code
        if C1 is not None and C2 is not None:
            # get iono-free value
            iono_free = self._compute_iono_free(C1, C2)
            v_obs_out.append(Observation(self.pr_if, iono_free))

        if L1 is not None and L2 is not None:
            # get iono-free value
            iono_free = self._compute_iono_free(L1, L2)
            v_obs_out.append(Observation(self.cp_if, iono_free))

        return v_obs_out
