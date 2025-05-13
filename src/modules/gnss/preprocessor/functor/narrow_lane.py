""" Module with the implementation of the Narrow-Lane Functor """
from math import sqrt

from src.io.config import config_dict
from src.common_log import MODEL_LOG, get_logger
from . import Functor
from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.observation import Observation


class NLFunctor(Functor):
    """
    The `NLFunctor` computes narrow-lane observations.

    The pseudorange and carrier phase narrow-lane combinations are computed with the following equations:

        PR_NL = f1 * PR1 / (f1 + f2) + f2 * PR2 / (f1 + f2)
        CP_NL = f1 * CP1 / (f1 + f2) + f2 * CP2 / (f1 + f2)

    where PRi and CPi are the pseudorange and carrier phase observables for frequency fi.

    Reference:
        [1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
            Springer Cham, 2017

    More information can be found in:
    https://gssc.esa.int/navipedia/index.php/Combination_of_GNSS_Measurements
    """
    def __init__(self, constellation: str, base_freq: DataType, second_freq: DataType):
        super().__init__()
        self.constellation = constellation
        self.base_freq = base_freq
        self.second_freq = second_freq
        self.pr1 = None
        self.pr2 = None
        self.pr_nl = None
        self.cp1 = None
        self.cp2 = None
        self.cp_nl = None
        self._compute_nl_std(constellation)

    def _compute_nl_std(self, constellation):
        """ Definition of auxiliary data """
        std_dict = config_dict.get_obs_std()[constellation]
        datatypes = list(std_dict.keys())

        for datatype in datatypes:
            if DataType.is_code(datatype) and datatype.freq == self.base_freq:
                self.pr1 = datatype
                self.cp1 = DataType.get_carrier_from_code(self.pr1)
            if DataType.is_code(datatype) and datatype.freq == self.second_freq:
                self.pr2 = datatype
                self.cp2 = DataType.get_carrier_from_code(self.pr2)

        if self.pr1 is not None and self.pr2 is not None:
            pr_std_1 = std_dict[self.pr1]
            cp_std_1 = std_dict[self.cp1]
            f1 = self.pr1.freq.freq_value

            pr_std_2 = std_dict[self.pr2]
            cp_std_2 = std_dict[self.cp2]
            f2 = self.pr2.freq.freq_value

            # Eq. (20.28) from [1]
            pr_std_nl = sqrt((f1**2 * pr_std_1**2 + f2**2 * pr_std_2**2) / (f1 + f2)**2)
            cp_std_nl = sqrt((f1**2 * cp_std_1**2 + f2**2 * cp_std_2**2) / (f1 + f2)**2)
            self.pr_nl = DataType.get_nl_datatype(self.pr1, self.pr2, constellation)
            self.cp_nl = DataType.get_nl_datatype(self.cp1, self.cp2, constellation)
            config_dict.update_obs_std(constellation, self.pr_nl, pr_std_nl)
            config_dict.update_obs_std(constellation, self.cp_nl, cp_std_nl)

            log = get_logger(MODEL_LOG)
            log.info(f"Computing narrow-lane observation std for constellation {constellation}. {self.pr1}: {pr_std_1}"
                     f" - {self.pr2}: {pr_std_2} -> {self.pr_nl}: {pr_std_nl} and {self.cp1}: {cp_std_1} - "
                        f"{self.cp2}: {cp_std_2} -> {self.cp_nl}: {cp_std_nl}")
        else:
            raise AttributeError(f"Error computing observation stds for constellation {constellation}: "
                                 f"Mismatch between observation std list {std_dict} and base and second frequencies "
                                 f"{self.base_freq} and {self.second_freq}")

    def _get_observations(self, v_obs_in) -> list[Observation]:
        """ Gets the observations matching each datatype from the input list of Observations """
        # get code and carrier for first and second frequencies
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

    @staticmethod
    def _compute_nl(obs1: Observation, obs2: Observation) -> float:
        """ Computes the narrow-lane combination (value) """
        f1 = obs1.datatype.freq.freq_value
        f2 = obs2.datatype.freq.freq_value
        # Eq. (20.27) of [1]
        nl_obs = (f1 * obs1.value + f2 * obs2.value) / (f1 + f2)
        return nl_obs

    def __call__(self, obs_data_in, epoch, sat):
        v_obs_in = obs_data_in.get_observables_at_epoch(epoch, sat)
        v_obs_out = []

        C1, C2, L1, L2 = self._get_observations(v_obs_in)

        # get narrow-lane code
        if C1 is not None and C2 is not None:
            nl_code = self._compute_nl(C1, C2)
            v_obs_out.append(Observation(self.pr_nl, nl_code))

        # get narrow-lane carrier
        if L1 is not None and L2 is not None:
            nl_carrier = self._compute_nl(L1, L2)
            v_obs_out.append(Observation(self.cp_nl, nl_carrier))

        return v_obs_out
