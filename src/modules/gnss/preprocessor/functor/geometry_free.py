""" Module with the implementation of the Geometry-Free Functor """

from src.io.config import config_dict
from src.common_log import MODEL_LOG, get_logger
from . import Functor
from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.observation import Observation


class GeometryFreeFunctor(Functor):
    """
    The `GeometryFreeFunctor` computes geometry-free observations.

    The pseudorange and carrier phase geometry-free combinations are computed with the following equations:

        GF_PR12 = PR2 - PR1
        GF_CP12 = CP1 - CP2

    where PRi and CPi are the pseudorange and carrier phase observables for frequency fi.
    Note the change in the order of the terms in the equations.

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
        self.pr_gf = None
        self.cp1 = None
        self.cp2 = None
        self.cp_gf = None

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
            self.pr_gf = DataType.get_gf_datatype(self.pr1, self.pr2, constellation)
            self.cp_gf = DataType.get_gf_datatype(self.cp1, self.cp2, constellation)
        else:
            raise AttributeError(f"Error computing geometry-free data for constellation {constellation}: "
                                 f"Mismatch between observation std list {std_dict} and base and second frequencies "
                                 f"{self.base_freq} and {self.second_freq}")

    def _get_observations(self, v_obs_in) -> list[Observation]:
        """ Gets the observations matching each datatype from the input list of Observations """
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

    def __call__(self, obs_data_in, epoch, sat):
        v_obs_in = obs_data_in.get_observables_at_epoch(epoch, sat)
        v_obs_out = []

        C1, C2, L1, L2 = self._get_observations(v_obs_in)

        # get iono free code
        if C1 is not None and C2 is not None:
            # get geometry-free value
            gf_obs = C2.value - C1.value
            v_obs_out.append(Observation(self.pr_gf, gf_obs))

        if L1 is not None and L2 is not None:
            # get geometry-free value
            gf_obs = L1.value - L2.value
            v_obs_out.append(Observation(self.cp_gf, gf_obs))

        return v_obs_out
