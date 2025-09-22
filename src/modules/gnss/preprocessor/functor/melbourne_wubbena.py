""" Module with the implementation of the Melbourne–Wubbena Functor """

from src.data_mng.gnss import ObservationData
from . import Functor
from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.observation import Observation


class MWFunctor(Functor):
    """
    The `MWFunctor` computes Melbourne–Wubbena observations.

    The Melbourne–Wubbena combination is formed from the difference of narrow-lane pseudorange observations
    and wide-lane carrier-phase observations:

        MW = CP_WL - PR_NL

    Reference:
        [1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
            Springer Cham, 2017

    More information can be found in:
    https://gssc.esa.int/navipedia/index.php/Combination_of_GNSS_Measurements
    """
    def __init__(self, constellation: str, nl_obs_data: ObservationData, wl_obs_data: ObservationData):
        super().__init__()
        self.constellation = constellation
        self.nl_obs_data = nl_obs_data
        self.wl_obs_data = wl_obs_data
        self.pr_nl = None
        self.cp_wl = None
        self.mw = None

        nl_types = nl_obs_data.get_types(constellation)
        wl_types = wl_obs_data.get_types(constellation)

        for nl in nl_types:
            if DataType.is_nl_code(nl):
                self.pr_nl = nl
                break
        for wl in wl_types:
            if DataType.is_wl_carrier(wl):
                self.cp_wl = wl
                break

        # Perform validations
        if self.pr_nl is None or self.cp_wl is None:
            raise AttributeError(f"Error computing Melbourne-Wubbena observation: "
                                 f"Mismatch between observation types {nl_types} and {wl_types} and expected "
                                 f"narrow-lane code {self.pr_nl} and wide-lane carrier {self.cp_wl}")
        if self.pr_nl.freq_number != self.cp_wl.freq_number:
            raise AttributeError(f"Error computing Melbourne-Wubbena observation: "
                                 f"Mismatch between frequencies {self.pr_nl.freq_number} and {self.cp_wl.freq_number} "
                                 f"of narrow-lane code {self.pr_nl} and wide-lane carrier {self.cp_wl}")
        # Fetch MW DataType
        self.mw = DataType.get_mw_datatype(self.pr_nl, self.cp_wl, constellation)

    @staticmethod
    def _compute_mw(pr_nl: Observation, cp_wl: Observation) -> float:
        """ Computes the melbourne-wubbena combination (value) """
        mw_obs = cp_wl.value - pr_nl.value
        return mw_obs

    def __call__(self, obs_data_in, epoch, sat):
        v_obs_out = []
        if sat.sat_system == self.constellation:
            pr_nl = self.nl_obs_data.get_observable_at_epoch(sat, epoch, self.pr_nl)
            cp_wl = self.wl_obs_data.get_observable_at_epoch(sat, epoch, self.cp_wl)

            mw_obs = self._compute_mw(pr_nl, cp_wl)
            v_obs_out.append(Observation(self.mw, mw_obs))

        return v_obs_out
