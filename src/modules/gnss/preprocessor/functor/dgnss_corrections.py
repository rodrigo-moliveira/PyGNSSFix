""" Module with the implementation of the D-GNSS Corrections """
import numpy as np

from src.data_mng.gnss import GnssStateSpace
from src.data_mng.gnss.geometry import SatelliteGeometry
from src.io.config import EnumTransmissionTime
from . import Functor
from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.observation import Observation


class DGNSSFunctor(Functor):
    """
    The `DGNSSFunctor` computes the differential GNSS observables, that is, the differential observations between the
    rover and the reference station. These are required by code based (D-GNSS) and phase based (RTK) processing
    algorithms.

    The differential pseudorange and carrier phase observations are computed with the following equations:

        PRC_PRi = rho - PRi
        PRC_CPi = rho - CPi

    where PRC_PRi and PRC_CPi are the pseudorange and carrier phase differential corrections for frequency fi.
    Moreover, rho is the true geometric range, and PRi and CPi are the raw PR and CP observations.
    Note that in the current software version, only code-based corrections are computed.

    The corrected pseudorange measurement applied to the rover is then computed as:
        * PR_corrected = PR_raw + PRC_PR, which is free of satellite-related clocks and biases

    More information can be found in:
        * [1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
            Springer Cham, 2017, Section 21.4.1
    """
    def __init__(self, true_position, sat_clocks, sat_orbits):
        super().__init__()
        self._true_position = true_position
        self._sat_clocks = sat_clocks
        self._sat_orbits = sat_orbits


    def __call__(self, obs_data_in, epoch, sat):
        v_obs_in = obs_data_in.get_observables_at_epoch(epoch, sat)
        v_obs_out = []

        geometry = SatelliteGeometry()
        state = GnssStateSpace(epoch=epoch)
        state.position = np.array(self._true_position)
        state.clock_bias = 0.0
        geometry.compute(sat, epoch, state, sat.sat_system, EnumTransmissionTime.GEOMETRIC, None,
                         self._sat_orbits, self._sat_clocks)
        true_range = geometry.true_range

        for obs in v_obs_in:
            prc_val = true_range - obs.value
            prc_datatype = DataType.get_dgnss_datatype(obs.datatype)
            v_obs_out.append(Observation(prc_datatype, prc_val))
        return v_obs_out
