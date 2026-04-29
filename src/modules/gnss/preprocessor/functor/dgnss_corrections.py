""" Module with the implementation of the DGNSS Corrections """
import numpy as np

from src.data_mng.gnss import GnssStateSpace
from src.data_mng.gnss.geometry import SatelliteGeometry
from src.io.config import config_dict, EnumTransmissionTime
from . import Functor
from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.observation import Observation


class DGNSSFunctor(Functor):
    """
    The `DGNSSFunctor` computes ...

    ...

    More information can be found in:
    TBC
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
