from datetime import timedelta

import numpy as np

from .clock_obs import compute_tx_time
from .ephemeride_propagator import EphemeridePropagator
from ..frames.frames import enu2azel, ecef2enu, cartesian2geodetic
from ...data_mng.container import Container
from ...errors import TimeSeriesError


class SatelliteGeometry(Container):
    __slots__ = ["transit_time",
                 "time_emission", "time_reception",
                 "true_range", "az", "el", "satellite_position",
                 "dt_rel_correction", "los"]

    def __init__(self):
        super().__init__()
        self.transit_time = 0
        self.time_emission = None
        self.time_reception = None
        self.true_range = 0
        self.az = 0
        self.el = 0
        self.satellite_position = None
        self.dt_rel_correction = 0
        self.los = None

    def __str__(self):
        _allAttrs = ""
        for atr in self.__slots__:
            _allAttrs += atr + "=" + str(getattr(self, atr)) + ", "
        return f'{type(self).__name__}({_allAttrs[0:-2]})'

    def __repr__(self):
        return str(self)

    def compute(self, rec_pos, epoch, rec_bias, nav_message, compute_tx, PR_obs):
        """
        compute satellite-related quantities (tropo, iono, transmission time, etc.) to be used in the PVT gnss_obs
        reconstruction equation, for a given satellite.

        Args:
            rec_pos (numpy.ndarray) : Receiver position
            epoch (src.data_types.basics.Epoch.Epoch) : epoch under evaluation
            rec_bias (float) : Receiver clock bias
            nav_message (src.data_types.containers.NavigationData.NavigationPointGPS) : navigation point for the
                                                                                        satellite under evaluation
            compute_tx (function) : function to compute the transmission time
            PR_obs (src.data_types.data_types.Observation.Observation) : Code gnss_obs to use in some computations
        """
        # get reception time in GNSS time system ( T_GNSS = T_receiver - t_(receiver_bias) )
        # TODO: o ISB vai ter que entrar aqui para a slave constellation
        # para a slave constellation: rec_bias += ISB
        time_reception = epoch + timedelta(seconds=-rec_bias)

        # algorithm to compute Transmission time (in GNSS time)
        time_emission, transit = compute_tx_time(model=compute_tx,
                                                 r_receiver=rec_pos,
                                                 t_reception=epoch,
                                                 dt_receiver=rec_bias,
                                                 nav_message=nav_message,
                                                 pseudorange_obs=PR_obs)

        # get and satellite position at RX ECEF frame
        pos_sat, dt_relative = EphemeridePropagator.compute_sat_nav_position_dt_rel(
            nav_message, time_emission.gnss_time, transit)

        # compute true range
        true_range = np.linalg.norm(pos_sat - rec_pos)

        # get satellite elevation and azimuth angles (from the receiver), ENU frame
        lat, long, h = cartesian2geodetic(rec_pos[0], rec_pos[1], rec_pos[2])
        enu_coord = ecef2enu(pos_sat[0], pos_sat[1], pos_sat[2], lat, long, h)
        az, el = enu2azel(*enu_coord)

        # line of sight
        los = [(rec_pos[i] - pos_sat[i]) / true_range for i in (0, 1, 2)]

        # save results in container
        self.transit_time = transit
        self.time_emission = time_emission
        self.time_reception = time_reception
        self.true_range = true_range
        self.az = az
        self.el = el
        self.satellite_position = pos_sat
        self.dt_rel_correction = dt_relative
        self.los = los


class SystemGeometry:
    def __init__(self, nav_data, obs_data):
        """
        Args:
            nav_data (src.data_types.containers.NavigationData.NavigationDataMap) : Navigation data map
            obs_data (src.data_types.containers.ObservationData.EpochData) : gnss_obs epoch data for
        """
        self._data = dict.fromkeys(obs_data.get_satellites())
        self.nav_data = nav_data
        self.obs_data = obs_data

    def _clean(self):
        # reinitialize self._data
        vSats = self.get_satellites()
        self._data.clear()
        self._data = dict.fromkeys(vSats)

    def items(self):
        return self._data.items()

    def get_satellites(self):
        return list(self._data.keys())

    def remove(self, sat):
        if sat in self._data:
            self._data.pop(sat)

    def get(self, attribute, sat):
        if sat in self._data:
            return getattr(self._data[sat], attribute)
        return None

    def compute(self, epoch, state, metadata):
        """
        compute satellite-related quantities (tropo, iono, transmission time, etc.) to be used in the PVT gnss_obs
        reconstruction equation, for all available satellites.

        Args:
            epoch (src.data_types.basics.Epoch.Epoch) : epoch under evaluation
            state (src.data_types.state_space.statevector.Position) : state vector
            metadata (dict) : function to compute the transmission time
        """
        self._clean()
        _to_remove = []
        sat_list = self.get_satellites()

        for sat in sat_list:
            geometry = SatelliteGeometry()

            try:
                # fetch navigation message for this satellite
                nav_message = self.nav_data.get_closest_message(sat, epoch)
            except TimeSeriesError:
                _to_remove.append(sat)
                continue

            # fetch pseudorange gnss_obs for this satellite at epoch (used in the compute_TX_time algorithm)
            observable_lst = self.obs_data.get_code_observables(sat)
            if len(observable_lst) == 0:
                _to_remove.append(sat)
                continue

            # compute geometry for this satellite
            geometry.compute(state.position, epoch, state.clock_bias, nav_message,
                             metadata["TX_TIME_ALG"], observable_lst[0])

            self._data[sat] = geometry

        for sat in _to_remove:
            self.remove(sat)

    def get_unit_line_of_sight(self, sat):
        """
        Computes the line of sight vector between the receiver and the satellite, used in the PVT geometry matrix.

        Args:
            sat (src.data_types.data_types.Satellite.Satellite) : satellite to compute the LOS vector

        Return:
            list [float, float, float] : Line of sight for [x, y, z] axis of ECEF frame
        """
        return self.get("los", sat)

    def __str__(self):
        return str(self._data)

    def __len__(self):
        return len(self.get_satellites())
