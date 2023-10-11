from datetime import timedelta

from .clock_obs import compute_tx_time
from .ephemeride_propagator import EphemeridePropagator
from src.data_types.gnss.data_type import DataType, L1
from ..frames.frames import enu2azel, ecef2enu, cartesian2geodetic
from ...data_mng.container import Container
from ...errors import NonExistentObservable


class SatelliteGeometry(Container):
    __slots__ = ["transit_time", "receiver_position",
                 "time_emission", "time_reception",
                 "true_range", "az", "el", "satellite_position",
                 "dt_rel_correction"]

    def __init__(self):
        super().__init__()
        self.transit_time = 0
        self.time_emission = None
        self.time_reception = None
        self.true_range = 0
        self.az = 0
        self.el = 0
        self.satellite_position = None
        self.receiver_position = None
        self.dt_rel_correction = 0

    def __str__(self):
        _allAttrs = ""
        for atr in self.__slots__:
            _allAttrs += atr + "=" + str(getattr(self, atr)) + ", "
        return f'{type(self).__name__}({_allAttrs[0:-2]})'

    def __repr__(self):
        return str(self)

    def compute(self, rec_pos, epoch, rec_bias, nav_message, compute_tx, PR_obs, relativistic_correction):
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
            relativistic_correction (bool) : whether to compute relativistic correction
        """
        # get reception time in GNSS time system ( T_GNSS = T_receiver - t_(receiver_bias) )
        time_reception = epoch + timedelta(seconds=-rec_bias)

        # algorithm to compute Transmission time (in GNSS time)
        time_emission, transit = compute_tx_time(model=compute_tx,
                                                 r_receiver=rec_pos,
                                                 t_reception=epoch,
                                                 dt_receiver=rec_bias,
                                                 nav_message=nav_message,
                                                 pseudorange_obs=PR_obs)

        # get rho_0 and satellite position at RX ECEF frame
        p_sat, true_range, dt_relative = EphemeridePropagator.get_sat_position_and_true_range(nav_message,
                                                                                              time_emission.gps_time,
                                                                                              transit, rec_pos,
                                                                                              relativistic_correction)

        # get satellite elevation and azimuth angles (from the receiver), ENU frame
        lat, long, h = cartesian2geodetic(rec_pos[0], rec_pos[1], rec_pos[2])
        enu_coord = ecef2enu(p_sat[0], p_sat[1], p_sat[2], lat, long, h)
        az, el = enu2azel(*enu_coord)

        # save results in container
        self.transit_time = transit
        self.time_emission = time_emission
        self.time_reception = time_reception
        self.true_range = true_range
        self.az = az
        self.el = el
        self.satellite_position = p_sat
        self.receiver_position = rec_pos
        self.dt_rel_correction = dt_relative


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

    def compute(self, epoch, receiver_position, receiver_clock, compute_tx,
                main_datatype, relativistic_correction):
        """
        compute satellite-related quantities (tropo, iono, transmission time, etc.) to be used in the PVT gnss_obs
        reconstruction equation, for all available satellites.

        Args:
            epoch (src.data_types.basics.Epoch.Epoch) : epoch under evaluation
            receiver_position (src.data_types.state_space.statevector.Position) : Receiver position
            receiver_clock (float) : Receiver clock bias
            compute_tx (function) : function to compute the transmission time
            main_datatype (DataType) : Code gnss_obs to use in some computations
            relativistic_correction (bool) : whether to compute relativistic correction
        """
        self._clean()
        _to_remove = []
        sat_list = self.get_satellites()

        for sat in sat_list:
            geometry = SatelliteGeometry()

            # fetch navigation message for this satellite
            nav_message = self.nav_data.get_sat_data_for_epoch(sat, epoch)

            # fetch pseudorange gnss_obs for this satellite at epoch (used in the compute_TX_time algorithm)
            try:
                main_observation = self.obs_data.get_observable(sat, main_datatype)
            except NonExistentObservable:
                # could not retrieve this observable
                _to_remove.append(sat)
                continue

            geometry.compute(receiver_position, epoch, receiver_clock, nav_message,
                             compute_tx, main_observation, relativistic_correction)

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

        receiver = self.get("receiver_position", sat)
        satellite = self.get("satellite_position", sat)
        true_range = self.get("true_range", sat)

        return [(receiver[i] - satellite[i]) / true_range for i in (0, 1, 2)]

    def __str__(self):
        return str(self._data)

    def __len__(self):
        return len(self.get_satellites())
