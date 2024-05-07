from datetime import timedelta
import numpy as np

from src.models.frames import enu2azel, ecef2enu, cartesian2geodetic
from src.data_mng import Container
from src.errors import TimeSeriesError
from src.models.gnss_models.navigation import compute_tx_time


class SatelliteGeometry(Container):
    """
    SatelliteGeometry class. This class is a container that stores geometry-related data for the
    reconstruction of GNSS observation equations. It stores the following data (for a single epoch and satellite):
        * `transit_time`: transmission time (time that signal takes to travel from the satellite to the receiver)
        * `time_emission`: epoch of signal transmission (time_emission = time_reception - transit_time)
        * `time_reception`: epoch of signal reception
        * `true_range`: geometric range
        * `az`: satellite azimuth as seen from the receiver frame
        * `el`: satellite elevation as seen from the receiver frame
        * `satellite_position`: satellite position in ECEF frame
        * `satellite_velocity`: satellite velocity in ECEF frame
        * `dt_rel_correction`: relativistic correction of satellite clock
        * `los`: line of sight vector from receiver to satellite
        * `tropo_map_wet`: map of tropospheric wet component
    """
    __slots__ = ["transit_time", "time_emission", "time_reception", "true_range", "az", "el",
                 "satellite_position", "satellite_velocity", "dt_rel_correction", "los", "tropo_map_wet"]

    def __init__(self):
        super().__init__()
        self.transit_time = 0
        self.time_emission = None
        self.time_reception = None
        self.true_range = 0
        self.az = 0
        self.el = 0
        self.tropo_map_wet = 0
        self.satellite_position = None
        self.satellite_velocity = None
        self.dt_rel_correction = 0
        self.los = None

    def __str__(self):
        _allAttrs = ""
        for atr in self.__slots__:
            _allAttrs += atr + "=" + str(getattr(self, atr)) + ", "
        return f'{type(self).__name__}({_allAttrs[0:-2]})'

    def __repr__(self):
        return str(self)

    def compute(self, sat, epoch, state, constellation, compute_tx, PR_obs, sat_orbits, sat_clocks):
        """
        compute satellite-related quantities (tropo, iono, transmission time, etc.) to be used in the PVT gnss_models
        reconstruction equation, for a given satellite.

        Args:
            sat(src.data_types.gnss.satellite.Satellite): the satellite to compute the geometry data
            epoch (src.data_types.basics.Epoch.Epoch) : epoch under evaluation
            state (State) : state
            constellation (Constellation) : constellation
            compute_tx (function) : function to compute the transmission time
            PR_obs (src.data_types.data_types.Observation.Observation) : Code gnss_models to use in some computations
            sat_orbits(src.data_mng.gnss.sat_orbit_data.SatelliteOrbits): `SatelliteOrbits` object with orbit data
            sat_clocks(src.data_mng.gnss.sat_clock_data.SatelliteClocks): `SatelliteClocks` object with clock data
        """
        # TODO: parei aqui..
        # TODO: nav_messages are no longer required.
        #   nav_header I need to evaluate, but most likely not, it is included in the sat_clocks.nav
        #   if it is None or empty, then set to 0 or something like that
        #   if the sat is missing, in sat_clocks or orbits, an exception must be raised...
        rec_pos = state.position
        rec_bias = state.get_clock_bias(constellation, nav_header.time_correction)

        # get reception time in GNSS time system ( T_GNSS = T_receiver - t_(receiver_bias) )
        time_reception = epoch + timedelta(seconds=-rec_bias)

        # algorithm to compute Transmission time (in GNSS time)
        time_emission, transit = compute_tx_time(model=compute_tx, r_receiver=rec_pos, t_reception=epoch,
                                                 dt_receiver=rec_bias, pseudorange_obs=PR_obs, sat_orbits=sat_orbits,
                                                 sat_clocks=sat_clocks, sat=sat)

        # get and satellite position at RX ECEF frame
        p_sat, v_sat, dt_relative = sat_orbits.compute_sat_nav_position_dt_rel(sat, time_emission, transit)

        # compute true range
        true_range = np.linalg.norm(p_sat - rec_pos)

        # TODO: apply shapiro correction here

        # get satellite elevation and azimuth angles (from the receiver), ENU frame
        lat, long, h = cartesian2geodetic(rec_pos[0], rec_pos[1], rec_pos[2])
        enu_coord = ecef2enu(p_sat[0], p_sat[1], p_sat[2], lat, long, h)
        az, el = enu2azel(*enu_coord)

        # line of sight
        los = np.array([(rec_pos[i] - p_sat[i]) / true_range for i in (0, 1, 2)])

        # save results in container
        self.transit_time = transit
        self.time_emission = time_emission
        self.time_reception = time_reception
        self.true_range = true_range
        self.az = az
        self.el = el
        self.satellite_position = p_sat
        self.satellite_velocity = v_sat
        self.dt_rel_correction = dt_relative
        self.los = los


class SystemGeometry:
    def __init__(self, obs_data, sat_clocks, sat_orbits):
        """
        Constructor of the SystemGeometry. This is a container that stores data (instances of `SatelliteGeometry`)
        for all available satellites and for a single epoch. This data is useful in the reconstruction equations
        of the GNSS observables.

        Args:
            obs_data (src.data_mng.gnss.observation_data.EpochData) : instance of `EpochData` (GNSS observable database
                for a single epoch)
            sat_clocks(src.data_mng.gnss.sat_clock_data.SatelliteClocks): `SatelliteClocks` object with clock database
            sat_orbits(src.data_mng.gnss.sat_orbit_data.SatelliteOrbits): `SatelliteOrbits` object with orbit database
        """
        self._data = dict.fromkeys(obs_data.get_satellites())
        self.obs_data = obs_data
        self.sat_clocks = sat_clocks
        self.sat_orbits = sat_orbits

    def _clean(self):
        vSats = self.get_satellites()
        self._data.clear()
        self._data = dict.fromkeys(vSats)

    def items(self):
        return self._data.items()

    def get_satellites(self):
        """Return list of available satellites"""
        return list(self._data.keys())

    def remove(self, sat):
        """Remove a satellite from the internal dict"""
        if sat in self._data:
            self._data.pop(sat)

    def get(self, attribute, sat):
        if sat in self._data:
            return getattr(self._data[sat], attribute)
        return None

    def set(self, attribute, value, sat):
        if sat in self._data:
            setattr(self._data[sat], attribute, value)
        else:
            raise TimeSeriesError(f"Satellite {sat} not in SystemGeometry data structure")

    def compute(self, epoch, state, metadata):
        """
        compute satellite-related quantities (tropo, iono, transmission time, etc.) to be used in the
        Observation reconstruction models

        Args:
            epoch (src.data_types.date.date.Epoch) : observation epoch for the computations
            state (src.data_mng.gnss.state_space.GnssStateSpace) : current GNSS state vector
            metadata (dict) : dictionary with metadata information for the models
        """
        self._clean()
        _to_remove = []
        sat_list = self.get_satellites()

        for sat in sat_list:
            geometry = SatelliteGeometry()

            # fetch pseudorange observables for this satellite at epoch
            observable_lst = self.obs_data.get_code_observables(sat)
            if len(observable_lst) == 0:
                _to_remove.append(sat)
                continue

            # compute geometry for this satellite
            try:
                geometry.compute(sat, epoch, state, sat.sat_system, metadata["TX_TIME_ALG"], observable_lst[0],
                                 self.sat_orbits, self.sat_clocks)
                self._data[sat] = geometry
            except:  # TODO: specify which exception
                _to_remove.append(sat)
                continue

        for sat in _to_remove:
            self.remove(sat)

    def get_unit_line_of_sight(self, sat):
        """
        Computes the line of sight vector between the receiver and the satellite, used in the PVT geometry matrix.

        Args:
            sat (src.data_types.data_types.Satellite.Satellite) : satellite to compute the LOS vector

        Return:
            numpy.ndarray : Line of sight for [x, y, z] axis of ECEF frame
        """
        return self.get("los", sat)

    def __str__(self):
        return str(self._data)

    def __len__(self):
        return len(self.get_satellites())
