""" Satellite Geometry Module.
This module computes several geometry-related data for the reconstruction of GNSS observation equations.
"""
from datetime import timedelta
import numpy as np

from src.constants import SPEED_OF_LIGHT
from src.io.config import EnumTransmissionTime, config_dict, EnumAlgorithmPNT
from src.models.frames import enu2azel, ecef2enu, cartesian2geodetic, dcm_e_i
from src.data_mng import Container
from src.common_log import MODEL_LOG, get_logger
from src.spicepy_wrapper import compute_sun_pos
from src.models.gnss_models import gnss_attitude, compute_tx_time

_warning_cache = set()


class SatelliteGeometry(Container):
    """
    `SatelliteGeometry` class, inherits from the `Container` class.
    This class is a container that stores geometry-related data for the reconstruction of GNSS observation
    equations. It stores the following data (for a single epoch and satellite):

    Attributes:
        transit_time(float): transmission time (time that signal takes to travel from the satellite to the receiver)
        time_emission(src.data_types.date.Epoch): epoch of signal transmission
            (time_emission = time_reception - transit_time)
        time_reception(src.data_types.date.Epoch): epoch of signal reception
        true_range(float): geometric range
        az(float): satellite azimuth as seen from the receiver frame
        el(float): satellite elevation as seen from the receiver frame
        satellite_position(numpy.ndarray): satellite position in ECEF frame
        satellite_velocity(numpy.ndarray): satellite velocity in ECEF frame
        dt_rel_correction(float): relativistic correction of satellite clock
        los(numpy.ndarray): line of sight vector from receiver to satellite
        tropo_map_wet(float): map of tropospheric wet component
        drift_rel_correction(float): relativistic clock drift correction of satellite clock
        dcm_b_e(numpy.ndarray): DCM from satellite body-fixed frame to ECEF frame (at TX time)
        nadir_sat(float): nadir angle of the satellite. This is the angle between the direction from satellite toward
            the center of the earth and the LOS vector
        azimuth_sat(float): azimuth angle of the satellite. This is the angle between the LOS vector projected in the
            XY-plane and the +Y-axis, measured clockwise toward +X when looking in the direction of -Z (toward deep
            space). Axes here are referred to the satellite body-fixed frame.
    """
    __slots__ = ["transit_time", "time_emission", "time_reception", "true_range", "az", "el",
                 "satellite_position", "satellite_velocity", "dt_rel_correction", "los", "tropo_map_wet",
                 "drift_rel_correction", "dcm_b_e", "nadir_sat", "azimuth_sat"]

    def __init__(self):
        """ Base Constructor with no arguments. The attributes are filled in the `compute` method.  """
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
        self.drift_rel_correction = 0
        self.los = None
        self.dcm_b_e = None
        self.nadir_sat = 0
        self.azimuth_sat = 0

    def __str__(self):
        _allAttrs = ""
        for atr in self.__slots__:
            _allAttrs += atr + "=" + str(getattr(self, atr)) + ", "
        return f'{type(self).__name__}({_allAttrs[0:-2]})'

    def __repr__(self):
        return str(self)

    def compute(self, sat, epoch, state, constellation, compute_tx, PR_obs, sat_orbits, sat_clocks):
        """
        compute satellite-related quantities (tropo, iono, transmission time, etc.) to be used in the reconstruction
        of GNSS observables, for a given satellite.

        Args:
            sat(src.data_types.gnss.Satellite): the satellite to compute the geometry data
            epoch (src.data_types.date.Epoch) : epoch under evaluation
            state (src.data_mng.gnss.state_space.GnssStateSpace) : state
            constellation (src.data_types.gnss.Constellation) : constellation
            compute_tx (function) : function to compute the transmission time
            PR_obs (src.data_types.gnss.observation.Observation) : pseudorange observation to use in some computations
            sat_orbits(src.data_mng.gnss.sat_orbit_data.SatelliteOrbits): `SatelliteOrbits` object with orbit data
            sat_clocks(src.data_mng.gnss.sat_clock_data.SatelliteClocks): `SatelliteClocks` object with clock data


        Reference:
            [1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
                Springer Cham, 2017
        """
        rec_pos = state.position
        time_correction = sat_clocks.nav_data.header.time_correction if sat_clocks.nav_data is not None else None
        rec_bias = state.get_clock_bias(constellation, time_correction) / SPEED_OF_LIGHT  # receiver clock bias in [s]

        # get reception time in GNSS time system ( T_GNSS = T_receiver - t_(receiver_bias) )
        time_reception = epoch + timedelta(seconds=-rec_bias)

        # algorithm to compute Transmission time (in GNSS time)
        time_emission, transit = compute_tx_time(model=compute_tx, r_receiver=rec_pos,
                                                 t_reception=epoch, v_receiver=state.velocity,
                                                 dt_receiver=rec_bias, pseudorange_obs=PR_obs, sat_orbits=sat_orbits,
                                                 sat_clocks=sat_clocks, sat=sat)

        # get satellite position and velocity at RX ECEF frame
        p_sat, v_sat, dt_relative, drift_relative, r_sat_tx = sat_orbits.compute_orbit_at_rx_time(sat, time_emission,
                                                                                                  transit)

        # Compute geometrical range
        if compute_tx != EnumTransmissionTime.NAPEOS:
            true_range = np.linalg.norm(p_sat - rec_pos)
        else:
            true_range = transit * SPEED_OF_LIGHT

        # TODO: apply shapiro correction here (Eq. (19.14) of ref [1])

        # get satellite elevation and azimuth angles (from the receiver), ENU frame
        lat, long, h = cartesian2geodetic(rec_pos[0], rec_pos[1], rec_pos[2])
        enu_coord = ecef2enu(p_sat[0], p_sat[1], p_sat[2], lat, long, h)
        az, el = enu2azel(*enu_coord)

        # line of sight (Eq. (21.21) of [1])
        los = np.array([(rec_pos[i] - p_sat[i]) / true_range for i in (0, 1, 2)])

        dcm_b_e = None
        nadir_sat = None
        azimuth_sat = None
        if config_dict.get("inputs", "cspice_kernels", "enable"):
            pvt_alg = config_dict.get("gnss_alg")
            if pvt_alg == EnumAlgorithmPNT.SPS:
                if "sps" not in _warning_cache:
                    log = get_logger(MODEL_LOG)
                    log.warning(f"Skipping GNSS Attitude computation in SPS Mode.")
                    _warning_cache.add("sps")
            else:
                try:
                    # GNSS Attitude: rotation matrix from satellite body-fixed frame to ECEF frame at RX time
                    _R = dcm_e_i(-transit)
                    sun_pos_tx = compute_sun_pos(time_emission)
                    dcm_tx = gnss_attitude.gnss_attitude(r_sat_tx, sun_pos_tx)
                    dcm_b_e = _R @ dcm_tx

                    # compute satellite nadir and azimuth angles (in the satellite body frame)
                    los_body = dcm_b_e.T @ los  # los in satellite body-fixed frame
                    azimuth_sat = np.arctan2(los_body[0], los_body[1])
                    if azimuth_sat < 0:
                        azimuth_sat += 2 * np.pi  # Ensure range [0, 2pi]
                    e_z = -p_sat / np.linalg.norm(p_sat)
                    nadir_sat = np.arccos(np.dot(los, e_z) / np.linalg.norm(los))
                    # NOTE: The azimuth angle of the receiver as seen from the satellite is the angle between the
                    # projected LOS vector in the XY-plane and the +Y-axis, measured clockwise toward +X when looking in
                    # the direction of -Z (toward deep space).
                except Exception as e:
                    if sat not in _warning_cache:
                        log = get_logger(MODEL_LOG)
                        log.warning(f"Failed to compute satellite attitude for sat {sat} due to: {e}")
                        _warning_cache.add(sat)
        else:
            if sat not in _warning_cache:
                log = get_logger(MODEL_LOG)
                log.info(f"Not computing satellite attitude for sat {sat} because CSpice is disabled.")
                _warning_cache.add(sat)

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
        self.drift_rel_correction = drift_relative
        self.dcm_b_e = dcm_b_e
        self.nadir_sat = nadir_sat
        self.azimuth_sat = azimuth_sat


class SystemGeometry:
    """ System Geometry class.
    This class is serves as a dataframe to store `SatelliteGeometry` objects for all available satellites.
    """
    def __init__(self, obs_data, sat_clocks, sat_orbits, phase_center, sat_bias):
        """
        Constructor of the SystemGeometry. This is a container that stores data (instances of `SatelliteGeometry`)
        for all available satellites and for a single epoch. This data is useful in the reconstruction equations
        of the GNSS observables.

        Args:
            obs_data (src.data_mng.gnss.observation_data.EpochData) : instance of `EpochData` (GNSS observable database
                for a single epoch)
            sat_clocks(src.data_mng.gnss.sat_clock_data.SatelliteClocks): `SatelliteClocks` object with clock database
            sat_orbits(src.data_mng.gnss.sat_orbit_data.SatelliteOrbits): `SatelliteOrbits` object with orbit database
            phase_center(src.data_mng.gnss.phase_center_mng.PhaseCenterManager): `PhaseCenterManager` object with
                phase center data
            sat_bias(src.data_mng.gnss.bias_manager.BiasManager): manager of satellite code and phase biases

        """
        self._data = dict.fromkeys(obs_data.get_satellites())
        self.obs_data = obs_data
        self.sat_clocks = sat_clocks
        self.sat_orbits = sat_orbits
        self.phase_center = phase_center
        self.sat_bias = sat_bias

    def _clean(self):
        vSats = self.get_satellites()
        self._data.clear()
        self._data = dict.fromkeys(vSats)

    def items(self):
        return self._data.items()

    def get_satellites(self):
        """ Return list of available satellites.

        Returns:
            list[src.data_types.gnss.Satellite]: list of all available satellites.
        """
        return list(self._data.keys())

    def remove(self, sat):
        """ Remove a satellite from the internal dataframe dict """
        if sat in self._data:
            self._data.pop(sat)

    def get(self, attribute, sat):
        """
        Returns the required attribute for the given satellite. The attribute is one of the attributes of the
        :py:class:`SatelliteGeometry` class.

        Args:
             attribute(str): the string name of the attribute to fetch
             sat(src.data_types.gnss.Satellite): the satellite to fetch the data

        Returns:
            Any: The value of the queried attribute. The method returns None if the satellite is not available in
                the internal dataframe dict.
        """
        if sat in self._data:
            return getattr(self._data[sat], attribute)
        return None

    def set(self, attribute, value, sat):
        """
        Sets the provided attribute for the given satellite. The attribute is one of the attributes of the
        :py:class:`SatelliteGeometry` class.

        Args:
             attribute(str): the string name of the attribute to set
             value(Any): the value of the attribute
             sat(src.data_types.gnss.Satellite): the satellite to fill in the attribute data

        Raises:
            KeyError: a `KeyError` exception is raised if the provided satellite is unavailable
        """
        if sat in self._data:
            setattr(self._data[sat], attribute, value)
        else:
            raise KeyError(f"Satellite {sat} not in SystemGeometry data structure")

    def compute(self, epoch, state, metadata):
        """
        compute satellite-related quantities (tropo, iono, transmission time, etc.) to be used in the
        Observation reconstruction models

        Args:
            epoch (src.data_types.date.Epoch) : observation epoch for the computations
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
                log = get_logger(MODEL_LOG)
                log.warning(f"Failed to compute geometry for sat {sat} at epoch {epoch} due to:"
                            f" no valid observations for this sat and epoch")
                _to_remove.append(sat)
                continue

            # compute geometry for this satellite
            try:
                geometry.compute(sat, epoch, state, sat.sat_system, metadata["TX_TIME_ALG"], observable_lst[0],
                                 self.sat_orbits, self.sat_clocks)
                self._data[sat] = geometry
            except Exception as e:
                log = get_logger(MODEL_LOG)
                log.warning(f"Failed to compute geometry for sat {sat} at epoch {epoch} due to:"
                            f" {e}")
                _to_remove.append(sat)
                continue

        for sat in _to_remove:
            self.remove(sat)

    def get_unit_line_of_sight(self, sat):
        """
        Computes the line of sight vector between the receiver and the satellite, used in the PVT geometry matrix.

        Args:
            sat (src.data_types.gnss.Satellite) : satellite to compute the LOS vector

        Returns:
            numpy.ndarray : Line of sight for [x, y, z] axis of ECEF frame
        """
        return self.get("los", sat)

    def __str__(self):
        return str(self._data)

    def __len__(self):
        """ Returns the number of available satellites """
        return len(self.get_satellites())
