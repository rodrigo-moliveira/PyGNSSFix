import numpy as np

from src.data_types.gnss.data_type import DataType
from src.io.config import config_dict
from src.io.config.enums import EnumOnOff, EnumIono, EnumModel
from src.models.frames import cartesian2geodetic
from src.models.gnss_obs.ionosphere.iono_klobuchar import iono_klobuchar
from src.models.gnss_obs.ionosphere.iono_ntcmg import NTCMG
from src.models.gnss_obs.clock_obs import broadcast_clock, nav_sat_clock_correction
from src.data_types.gnss.observation import Observation
from src import constants


class PseudorangeReconstructor:
    def __init__(self, system_geometry, metadata, state, nav_header):
        self._metadata = metadata
        self._state = state
        self._system_geometry = system_geometry
        self._nav_header = nav_header

    def compute(self, nav_message, sat, epoch, datatype):
        iono = 0.0

        az = self._system_geometry.get("az", sat)  # satellite azimuth from receiver
        el = self._system_geometry.get("el", sat)  # satellite elevation from receiver
        [lat, long, height] = cartesian2geodetic(*self._state.position)  # user lat, long and height

        # true range
        true_range = self._system_geometry.get("true_range", sat)

        # user clock in meters (with proper ISB applied, if necessary)
        dt_rec = self._state.get_clock_bias(sat.sat_system, self._nav_header.time_correction) * constants. \
            SPEED_OF_LIGHT

        # satellite clock
        dt_sat, _ = broadcast_clock(nav_message.af0,
                                    nav_message.af1,
                                    nav_message.af2,
                                    nav_message.toc.gnss_time[1],
                                    self._system_geometry.get("time_emission", sat).gnss_time[1])

        # correct satellite clock for relativistic corrections
        if self._metadata["REL_CORRECTION"] == EnumOnOff.ENABLED:
            dt_sat += self._system_geometry.get("dt_rel_correction", sat)

        # correct satellite clock for BGDs
        dt_sat = nav_sat_clock_correction(dt_sat, datatype, nav_message)

        # ionosphere (a-priori correction)
        if not DataType.is_iono_free_code(datatype) and not DataType.is_iono_free_smooth_code(datatype):
            if self._metadata["IONO"][sat.sat_system] == EnumIono.KLOBUCHAR:
                iono = iono_klobuchar(lat, long, el, az, self._nav_header.iono_corrections["GPSA"],
                                      self._nav_header.iono_corrections["GPSB"], epoch.gnss_time[1],
                                      datatype.freq)
            elif self._metadata["IONO"][sat.sat_system] == EnumIono.NTCMG:
                ut1 = epoch.change_scale("UT1")
                iono = NTCMG.calculate_ionospheric_contribution(ut1, lat, long, el, az,
                                                                self._nav_header.iono_corrections["GAL"], datatype.freq)
        # iono estimated correction dI
        dI = 0.0
        if self._metadata["MODEL"][sat.sat_system] == EnumModel.DUAL_FREQ:
            try:
                factor = (self._metadata["CODES"][sat.sat_system][0].freq.freq_value /
                          datatype.freq.freq_value) ** 2
                dI = factor * self._state.iono[sat]
            except KeyError:
                pass

        # troposphere
        tropo, map_wet = self._metadata["TROPO"].compute_tropo_delay(lat, long, height, el, epoch, self._state)
        self._system_geometry.set("tropo_map_wet", map_wet, sat)

        # finally, construct obs
        obs = true_range + dt_rec - dt_sat * constants.SPEED_OF_LIGHT + iono + tropo + dI
        return Observation(datatype, obs)

    def get_unit_line_of_sight(self, sat):
        return self._system_geometry.get_unit_line_of_sight(sat)

    def get_obs_std(self, sat, datatype):
        elevation_mask = config_dict.get("model", sat.sat_system, "elevation_mask")
        el_std = 1.0
        if elevation_mask:
            elevation = self._system_geometry.get("el", sat)
            el_std = np.e ** (-elevation)
        try:
            obs_std = config_dict.get_obs_std()[sat.sat_system][datatype]
        except KeyError:
            # TODO add logger message
            obs_std = 1.0
        return el_std * obs_std


class DopplerReconstructor:
    def __init__(self, system_geometry, metadata, state):
        self._metadata = metadata
        self._state = state
        self._system_geometry = system_geometry

    def compute(self, nav_message, sat, epoch, datatype):
        return 0.0

    def get_unit_line_of_sight(self, sat):
        return self._system_geometry.get_unit_line_of_sight(sat)
