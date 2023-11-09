import numpy as np

from src.io.config import config_dict
from src.io.config.enums import EnumOnOff, EnumIono, EnumTropo
from src.models.frames import cartesian2geodetic
from src.models.gnss_obs.iono_klobuchar import iono_klobuchar
from src.models.gnss_obs.iono_ntcmg import NTCMG
from src.models.gnss_obs.clock_obs import broadcast_clock, nav_sat_clock_correction
from src.data_types.gnss.observation import Observation
from src import constants
from src.models.gnss_obs.tropo_saastamoinen import tropo_saastamoinen


class ObservationReconstruction:
    def __init__(self, system_geometry, metadata, state, nav_header):
        self._metadata = metadata
        self._state = state
        self._system_geometry = system_geometry
        self._nav_header = nav_header

    def compute(self, nav_message, sat, epoch, datatype):
        iono = 0.0
        tropo = 0.0

        receiver_position = self._system_geometry.get("receiver_position", sat)
        az = self._system_geometry.get("az", sat)  # satellite azimuth from receiver
        el = self._system_geometry.get("el", sat)  # satellite elevation from receiver
        [lat, long, height] = cartesian2geodetic(*receiver_position)  # user lat, long and height

        # true range
        true_range = self._system_geometry.get("true_range", sat)

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

        # ionosphere
        if not config_dict.is_iono_free():
            if self._metadata["IONO"][sat.sat_system] == EnumIono.KLOBUCHAR:
                iono = iono_klobuchar(lat, long, el, az, self._nav_header.iono_corrections["GPSA"],
                                      self._nav_header.iono_corrections["GPSB"], epoch.gnss_time[1],
                                      datatype.freq)
            elif self._metadata["IONO"][sat.sat_system] == EnumIono.NTCMG:
                ut1 = epoch.change_scale("UT1")
                iono = NTCMG.calculate_ionospheric_contribution(ut1, lat, long, el, az,
                                                                self._nav_header.iono_corrections["GAL"], datatype.freq)

        # troposphere
        if self._metadata["TROPO"][sat.sat_system] == EnumTropo.SAASTAMOINEM:
            tropo = tropo_saastamoinen(height, lat, epoch.doy, el)

        # finally, construct obs
        obs = true_range - dt_sat * constants.SPEED_OF_LIGHT + iono + tropo

        return Observation(datatype, obs)

    def get_unit_line_of_sight(self, sat):
        return self._system_geometry.get_unit_line_of_sight(sat)

    def get_weight(self, sat):
        # TODO: need to add here the user defined sigmas as a multiplication factor
        # "obs_std", and can add the possibility of this mask as well.
        elevation = self._system_geometry.get("el", sat)
        sigma_elevation = np.e ** (-elevation)
        w = (1 / sigma_elevation) ** 2

        return w
