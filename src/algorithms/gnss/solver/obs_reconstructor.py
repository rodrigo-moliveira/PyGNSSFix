from src.io.config import config_dict
from src.io.config.enums import EnumOnOff, EnumIono, EnumTropo
from src.models.frames import cartesian2geodetic
from src.models.observation.atmosphere_obs import iono_klobuchar, tropo_saastamoinen
from src.models.observation.clock_obs import gps_broadcast_clock
from src.data_types.gnss.data_type import L1, DataType
from src.data_types.gnss.observation import Observation
from src import constants


class ObservationReconstruction:
    def __init__(self, system_geometry, tropo, iono, nav_header, relativistic_correction):
        self._model = {"tropo": tropo,  # a priori troposphere model
                       "iono": iono,  # a priori ionosphere model
                       "relativistic_correction": relativistic_correction
                       }
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
        dt_sat, _ = gps_broadcast_clock(nav_message.af0,
                                        nav_message.af1,
                                        nav_message.af2,
                                        nav_message.toc.gps_time[1],
                                        self._system_geometry.get("time_emission", sat).gps_time[1])

        # correct satellite clock for relativistic corrections
        if self._model["relativistic_correction"] == EnumOnOff.ENABLED:
            dt_sat += self._system_geometry.get("dt_rel_correction", sat)

        # correct for satellite clock for TGD (TGD is 0 for iono free observables (in GPS SPS only...))
        if not DataType.is_iono_free_smooth_code(datatype) and not DataType.is_iono_free_code(datatype):
            tgd = (L1.freq_value / datatype.freq.freq_value) ** 2 * self._system_geometry.get("tgd", sat)
            dt_sat -= tgd

        # ionosphere
        if not config_dict.is_iono_free():
            if self._model["iono"][sat.sat_system] == EnumIono.KLOBUCHAR:
                time_reception = self._system_geometry.get("time_reception", sat)
                iono = iono_klobuchar(lat, long, el, az, self._nav_header.iono_corrections["GPSA"],
                                      self._nav_header.iono_corrections["GPSB"], time_reception.gps_time[1])

                # fix I for non L1 users (L2 or L5)
                if datatype.freq != L1:
                    iono = (L1.freq_value / datatype.freq.freq_value) ** 2 * iono

        # troposphere
        if self._model["tropo"][sat.sat_system] == EnumTropo.SAASTAMOINEM:
            tropo = tropo_saastamoinen(height, lat, epoch.doy, el)

        # finally, construct obs
        obs = true_range - dt_sat * constants.SPEED_OF_LIGHT + iono + tropo

        return Observation(datatype, obs)
