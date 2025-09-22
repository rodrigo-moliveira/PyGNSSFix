""" IONEX ionospheric module from Global Ionosphere Maps """

import math
from src.constants import DEG2RAD


class IonoIonex:
    """
    Ionospheric model from precise Global Ionosphere Maps (GIM) in IONEX format.

    More information can be found in :
        [1]: http://ftp.aiub.unibe.ch/ionex/draft/ionex11.pdf
        [2] RTKLIB ver. 2.4.2 Manual, Section E.5, p. 151-152. (https://www.rtklib.com/prog/manual_2.4.2.pdf)
    """
    def __init__(self, gim_data):
        """
        Constructor of the IonoIonex model.

        Args:
            gim_data(src.data_mng.gnss.global_iono_map.GlobalIonoMap): Global Ionosphere Map data
        """
        self.gim_data = gim_data

    def __str__(self):
        """ String representation of the model """
        return "IONEX Ionospheric Model"

    def compute(self, ut1, user_lat, user_long, sv_el, sv_az, freq):
        """

        Main function of the ionosphere from the IONEX GIM model.

        Args:
            ut1(src.data_types.date.Epoch): required epoch in UT1 to compute the iono delay
            user_lat(float): user latitude in [rad]
            user_long(float): user longitude in [rad]
            sv_el(float): satellite elevation in [rad]
            sv_az(float): satellite azimuth in [rad]
            freq(src.data_types.gnss.data_type.DataType): Frequency band required for the iono delay
        Returns:
            float : ionosphere correction for the selected frequency [m]
        """
        # compute pierce point
        ipp_latlon, slant_factor = self.calculate_pierce_point_lat_lon([user_lat, user_long], [sv_az, sv_el], )

        # interpolate VTEC from GIM data for the pierce point and epoch
        vtec = self.gim_data.interpolate(ut1, ipp_latlon[0], ipp_latlon[1], sun_fixed=True)

        # compute the ionospheric delay, in meters
        iono = slant_factor * 40.30E16 / freq.freq_value / freq.freq_value * vtec

        return iono

    def calculate_pierce_point_lat_lon(self, pos, azel):
        """
        Compute the ionospheric pierce point latitude and longitude.

        This implementation is based on the algorithm described in reference [2] (RTKLIB manual).

        Args:
            pos(list)  : list of two elements [latitude (rad), longitude (rad)]
            azel(list) : list of two elements [azimuth (rad), elevation angle (rad)]

        Returns:
            list: [latitude (rad), longitude (rad)] coordinates of the ionospheric pierce point
            float: slant factor
        """
        Re = self.gim_data.header.base_radius * 1E3  # Earth radius in meters
        h_iono = self.gim_data.header.hgt[0] * 1E3  # Ionospheric reference height in meters

        cos_az = math.cos(azel[0])
        rp = Re / (Re + h_iono) * math.cos(azel[1])
        ap = math.pi / 2.0 - azel[1] - math.asin(rp)
        sin_ap = math.sin(ap)
        tan_ap = math.tan(ap)

        # Calculate ionospheric pierce point latitude
        ipp_lat = math.asin(
            math.sin(pos[0]) * math.cos(ap) + math.cos(pos[0]) * sin_ap * cos_az
        )

        # Check conditions to adjust longitude
        if ((pos[0] > 70*DEG2RAD and tan_ap * cos_az > math.tan(math.pi / 2.0 - pos[0])) or
                (pos[0] < -70*DEG2RAD and -tan_ap * cos_az > math.tan(math.pi / 2.0 + pos[0]))):
            ipp_lon = pos[1] + math.pi - math.asin(sin_ap * math.sin(azel[0]) / math.cos(ipp_lat))
        else:
            ipp_lon = pos[1] + math.asin(sin_ap * math.sin(azel[0]) / math.cos(ipp_lat))

        ipp_latlon = [ipp_lat, ipp_lon]

        # compute slant factor
        slant_factor = 1.0 / math.sqrt(1.0 - rp ** 2)

        return ipp_latlon, slant_factor
