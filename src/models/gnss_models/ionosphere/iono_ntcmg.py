"""NTCM-G ionospheric module for Galileo
"""

import numpy as np
from dataclasses import dataclass

from math import sqrt, asin, cos, sin, degrees, radians, exp
from src.constants import PI

# Constants for this model (note that the Earth radius has this specific value for this model) 
Re_km = 6371
hI_km = 450
lat_north_pole_deg = 79.74
lon_north_pole_deg = -71.78


# Values of NTCM-G coefficients
@dataclass
class NTCMG_coefs:
    k1: float = 0.92519
    k2: float = 0.16951
    k3: float = 0.00443
    k4: float = 0.06626
    k5: float = 0.00899
    k6: float = 0.21289
    k7: float = -0.15414
    k8: float = -0.38439
    k9: float = 1.14023
    k10: float = 1.20556
    k11: float = 1.41808
    k12: float = 0.13985


class NTCMG:
    """
    NTCM-G ionospheric model for Galileo users

    More information can be found in :
        https://www.gsc-europa.eu/sites/default/files/NTCM-G_Ionospheric_Model_Description_-_v1.0.pdf
    """

    @staticmethod
    def calculate_azpar(effec_iono: np.array):
        """Compute Azpar from the `effec_iono` broadcast parameters
        Args:
            effec_iono(list): Effective ionisation level coefficients a0 a1 a2, length 3
        Return:
            float: azpar
        """
        a0 = effec_iono[0]
        a1 = effec_iono[1]
        a2 = effec_iono[2]
        azpar = sqrt(a0 ** 2 + 1633.33 * (a1 ** 2) + 4802000 * (a2 ** 2) + 3266.67 * a0 * a2)
        return azpar

    @staticmethod
    def calculate_pierce_point_lat_lon(sat_el_rad: float, sat_az_rad: float, user_lat_rad: float,
                                       user_lon_rad: float):
        """Compute the geographic latitude and longitude of the ionospheric pierce point IPP
        """
        # Calculate the earths central angle delta_pp between the user position and the earth projection of the piercing
        # point
        lambda_pp_rad = PI / 2 - sat_el_rad - asin((Re_km / (Re_km + hI_km)) * cos(sat_el_rad))

        # Calculate pierce point latitude
        lat_pp_rad = asin(
            sin(user_lat_rad) * cos(lambda_pp_rad) + cos(user_lat_rad) * sin(lambda_pp_rad) * cos(sat_az_rad))

        # Calculate pierce point longitude
        lon_pp_rad = user_lon_rad + asin(sin(lambda_pp_rad) * sin(sat_az_rad) / cos(lat_pp_rad))

        return lat_pp_rad, lon_pp_rad

    @staticmethod
    def calculate_local_time(lon_pp_rad: float, universal_time_UT_hours: float):
        """Computes the local time LT given the `pierce point` longitude and UT
        """
        return universal_time_UT_hours + degrees(lon_pp_rad) / 15

    @staticmethod
    def calculate_sun_declination(doy: float):
        """Computes the sun's declination for the applicable day
        """
        return radians(23.44 * sin(radians(0.9856 * (doy - 80.7))))

    @staticmethod
    def calculate_solar_zenith_angle_dependence(lat_pp_rad: float, sun_declination_rad: float):
        """Computes the dependency of TEC with the solar zenith angle X
        """
        cos_solar_zenith_double_star = cos(lat_pp_rad - sun_declination_rad) + 0.4
        cos_solar_zenith_triple_star = cos(lat_pp_rad - sun_declination_rad) - (2 / PI) * lat_pp_rad * sin(
            sun_declination_rad)

        return cos_solar_zenith_double_star, cos_solar_zenith_triple_star

    @staticmethod
    def calculate_geomagnetic_lat_from_geographic(lat_pp_rad: float, lon_pp_rad: float):
        """Computes the geomagnetic latitude of the `pierce point` based on a magnetic dipole model
        """
        lat_north_pole_rad = radians(lat_north_pole_deg)
        lon_north_pole_rad = radians(lon_north_pole_deg)
        geo_lat_rad = asin(sin(lat_pp_rad) * sin(lat_north_pole_rad) + cos(lat_pp_rad) * cos(lat_north_pole_rad) *
                           cos(lon_pp_rad - lon_north_pole_rad))
        return geo_lat_rad

    @staticmethod
    def compute_local_time_dependency_F1(lat_pp_rad: float, lon_pp_rad: float, universal_time_hour: float):
        """Computes the local time dependency
        """
        cos_double_star, cos_triple_star = NTCMG.calculate_solar_zenith_angle_dependence(lat_pp_rad, lon_pp_rad)
        LT = NTCMG.calculate_local_time(lon_pp_rad, universal_time_hour)
        VD = 2 * PI * (LT - 14) / 24  # Diurnal angular phase
        VSD = 2 * PI * LT / 12  # Semi-diurnal
        VTD = 2 * PI * LT / 8  # ter-diurnal

        F1 = cos_triple_star + cos_double_star * (NTCMG_coefs.k1 * cos(VD) + NTCMG_coefs.k2 * cos(VSD) +
                                                  NTCMG_coefs.k3 * sin(VSD) + NTCMG_coefs.k4 * cos(VTD) +
                                                  NTCMG_coefs.k5 * sin(VTD))

        return F1

    @staticmethod
    def compute_season_dependency_F2(doy: int):
        """Computes the season dependency term
        """
        VA = 2 * PI * (doy - 18) / 365.25
        VSA = 4 * PI * (doy - 6) / 365.25
        F2 = 1 + NTCMG_coefs.k6 * cos(VA) + NTCMG_coefs.k7 * cos(VSA)

        return F2

    @staticmethod
    def compute_geomagnetic_field_dependency_F3(lat_pp_rad: float, lon_pp_rad: float):
        """Computes the geomagnetic field dependency term
        """
        geo_lat_pp_rad = NTCMG.calculate_geomagnetic_lat_from_geographic(lat_pp_rad, lon_pp_rad)
        F3 = 1 + NTCMG_coefs.k8 * cos(geo_lat_pp_rad)

        return F3

    @staticmethod
    def compute_equatorial_anomaly_dependency_F4(lat_pp_rad: float, lon_pp_rad: float):
        """Computes the equatorial anomaly dependency term
        """
        geo_lat_pp_deg = degrees(NTCMG.calculate_geomagnetic_lat_from_geographic(lat_pp_rad, lon_pp_rad))
        oc1_2 = 144
        oc2_2 = 169
        EC1 = - ((geo_lat_pp_deg - 16) ** 2) / (2 * oc1_2)
        EC2 = - ((geo_lat_pp_deg + 10) ** 2) / (2 * oc2_2)
        F4 = 1 + NTCMG_coefs.k9 * exp(EC1) + NTCMG_coefs.k10 * exp(EC2)

        return F4

    @staticmethod
    def compute_solar_activity_dependency_F5(azpar: float):
        """Computes the solar activity dependency term
        """
        F5 = NTCMG_coefs.k11 + NTCMG_coefs.k12 * azpar

        return F5

    @staticmethod
    def compute_vtec(lat_pp_rad: float, lon_pp_rad: float, ut1, azpar: float):
        """Computes the VTEC
        """
        doy = ut1.timetuple().tm_yday

        universal_time_hour = ut1.hour + ut1.minute / 60 + ut1.second / 3600  # Hour and decimal
        F1 = NTCMG.compute_local_time_dependency_F1(lat_pp_rad, lon_pp_rad, universal_time_hour)
        F2 = NTCMG.compute_season_dependency_F2(doy)
        F3 = NTCMG.compute_geomagnetic_field_dependency_F3(lat_pp_rad, lon_pp_rad)
        F4 = NTCMG.compute_equatorial_anomaly_dependency_F4(lat_pp_rad, lon_pp_rad)
        F5 = NTCMG.compute_solar_activity_dependency_F5(azpar)
        return F1 * F2 * F3 * F4 * F5

    @staticmethod
    def vtec_to_stec_mapping(sat_el_rad: float):
        """Computes the VTEC
        """
        sinz = (Re_km / (Re_km + hI_km)) * sin(0.9782 * ((PI / 2) - sat_el_rad))
        return 1 / sqrt(1 - sinz ** 2)

    @staticmethod
    def compute(ut1, gal_param, user_lat, user_long, sv_el, sv_az, freq):
        """
        Main function of the NTCM-G model. It computes the ionosphere correction using the a priori
        NTCM-G Model
        refs:
            * https://www.gsc-europa.eu/sites/default/files/NTCM-G_Ionospheric_Model_Description_-_v1.0.pdf

        Args:
            ut1(src.data_types.date.Epoch): required epoch in UT1 to compute the iono delay
            gal_param (list) : Effective Ionisation Level coefficients broadcast in the Galileo navigation message
            user_lat(float): user latitude in [rad]
            user_long(float): user longitude in [rad]
            sv_el(float): satellite elevation in [rad]
            sv_az(float): satellite azimuth in [rad]
            freq(src.data_types.gnss.data_type.DataType): Frequency band required for the iono delay
        Returns:
            float : ionosphere correction [m]
        """
        # Compute Azpar
        azpar = NTCMG.calculate_azpar(gal_param)

        # Calculate ionospheric pierce point location (lat_pp, lon_pp) for user-to-sat link at 450km height
        lat_pp, lon_pp = NTCMG.calculate_pierce_point_lat_lon(sv_el, sv_az, user_lat, user_long)

        # Call NTCM G to calculate VTEC at pierce point loc and local time LT
        vtec = NTCMG.compute_vtec(lat_pp, lon_pp, ut1, azpar)

        # Ionospheric mapping function MF
        # VTEC to STEC using the ionospheric mapping function
        stec = vtec * NTCMG.vtec_to_stec_mapping(sv_el)

        # There is still a need to multiply 40.3/f**2 by the stec to retrieve the correction
        return ((40.3 * 1e16) / (freq.freq_value ** 2)) * stec
