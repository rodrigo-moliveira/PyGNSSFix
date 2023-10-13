import numpy as np
from dataclasses import dataclass

from math import sqrt, asin, cos, sin, degrees, radians, exp
from src.constants import PI
from src.models.frames import cartesian2geodetic

# Constants for this model (note that the Earth radius has this specific value for this model) 
Re_km = 6371
hI_km = 450
lat_north_pole_deg = 79.74
lon_north_pole_deg = -71.78


# NTCM-gG class to calculate the ionospheric correction for galileo users
# More information can be found in :
# https://www.gsc-europa.eu/sites/default/files/NTCM-G_Ionospheric_Model_Description_-_v1.0.pdf


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

    @staticmethod
    def calculate_azpar(effec_iono: np.array):
        """
        effec_iono: Effective ionisation level coefficients a0 a1 a2
        azpar
        """
        a0 = effec_iono[0]
        a1 = effec_iono[1]
        a2 = effec_iono[2]
        azpar = sqrt(a0 ** 2 + 1633.33 * (a1 ** 2) + 4802000 * (a2 ** 2) + 3266.67 * a0 * a2)
        return azpar

    @staticmethod
    def calculate_pierce_point_lat_lon(sat_el_rad: float, sat_az_rad: float, user_lat_rad: float,
                                       user_lon_rad: float):
        """
        Description: Compute the geographic latitude and longitude of the ionospheric pierce point IPP
        sat_el_rad: Satellite elevation in radians
        sat_az_rad: Satellite azimuth in radians
        user_lat_rad: user receiver latitude in radians
        user_lon_rad: user receiver longitude in radians
        return lat_pp_rad, lon_pp_rad
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
        """
        Description: Computes the local time LT given the pierce point longitude and UT
        @param lon_pp_rad: Pierce point longitude in radians
        @param universal_time_UT_hours: Hours and decimals
        @return local time LT in hours and decimals
        """
        return universal_time_UT_hours + degrees(lon_pp_rad) / 15

    @staticmethod
    def calculate_sun_declination(doy):
        """
        Description: Computes the sun's declination for the applicable day
        @param doy: Day of year
        @return sun's declination in radians
        """
        return radians(23.44 * sin(radians(0.9856 * (doy - 80.7))))

    @staticmethod
    def calculate_solar_zenith_angle_dependence(lat_pp_rad: float, sun_declination_rad: float):
        """
        Description: Computes the dependency of TEC with the solar zenith angle X
        @param lat_pp_rad: Pierce point latitude in radians
        @param sun_declination_rad: Sun's declination in radians
        @return two terms for the solar zenith angle dependence
        """
        cos_solar_zenith_double_star = cos(lat_pp_rad - sun_declination_rad) + 0.4
        cos_solar_zenith_triple_star = cos(lat_pp_rad - sun_declination_rad) - (2 / PI) * lat_pp_rad * sin(
            sun_declination_rad)

        return cos_solar_zenith_double_star, cos_solar_zenith_triple_star

    @staticmethod
    def calculate_geomagnetic_lat_from_geographic(lat_pp_rad: float, lon_pp_rad: float):
        """
        Description: Computes the geomagnetic latitude of the pierce based on a magnetic dipole model
        @param lat_pp_rad: Pierce point latitude in radians
        @param lon_pp_rad: Pierce point longitude in radians
        @return geomagnetic latitude in radians
        """
        lat_north_pole_rad = radians(lat_north_pole_deg)
        lon_north_pole_rad = radians(lon_north_pole_deg)
        geo_lat_rad = asin(sin(lat_pp_rad) * sin(lat_north_pole_rad) + cos(lat_pp_rad) * cos(lat_north_pole_rad) *
                           cos(lon_pp_rad - lon_north_pole_rad))
        return geo_lat_rad

    @staticmethod
    def compute_local_time_dependency_F1(lat_pp_rad: float, lon_pp_rad: float, universal_time_hour: float):
        """
        Description: Computes the local time dependency
        @param lat_pp_rad: Pierce point latitude in radians
        @param lon_pp_rad: Pierce point longitude in radians
        @param universal_time_hour: universal time, hours and decimals
        @return F1
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
        """
        Description: Computes the season dependency term
        @param doy: day of year
        @return F2
        """
        VA = 2 * PI * (doy - 18) / 365.25
        VSA = 4 * PI * (doy - 6) / 365.25
        F2 = 1 + NTCMG_coefs.k6 * cos(VA) + NTCMG_coefs.k7 * cos(VSA)

        return F2

    @staticmethod
    def compute_geomagnetic_field_dependency_F3(lat_pp_rad: float, lon_pp_rad: float):
        """
        Description: Computes the geomagnetic field dependency term
        @param lat_pp_rad: Pierce point latitude in radians
        @param lon_pp_rad: Pierce point longitude in radians
        @return F3
        """
        geo_lat_pp_rad = NTCMG.calculate_geomagnetic_lat_from_geographic(lat_pp_rad, lon_pp_rad)
        F3 = 1 + NTCMG_coefs.k8 * cos(geo_lat_pp_rad)

        return F3

    @staticmethod
    def compute_equatorial_anomaly_dependency_F4(lat_pp_rad: float, lon_pp_rad: float):
        """
        Description: Computes the equatorial anomaly dependency term
        @param lat_pp_rad: Pierce point latitude in radians
        @param lon_pp_rad: Pierce point longitude in radians
        @return F4
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
        """
        Description: Computes the solar activity dependency term
        @param azpar: proxy measure of the solar activity level
        @return F5
        """
        F5 = NTCMG_coefs.k11 + NTCMG_coefs.k12 * azpar

        return F5

    @staticmethod
    def compute_vtec(lat_pp_rad: float, lon_pp_rad: float, utc, azpar: float):
        """
        Description: Computes the VTEC
        @param lat_pp_rad: Pierce point latitude in radians
        @param lon_pp_rad: Pierce point longitude in radians
        @param utc: measure timestamp
        @param azpar: azpar
        @return VTEC
        """
        doy = utc.timetuple().tm_yday
        universal_time_hour = utc.hour + utc.minute / 60 + utc.second / 3600  # Hour and decimal
        F1 = NTCMG.compute_local_time_dependency_F1(lat_pp_rad, lon_pp_rad, universal_time_hour)
        F2 = NTCMG.compute_season_dependency_F2(doy)
        F3 = NTCMG.compute_geomagnetic_field_dependency_F3(lat_pp_rad, lon_pp_rad)
        F4 = NTCMG.compute_equatorial_anomaly_dependency_F4(lat_pp_rad, lon_pp_rad)
        F5 = NTCMG.compute_solar_activity_dependency_F5(azpar)
        return F1 * F2 * F3 * F4 * F5

    @staticmethod
    def vtec_to_stec_mapping(sat_el_rad: float):
        """
        Description: Computes the VTEC
        @param sat_el_rad: Satellite elevation in radians
        """
        sinz = (Re_km / (Re_km + hI_km)) * sin(0.9782 * ((PI / 2) - sat_el_rad))
        return 1 / sqrt(1 - sinz ** 2)

    @staticmethod
    def calculate_ionospheric_contribution(ut1, effec_iono: np.array, rec_pos: np.array,
                                           el_rad: float, az_rad: float, frequency: float):
        """
        main function
        Description: Computes the VTEC
        @param ut1: Timestamp in UT1
        @param effec_iono: Effective ionisation level coefficients a0 a1 a2
        @param rec_pos: Receiver position
        @param el_rad: Satellite elevation
        @param az_rad: Satellite azimuth
        @param frequency: Measure frequency
        @return Ionosphere delay
        """
        # INPUT OR CALL calculate_azpar to get azpar
        azpar = NTCMG.calculate_azpar(effec_iono)

        # [X] Calculate satellite elevation E and azimuth angles
        # Calculate ionospheric pierce point location (lat_pp, lon_pp) for user-to-sat link at 450km height
        llh = cartesian2geodetic(tuple(rec_pos))
        user_lat_rad = llh[0]
        user_lon_rad = llh[1]
        lat_pp_rad, lon_pp_rad = NTCMG.calculate_pierce_point_lat_lon(el_rad, az_rad, user_lat_rad, user_lon_rad)

        # Call NTCM G to calculate VTEC at pierce point loc and local time LT
        vtec = NTCMG.compute_vtec(lat_pp_rad, lon_pp_rad, ut1, azpar)

        # Ionospheric mapping function MF
        # VTEC to STEC using the ionospheric mapping function
        stec = vtec * NTCMG.vtec_to_stec_mapping(el_rad)

        # There is still a need to multiply 40.3/f**2 by the stec to retrieve the correction
        return ((40.3 * 1e16) / (frequency ** 2)) * stec
