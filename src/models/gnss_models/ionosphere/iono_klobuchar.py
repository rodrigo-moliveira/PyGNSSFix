"""Klobuchar ionospheric module for GPS
"""
from math import sin, cos

from src import constants
from src.data_types.gnss.data_type import L1


class IonoKlobuchar:
    """
    Implements the Ionospheric Klobuchar model for GPS satellites
    """

    @staticmethod
    def compute(gps_sow, alfa, beta, user_lat, user_long, sv_el, sv_az, freq):
        """
        This function computes the ionosphere correction using the a priori Klobuchar Ionospheric Model
        refs:
            * section 20.3.3.5.2.5 of [NAVSTAR GPS Space Segment/Navigation User Interfaces]
            * https://gssc.esa.int/navipedia/index.php/Klobuchar_Ionospheric_Model

        GPS satellites broadcast the parameters of the Klobuchar ionospheric model for single frequency users.
        The Klobuchar model was designed to minimise user computational complexity and user computer storage as far
        as to keep a minimum number of coefficients to transmit on satellite-user link.
        The model is estimated to reduce about the 50% RMS ionospheric range error worldwide

        Args:
            gps_sow (float) : seconds of GPS week for the epoch at evaluation
            alfa (list) : list of alfa parameters, length 4
            beta (list) : list of beta parameters, length 4
            user_lat(float): user latitude in [rad]
            user_long(float): user longitude in [rad]
            sv_el(float): satellite elevation in [rad]
            sv_az(float): satellite azimuth in [rad]
            freq(src.data_types.gnss.data_type.DataType): Frequency band required for the iono delay
        Returns:
            float : ionosphere correction [m]
        """

        # Get angles in semicircles
        sv_el_semi = sv_el / constants.PI
        user_lat_semi = user_lat / constants.PI
        user_long_semi = user_long / constants.PI

        # Calculate the earth-centred angle (elevation in semicircles)
        psi = 0.0137 / (sv_el_semi + 0.11) - 0.022

        # Compute and fix the latitude of the Ionospheric Pierce Point(IPP)
        lat_IPP = user_lat_semi + psi * cos(sv_az)
        if lat_IPP > 0.416:
            lat_IPP = 0.416
        elif lat_IPP < -0.416:
            lat_IPP = -0.416

        # Compute the longitude of the IPP
        long_IPP = user_long_semi + (psi * sin(sv_az)) / cos(lat_IPP * constants.PI)

        # Find the geomagnetic latitude of the IPP
        lat_m = lat_IPP + 0.064 * cos((long_IPP - 1.617) * constants.PI)

        # Find the local time at the IPP
        t = constants.SECONDS_IN_DAY / 2 * long_IPP + gps_sow
        t = t % constants.SECONDS_IN_DAY

        # Compute the amplitude of ionospheric delay.
        A_I = alfa[0] + alfa[1] * lat_m + alfa[2] * (lat_m ** 2) + alfa[3] * (lat_m ** 3)
        if A_I < 0:
            A_I = 0

        # Compute the period of ionospheric delay
        P_I = beta[0] + beta[1] * lat_m + beta[2] * (lat_m ** 2) + beta[3] * (lat_m ** 3)
        if P_I < 72000:
            P_I = 72000

        # Compute the phase of ionospheric delay
        X_I = 2 * constants.PI * (t - 50400) / P_I

        # Compute the slant factor (elevation in semicircles).
        F = 1.0 + 16.0 * (0.53 - sv_el_semi) ** 3

        # Compute the ionospheric time delay for L1 frequency/band [s]
        if abs(X_I) > 1.57:
            iono = 5E-9 * F
        else:
            iono = (5E-9 + A_I * (1 - (X_I ** 2) / 2 + (X_I ** 4) / 24)) * F

        # get ionosphere in meters for L1 frequency/band
        iono = iono * constants.SPEED_OF_LIGHT

        # convert from L1 to the user frequency
        iono = (L1.freq_value / freq.freq_value) ** 2 * iono

        return iono
