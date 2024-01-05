# Constants for the tropospheric Saastamoinen model
import numpy as np
from math import cos

from src import constants


class SaastamoinenTropo:
    LimLat = (15, 30, 45, 60, 75)

    #                   P0(mbar),T0(K),e0(mbar),beta(K/m),lambda0
    P_mean = np.array([[1013.25, 299.65, 26.31, 6.30e-3, 2.77],
                       [1017.25, 294.15, 21.79, 6.05e-3, 3.15],
                       [1015.75, 283.15, 11.66, 5.58e-3, 2.57],
                       [1011.75, 272.15, 6.78, 5.39e-3, 1.81],
                       [1013.00, 263.65, 4.11, 4.53e-3, 1.55]])

    #                     DP(mbar),DT(K),De(mbar),Db(K/m),dl
    P_season = np.array([[0.00, 0.00, 0.00, 0.00e-3, 0.00],
                         [-3.75, 7.00, 8.85, 0.25e-3, 0.33],
                         [-2.25, 11.00, 7.24, 0.32e-3, 0.46],
                         [-1.75, 15.00, 5.36, 0.81e-3, 0.74],
                         [-0.50, 14.50, 3.39, 0.62e-3, 0.30]])

    @staticmethod
    def compute(h, lat, doy):
        """
        This function computes the tropospheric correction of pseudorange measurements for online users,
        using the a priori Saastamoinen Model

        This code has been adapted from:
            Open source PANG-NAV software (https://geodesy.noaa.gov/gps-toolbox/PANG-NAV.htm), tropo_correction0.m script,
            implementing the Saastamoinen Model


        Args:
            h (float) : user altitude (geodetic coordinate)     [m]
            lat (float) : user latitude (geodetic coordinate)   [rad]
            doy (float) : Day of the year (from 1 to 365)       [1 - 365]

        Return:
            float : tropospheric correction [m]
        """

        # convert lat to degrees
        lat = constants.RAD2DEG * lat

        if lat < 0:
            D_star = 211
        else:
            D_star = 28

        if lat > SaastamoinenTropo.LimLat[-1]:  # if lat > 75 [deg]
            P0 = SaastamoinenTropo.P_mean[-1]
            DP = SaastamoinenTropo.P_season[-1]

        elif lat < SaastamoinenTropo.LimLat[0]:  # if lat < 15 [deg]
            P0 = SaastamoinenTropo.P_mean[0]
            DP = SaastamoinenTropo.P_season[0]

        else:
            _iis = [lat_i < lat for lat_i in SaastamoinenTropo.LimLat]
            _iis.reverse()
            _i = 4 - _iis.index(True)
            m = (lat - SaastamoinenTropo.LimLat[_i]) / (SaastamoinenTropo.LimLat[_i + 1] - SaastamoinenTropo.LimLat[_i])
            P0 = SaastamoinenTropo.P_mean[_i] + (SaastamoinenTropo.P_mean[_i + 1] - SaastamoinenTropo.P_mean[_i]) * m
            DP = SaastamoinenTropo.P_season[_i] + (SaastamoinenTropo.P_season[_i + 1] -
                                                   SaastamoinenTropo.P_season[_i]) * m

        Par = P0 - DP * cos(2 * np.pi * (doy - D_star) / 365.25)

        P, T, e, beta, Lambda = Par[:]

        D_z_dry, D_z_wet = SaastamoinenTropo.saasthyd(P, lat, e, T, h)

        return D_z_dry, D_z_wet

    @staticmethod
    def saasthyd(p, lat, e, T, hell):
        """
        This subroutine determines the zenith hydrostatic delay based on the
        equation by Saastamoinen (1972) as refined by Davis et al. (1985)

        Reference:
        Saastamoinen, J., Atmospheric correction for the troposphere and
        stratosphere in radio ranging of satellites. The use of artificial
        satellites for geodesy, Geophys. Monogr. Ser. 15, Amer. Geophys. Union,
        pp. 274-251, 1972.
        Davis, J.L, T.A. Herring, I.I. Shapiro, A.E.E. Rogers, and G. Elgered,
        Geodesy by Radio Interferometry: Effects of Atmospheric Modeling Errors
        on Estimates of Baseline Length, Radio Science, Vol. 20, No. 6,
        pp. 1593-1607, 1985.

        input parameters:
        p:     pressure in hPa
        lat:  ellipsoidal latitude in radians
        hell:  ellipsoidal height in m

        output parameters:
        zhd:  zenith hydrostatic delay in meter
        """
        # calculate denominator f
        f = 1 - 0.00266 * cos(2 * lat) - 0.00000028 * hell

        # calculate the zenith hydrostatic delay
        D_z_dry = 0.0022768*p/f

        # calculate the zenith wet delay
        D_z_wet = (0.002277 * (1255 / T + 0.05) * e) / f

        return D_z_dry, D_z_wet
