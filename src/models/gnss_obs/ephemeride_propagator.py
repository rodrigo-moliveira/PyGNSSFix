from math import sqrt, cos, sin

import numpy as np

from src import constants
from src.models.frames.frames import dcm_e_i, M2E, E2v


def fix_gnss_week_crossovers(time_diff: float) -> float:
    """
    Repairs over and underflow of GPS time, that is, the time difference must account for beginning or end of week
    crossovers.

    The time difference (time_diff) is the difference between a given GNSS epoch time t and toc:
        time_diff = t - toc
    time_diff is used, for example, in the computation of the clock bias, given the navigation clock model

    According to section [20.3.3.3.3.1] of **REF[3]**

    Args:
        time_diff (float) : time difference to be fixed
    Return:
        float : fixed time difference
    """
    half_week = 302400

    if time_diff > half_week:
        time_diff = time_diff - 2 * half_week

    elif time_diff < -half_week:
        time_diff = time_diff + 2 * half_week

    return time_diff


class EphemeridePropagator:

    @staticmethod
    def compute_sat_nav_position_dt_rel(nav_message, time_emission, transit) ->\
            tuple[np.array, float]:
        """
        Computes:
            * the satellite ephemeride
            * the clock relativistic correction
            * the true range rho
        at the requested epoch, and given the closest (valid) navigation data point and transit time
        the corresponding

        Args:
            nav_message (src.data_types.containers.NavigationData.NavigationPointGPS) : navigation data point object
            time_emission (src.data_types.basics.Epoch.Epoch) : Signal emission time (wrt GPS time system)
            transit (float) : Computed transit time in seconds. Used to rotate the computed satellite ephemeride to the
                            ECEF frame at reception time

        Returns:
            tuple [Position, ~numpy.array, float] : returns the computed satellite position at the request epoch, as
                                                    well as the corresponding clock relativistic corrections

        """
        # satellite coordinates in ECEF frame defined at TX time, relativistic correction for satellite clock
        r_sat, dt_relative = EphemeridePropagator.compute_nav_sat_pos(nav_message, time_emission)

        # rotation matrix from ECEF TX to ECEF RX (taking into consideration the signal transmission time)
        _R = dcm_e_i(-transit)

        # get satellite position vector at ECEF frame defined at RX time (to be compared with receiver position)
        p_sat = _R @ r_sat

        return p_sat, dt_relative

    @staticmethod
    def compute_nav_sat_pos(nav_message, epoch):
        """
        Implements the updating of GPS ephemerides (position) and the transformation to ECEF frame

        table 20-III [sec 20.3.3.4.3] of **REF[3]**
        """
        GM = constants.MU_WGS84 if nav_message.constellation == "GPS" else constants.MU_WGS84

        # unpack week and seconds of week
        _, sow = epoch

        # fetch navigation inputs
        M0 = getattr(nav_message, "M0")
        sqrtA = getattr(nav_message, "sqrtA")
        deltaN = getattr(nav_message, "deltaN")
        eccentricity = getattr(nav_message, "eccentricity")
        omega = getattr(nav_message, "omega")
        RAANDot = getattr(nav_message, "RAANDot")
        RAAN0 = getattr(nav_message, "RAAN0")
        cuc = getattr(nav_message, "cuc")
        cus = getattr(nav_message, "cus")
        crc = getattr(nav_message, "crc")
        crs = getattr(nav_message, "crs")
        i0 = getattr(nav_message, "i0")
        iDot = getattr(nav_message, "iDot")
        cic = getattr(nav_message, "cic")
        cis = getattr(nav_message, "cis")
        toe = getattr(nav_message, "toe")[1]  # to get seconds of week for TOE

        # semi major axis
        A = sqrtA * sqrtA
        A_3 = A * A * A

        # mean motion
        n = sqrt(GM / A_3)

        # time from ephemeris reference epoch (correct for beginning / end of week crossovers)
        dt = sow - toe
        dt = fix_gnss_week_crossovers(dt)

        # corrected mean motion
        n = n + deltaN

        # mean anomaly at epoch
        M = M0 + n * dt

        # eccentric anomaly (Kepler equation)
        E = M2E(eccentricity, M)

        # true anomaly
        v = E2v(eccentricity, E)

        # argument of latitude
        u = v + omega

        # corrections
        u_correction = cuc * cos(2 * u) + cus * sin(2 * u)
        radius_correction = crc * cos(2 * u) + crs * sin(2 * u)
        inclination_correction = cic * cos(2 * u) + cis * sin(2 * u)

        # apply corrections
        u = u + u_correction
        radius = A * (1 - eccentricity * cos(E)) + radius_correction
        i = i0 + inclination_correction + iDot * dt

        # SV position in orbital plane
        x_orbital = radius * cos(u)
        y_orbital = radius * sin(u)

        # corrected RAAN
        RAAN = RAAN0 + (RAANDot - constants.EARTH_ROTATION) * dt - constants.EARTH_ROTATION * toe

        # ECEF coordinates
        x_ECEF = x_orbital * cos(RAAN) - y_orbital * cos(i) * sin(RAAN)
        y_ECEF = x_orbital * sin(RAAN) + y_orbital * cos(i) * cos(RAAN)
        z_ECEF = y_orbital * sin(i)

        # TODO: add computation of sat velocities here

        # get StateVector object
        position = np.array([x_ECEF, y_ECEF, z_ECEF])

        # compute relativistic correction Eq 5.19 of **REF[1]**
        rel_correction = -2 * sqrt(GM) * sqrtA / constants.SPEED_OF_LIGHT ** 2 * eccentricity * sin(E)

        return position, rel_correction
