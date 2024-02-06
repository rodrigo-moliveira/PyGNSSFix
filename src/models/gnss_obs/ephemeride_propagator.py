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
    def compute_sat_nav_position_dt_rel(nav_message, time_emission, transit) -> \
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
        r_sat, v_sat, dt_relative = EphemeridePropagator.compute_nav_sat_pos(nav_message, time_emission)

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
        GM = constants.MU_WGS84 if nav_message.constellation == "GPS" else constants.MU_GTRF

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
        E_dot = n / (1.0 - eccentricity * cos(E))

        # true anomaly
        v = E2v(eccentricity, E)
        v_dot = sin(E) * E_dot * (1.0 + eccentricity * cos(v)) / (sin(v) * (1.0 - eccentricity * cos(E)))

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
        u_k_dot = v_dot + 2.0 * (cus * cos(2.0 * u) - cuc * sin(2.0 * u)) * v_dot
        r_k_dot = A * eccentricity * sin(E) * n / (1.0 - eccentricity * cos(E)) + 2.0 * (crs * cos(2.0 * u) - crc *
                                                                                         sin(2.0 * u)) * v_dot
        i_k_dot = iDot + (cis * cos(2.0 * u) - cic * sin(2.0 * u)) * 2.0 * v_dot

        # SV position in orbital plane
        x_orbital = radius * cos(u)
        y_orbital = radius * sin(u)
        x_kp_dot = r_k_dot * cos(u) - y_orbital * u_k_dot
        y_kp_dot = r_k_dot * sin(u) + x_orbital * u_k_dot

        # corrected longitude of ascending node
        OMEGA_k_dot = (RAANDot - constants.EARTH_ROTATION)
        OMEGA_k = RAAN0 + OMEGA_k_dot * dt - constants.EARTH_ROTATION * toe

        # ECEF coordinates
        x_ECEF = x_orbital * cos(OMEGA_k) - y_orbital * cos(i) * sin(OMEGA_k)
        y_ECEF = x_orbital * sin(OMEGA_k) + y_orbital * cos(i) * cos(OMEGA_k)
        z_ECEF = y_orbital * sin(i)

        vx_ECEF = (x_kp_dot - y_orbital * cos(i) * OMEGA_k_dot) * cos(OMEGA_k) - \
                  (x_orbital * OMEGA_k_dot + y_kp_dot * cos(i) - y_orbital * sin(i) * i_k_dot) * sin(OMEGA_k)

        vy_ECEF = (x_kp_dot - y_orbital * cos(i) * OMEGA_k_dot) * sin(OMEGA_k) + \
                  (x_orbital * OMEGA_k_dot + y_kp_dot * cos(i) - y_orbital * sin(i) * i_k_dot) * cos(OMEGA_k)

        vz_ECEF = y_kp_dot * sin(i) + y_orbital * cos(i) * i_k_dot

        # get StateVector object
        position = np.array([x_ECEF, y_ECEF, z_ECEF])
        velocity = np.array([vx_ECEF, vy_ECEF, vz_ECEF])

        # compute relativistic correction Eq 5.19 of **REF[1]**
        rel_correction = -2 * sqrt(GM) * sqrtA / constants.SPEED_OF_LIGHT ** 2 * eccentricity * sin(E)

        return position, velocity, rel_correction
