from math import sqrt, cos, sin

import numpy

from PositioningSolver.src.math_utils.Constants import Constant
from PositioningSolver.src.gnss.state_space.utils import M2E, E2v, matrix_ECEF2ECI
from PositioningSolver.src.gnss.state_space.gnss_state import PositionGNSS


def correct_gps_week_crossovers(time_diff: float) -> float:
    """
    Repairs over and underflow of GPS time, that is, the time difference
    must account for beginning or end of week crossovers.

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
    def get_sat_position_and_true_range(nav_message, time_emission, transit, rec_position, relativistic_correction) ->\
            tuple[PositionGNSS, numpy.array, float]:
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
            rec_position (Position) : the receiver position
            relativistic_correction (bool) : whether or not to compute the relativistic correction

        Returns:
            tuple [Position, ~numpy.array, float] : returns the computed satellite position at the request epoch, as
                                                    well as the corresponding clock relativistic corrections

        """
        # satellite coordinates in ECEF frame defined at TX time, relativistic correction for satellite clock
        r_sat, dt_relative = EphemeridePropagator.compute(nav_message, time_emission, relativistic_correction)

        # rotation matrix from ECEF TX to ECEF RX (taking into consideration the signal transmission time)
        _R = matrix_ECEF2ECI(-transit)

        # get satellite position vector at ECEF frame defined at RX time (to be compared with receiver position)
        p_sat = _R @ r_sat

        # compute true range
        rec_position.form = "cartesian"
        rho_0 = numpy.linalg.norm(p_sat - rec_position)

        return p_sat, rho_0, dt_relative

    @staticmethod
    def compute(nav_message, epoch, relativistic_correction, constellation="GPS") -> tuple[PositionGNSS, float]:
        """
        Computes the satellite ephemeride at the requested epoch, given the closest (valid) navigation data point

        Args:
            constellation (src.data_types.data_types.SatelliteSystem) : the constellation to process
            nav_message (src.data_types.containers.NavigationData.NavigationPointGPS) : navigation data point object
            epoch (src.data_types.basics.Epoch.Epoch) : Epoch to compute the ephemerides. For a correct implementation,
                                                        should be in GPS time (not SV nor receiver time)
            relativistic_correction (bool) : whether or not to compute the relativistic correction

        Returns:
            tuple [Position, float] : returns the computed satellite position at the request epoch, as well as the
                                      corresponding clock relativistic corrections

        """
        return EphemeridePropagator._compute_ephemeride_GPS(nav_message, epoch, relativistic_correction)
        # if constellation == "GPS":
        #    return EphemeridePropagator._compute_ephemeride_GPS(nav_message, epoch, relativistic_correction)
        # else:
        #    raise NotImplementedError(f"Galileo ephemeride propagator not yet implemented. Only GPS is currently "
        #                              f"possible.")

    @staticmethod
    def _compute_ephemeride_GPS(nav_message, epoch, relativistic_correction):
        """
        Implements the updating of GPS ephemerides (position) and the transformation to ECEF frame

        table 20-III [sec 20.3.3.4.3] of **REF[3]**
        """

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
        toe = getattr(nav_message, "toe")
        # toc = getattr(navPoint, "toc")
        # af0 = getattr(navPoint, "af0")
        # af1 = getattr(navPoint, "af1")
        # af2 = getattr(navPoint, "af2")
        # TGD = getattr(navPoint, "TGD")

        # semi major axis
        A = sqrtA * sqrtA
        A_3 = A * A * A

        # mean motion
        n = sqrt(Constant.MU / A_3)

        # time from ephemeris reference epoch (correct for beginning / end of week crossovers)
        dt = epoch - toe
        dt = correct_gps_week_crossovers(dt)

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
        RAAN = RAAN0 + (RAANDot - Constant.EARTH_ROTATION) * dt - Constant.EARTH_ROTATION * toe.seconds

        # ECEF coordinates
        x_ECEF = x_orbital * cos(RAAN) - y_orbital * cos(i) * sin(RAAN)
        y_ECEF = x_orbital * sin(RAAN) + y_orbital * cos(i) * cos(RAAN)
        z_ECEF = y_orbital * sin(i)

        # get StateVector object
        position = PositionGNSS([x_ECEF, y_ECEF, z_ECEF], "ECEF", "cartesian")

        # compute relativistic correction
        rel_correction = 0
        if relativistic_correction:
            # Eq 5.19 of **REF[1]**
            rel_correction = -2 * sqrt(Constant.MU) * sqrtA / Constant.SPEED_OF_LIGHT ** 2 * eccentricity * sin(E)
        return position, rel_correction
