""" Module with several utility functions for performing navigation (computation of satellite clock,
orbit ephemerides, and GNSS time)
"""

from datetime import timedelta
import numpy as np
from numpy import sqrt, cos, sin, linalg, array

from src import constants
from src.data_types.gnss import get_data_type
from src.io.config import EnumTransmissionTime
from src.common_log import get_logger, MODEL_LOG
from src.models.frames import dcm_e_i, M2E, E2v


_ggto_cache = {}  # local cache of GGTO computations
L1 = get_data_type("L1", "GPS")
L2 = get_data_type("L2", "GPS")

E1 = get_data_type("E1", "GAL")
E5a = get_data_type("E5a", "GAL")
E5b = get_data_type("E5b", "GAL")

# static control variables for issuing warnings
_warning_cache = {
    "_ggto_warning": False,
    "_gps_tgd_warning": False,
    "_gal_bgd_fnav_warning": False,
    "_gal_bgd_inav_warning": False
}


def fix_gnss_week_crossovers(time_diff: float) -> float:
    """
    Repairs over and underflow of GNSS time, that is, the time difference must account for beginning or end of week
    crossovers.

    The time difference (time_diff) is the difference between a given GNSS epoch time t and toc:
        time_diff = t - toc
    time_diff is used, for example, in the computation of the clock bias, given the navigation clock model

    Reference:
        [1] NAVSTAR GPS Space Segment/Navigation User Interfaces (IS-GPS-200). May 2021. Section 20.3.3.3.3.1

    Args:
        time_diff (float) : time difference to be fixed
    Returns:
        float : fixed time difference
    """
    half_week = 302400

    if time_diff > half_week:
        time_diff = time_diff - 2 * half_week

    elif time_diff < -half_week:
        time_diff = time_diff + 2 * half_week

    return time_diff


def compute_tx_time(model=None, **kwargs):
    """
    Main function to compute the transmission time of the satellite-receiver GNSS link.
    There are two alternatives, depending on the user configuration:
        * Geometric algorithm: Purely Geometric Algorithm from section 5.1.1.1 of [1]
        * Pseudorange-based algorithm: A Pseudorange-Based Algorithm from section 5.1.1.1 of [1]
        * NAPEOS algorithm: the transmission time and true range are computed using satellite position and velocity
            data evaluated at the reception time (see [2])

    Reference:
        [1] ESA GNSS DATA PROCESSING, Volume I: Fundamentals and Algorithms, J. Sanz Subirana,
        J.M. Juan Zornoza and M. Hernández-Pajares
        [2] NAPEOS Mathematical Models and Algorithms, DOPS-SYS-TN-0100-OPS-GN, 5-NOV-2009

    Args:
        model(EnumTransmissionTime): enumeration that evaluates to `EnumTransmissionTime.GEOMETRIC` or
            `EnumTransmissionTime.PSEUDORANGE`
        kwargs(dict): depending on the model, different items are provided in `kwargs`:
            * Geometric algorithm: `r_receiver`, `t_reception`, `dt_receiver`, `sat`, `sat_orbits`.
            * Pseudorange-based algorithm: `pseudorange_obs`, `t_reception`, `sat`, `sat_clocks`.
    Returns:
        tuple [src.data_types.date.Epoch, float] : computed TX epoch, computed transit time
    """
    if model == EnumTransmissionTime.GEOMETRIC:
        return tx_time_geometric(**kwargs)
    elif model == EnumTransmissionTime.PSEUDORANGE:
        return tx_time_pseudorange(**kwargs)
    elif model == EnumTransmissionTime.NAPEOS:
        return tx_time_napeos(**kwargs)


def tx_time_geometric(r_receiver=None, t_reception=None, dt_receiver=None, sat=None, sat_orbits=None, **_):
    """
    Purely Geometric Algorithm to compute:
        * transit time (propagation time from signal transmission to reception)
        * Transmission time in GNSS time system

    Fundamental equation
        t(emission)^{receiver} = t(reception)^{receiver} - tau

    Procedure:
        1. `tau = 0`
        2. get satellite position `r_satellite` at `t = t(emission)^{receiver} = t(reception)^{receiver} - tau`
        3. get pseudorange `rho = norm(R(-tau) @ r_satellite - r_receiver)`
        4. `tau = rho / c`
        5. go to step 2 and iterate until convergence
    where `R()` is the rotation matrix from ECEF to ECI, and `c` is the speed of light

    Args:
        r_receiver (numpy.ndarray): position of receiver at reception time `t(reception)^{receiver}`
        t_reception (src.data_types.date.Epoch): time of signal reception, measured by receiver
            `t(reception)^{receiver}`, i.e., RINEX obs time tag
        sat(src.data_types.gnss.Satellite): the satellite to compute the transmission time
        sat_orbits(src.data_mng.gnss.sat_orbit_data.SatelliteOrbits): `SatelliteOrbits` object with orbit data
        dt_receiver (float): current estimate of the receiver clock bias (seconds)
    Returns:
        tuple [src.data_types.date.Epoch, float] : computed TX epoch, computed transit time
    """
    max_iter = 10
    residual_th = 1E-8
    rho_previous = 0
    residual = 1
    N = 0

    # 1. Initialize tau
    tau = 0

    while residual > residual_th and N < max_iter:
        # 2. Get satellite coordinates
        t = t_reception + timedelta(seconds=-tau)
        r_satellite, _, _, _ = sat_orbits.get_orbit(sat, t)

        # 3. Compute pseudorange (in ECEF frame associated to t_receiver epoch)
        _R = dcm_e_i(-tau)
        rho = linalg.norm(_R @ r_satellite - r_receiver)
        # LOS = (_R @ r_satellite - r_receiver) / rho

        # 4. Compute new tau [s]
        tau = rho / constants.SPEED_OF_LIGHT

        # 5. Check convergence
        residual = abs(rho - rho_previous)
        rho_previous = rho
        N += 1

    # compute transmission time, using Eq. (5.9)
    T_emission = t_reception + timedelta(seconds=-tau - dt_receiver)
    return T_emission, tau


def tx_time_pseudorange(pseudorange_obs=None, t_reception=None, sat=None, sat_clocks=None, **_):
    """
    Pseudorange-Based Algorithm to compute:
        * transit time (propagation time from signal transmission to reception)
        * Transmission time in GNSS time system

    Fundamental equation
        `pseudorange = c (t(reception)^{receiver} - t(emission)^{satellite})`
            -> tau = t(reception)^{receiver} - t(emission)^{satellite}

    Args:
        pseudorange_obs (src.data_types.gnss.observation.Observation): instance of `Observation` with the
            pseudorange measured observable
        t_reception (src.data_types.date.Epoch): time of signal reception, measured by receiver
            `t(reception)^{receiver}`, i.e., RINEX obs time tag
        sat(src.data_types.gnss.Satellite): the satellite to compute the transmission time
        sat_clocks(src.data_mng.gnss.sat_clock_data.SatelliteClocks): `SatelliteClocks` object with clock data
    Returns:
        tuple [src.data_types.date.Epoch, float] : computed TX epoch, computed transit time
    """
    tau = pseudorange_obs.value / constants.SPEED_OF_LIGHT
    t_emission = t_reception + timedelta(seconds=-tau)  # t(emission)^{satellite} = t(reception)^{receiver} - tau

    dt_sat_, _ = sat_clocks.get_clock(sat, t_emission)

    # correct satellite clock for BGDs
    nav_message = sat_clocks.nav_data.get_closest_message(sat, t_emission)
    dt_sat = nav_sat_clock_correction(dt_sat_, pseudorange_obs.datatype, nav_message)

    # compute transmission time, using Eq. (5.6)
    T_emission = t_emission + timedelta(seconds=-dt_sat)  # t(emission)^{GPS} = t(emission)^{satellite} - dt_sat
    return T_emission, tau


def tx_time_napeos(r_receiver=None, v_receiver=None, dt_receiver=None, sat=None, t_reception=None,
                   sat_orbits=None, **_):
    """
    NAPEOS algorithm: the transmission time and true range are computed using only satellite position and velocity
            data evaluated at the reception time (see [1])

    Reference:
        [1] NAPEOS Mathematical Models and Algorithms, DOPS-SYS-TN-0100-OPS-GN, 5-NOV-2009

    Args:
        r_receiver (numpy.ndarray): position of receiver at reception time `t(reception)^{receiver}`
        v_receiver (numpy.ndarray): velocity of receiver at reception time `t(reception)^{receiver}`
        dt_receiver (float): current estimate of the receiver clock bias (seconds)
        t_reception (src.data_types.date.Epoch): time of signal reception, measured by receiver
            `t(reception)^{receiver}`, i.e., RINEX obs time tag
        sat(src.data_types.gnss.Satellite): the satellite to compute the transmission time
        sat_orbits(src.data_mng.gnss.sat_orbit_data.SatelliteOrbits): `SatelliteOrbits` object with orbit data
    Returns:
        tuple [src.data_types.date.Epoch, float] : computed TX epoch, computed transit time
    """
    MU = constants.MU_WGS84 if sat.sat_system.is_gps() else constants.MU_GTRF

    # get satellite position and velocity at TX ECEF frame
    pos_sat, vel_sat, _, _ = sat_orbits.get_orbit(sat, t_reception)
    vel_sat += np.cross(constants.EARTH_ANGULAR_RATE, pos_sat)  # convert to inertial
    pos_rec = r_receiver
    vel_rec = v_receiver + np.cross(constants.EARTH_ANGULAR_RATE, pos_rec)  # convert to inertial

    sat_radius = np.linalg.norm(pos_sat)
    pos_dif = pos_sat - pos_rec
    vel_dif = vel_sat - vel_rec

    geom_range = np.linalg.norm(pos_dif)
    unit_los = pos_dif / geom_range
    range_rate = np.dot(unit_los, vel_dif)

    # corrections
    nGravitationalEffect = np.dot(pos_sat, unit_los) * MU / (sat_radius ** 3)
    nSecondOrderRangeRateEffect = -(range_rate * range_rate) / geom_range
    nSecondOrderVelocityEffect = np.dot(vel_dif, vel_dif) / geom_range
    nTotalSecondOrderEffect = nSecondOrderRangeRateEffect + nSecondOrderVelocityEffect + nGravitationalEffect
    nTimeTagError = -dt_receiver
    nTimeTagCorrection = nTimeTagError * (range_rate + nTimeTagError * nTotalSecondOrderEffect / 2.0)
    nZeroOrderFlightTime = geom_range / constants.SPEED_OF_LIGHT
    nVelocityCorrection = -nZeroOrderFlightTime * np.dot(vel_sat, unit_los)
    nAccelCorrection = -(nZeroOrderFlightTime ** 2) * nGravitationalEffect / 2.0
    nTimeCorrection = nVelocityCorrection + nAccelCorrection
    # nRelCorrection = 2 * np.dot(pos_sat, vel_sat) / constants.SPEED_OF_LIGHT
    # true_range = geom_range - nTimeTagCorrection + nTimeCorrection  # + nRelCorrection
    tau = (geom_range - nTimeTagCorrection + nTimeCorrection) / constants.SPEED_OF_LIGHT

    T_emission = t_reception + timedelta(seconds=-tau - dt_receiver)
    return T_emission, tau


def broadcast_clock(af0: float, af1: float, af2: float, toc: float, tx_raw: float) -> tuple[float, float]:
    """
    Algorithm to compute the broadcast satellite clock bias and clock drift (GPS and GAL), using
    the provided broadcast navigation clock parameters.

    The polynomial defined in the following allows the user to determine the effective SV PRN code phase offset
    referenced to the phase center of the antennas (delta t sv) with respect to GNSS system time (t) at the time of
    inputs transmission.

    This estimated correction accounts for the deterministic SV clock error characteristics of bias, drift and
    aging, as well as for the SV implementation characteristics of group delay bias and mean differential group delay.
    -> These coefficients do not include corrections for relativistic effects.

    This algorithm is detailed in [1]

    Reference:
        [1] NAVSTAR GPS Space Segment/Navigation User Interfaces (IS-GPS-200). May 2021. Section 20.3.3.3.3.1


    The implemented equations are:
        Δt(SV) = af0 + af1(t - toc) + af2(t - toc)^2
        d(Δt(SV))/dt = af1 + 2 * af2(t - toc)

    Args:
        af0 (float) : clock model constant parameter A0 [sec]
        af1 (float) : clock model first order parameter A1 [sec/sec]
        af2 (float) : clock model second order parameter A2 [sec/sec^2]
        toc (float) : time of clock (in seconds of week). This is the reference epoch for the clock model polynomial
        tx_raw (float) : epoch (in seconds of week) to compute the bias and drift (time of signal transmission)
    Returns:
        tuple [float, float] : clock bias[s], clock drift[ ]
    """
    dt = fix_gnss_week_crossovers(tx_raw - toc)
    clock_bias = af0 + af1 * dt + af2 * dt * dt
    clock_drift = af1 + 2 * af2 * dt

    return clock_bias, clock_drift


def nav_sat_clock_correction(sat_clock, datatype, nav_message):
    """
    This function corrects the satellite clock using broadcast BGDs (in case of GAL) or TGD (in case of GPS)

    The broadcast satellite clocks and the ones provided in RINEX Clock files (precise clocks) correspond to the
    ionosphere-free combination of dual-frequency code observations, that is, these clocks contain the code bias
    for the iono-free combination.

    This function uses the BGDs/TGDs to convert the input satellite clock for the provided datatype:
        * GPS - in the case of GPS, the broadcast clocks correspond to the ionosphere-free combination of L1(P)+L2
        * GAL - in the case of GAL, the broadcast clocks correspond to the ionosphere-free combination of:
            1. E1+E5a in case of F/NAV messages
            2. E1+E5b in case of I/NAV messages

    Note:
        In case of GPS users employing the L1 C/A code instead of the P1 code, another correction should be taken
        into account, which accounts for the difference in hardware biases between the P1 and C/A.
        This correction is however not transmitted in the GPS legacy navigation message.
        It must be fetched from IGS DCB products
        This is **currently** not considered in this algorithm!

    Args:
        sat_clock(float): the iono-free satellite clock bias
        datatype(src.data_types.gnss.data_type.DataType)
        nav_message (src.data_mng.gnss.navigation_data.NavigationPoint): instance of `NavigationPoint` used
            to fetch to BGDs/TGD
    Returns:
        float: the corrected satellite clock for the provided datatype
    """
    bgd = get_bgd_correction(datatype, nav_message)
    return sat_clock - bgd


def get_bgd_correction(datatype, nav_message):
    """
    This function fetches the appropriate BGD/TGD correction for the provided datatype

    Note:
        This function is only called assuming standard SPS Mode, that is, the only BGDs available
        are the ones from the navigation message. No extra DCB information is used

    Args:
        datatype(src.data_types.gnss.data_type.DataType)
        nav_message (src.data_mng.gnss.navigation_data.NavigationPoint): instance of `NavigationPoint` used
            to fetch the BGDs/TGDs
    Returns:
        float: the appropriate TGD/BGD correction for the provided datatype
    """
    #
    if nav_message.constellation == "GPS" and datatype.constellation == "GPS":
        if datatype.freq_number == 1 or datatype.freq_number == 2:
            scale = (L1.freq_value / datatype.freq.freq_value) ** 2
            return scale * nav_message.TGD
        elif datatype.freq_number == 12:
            # the broadcast clock refers to the L1-L2 iono-free clock -> no correction needed
            return 0.0
        else:
            if not _warning_cache["_gps_tgd_warning"]:
                log = get_logger(MODEL_LOG)
                log.warning(f"Unable to get GPS TGD for datatype {datatype}. Only L1/L2 corrections are available")
                _warning_cache["_gps_tgd_warning"] = True
            return 0.0
    elif nav_message.constellation == "GAL" and datatype.constellation == "GAL":
        nav_type = nav_message.nav_type

        if nav_type == "FNAV":
            if datatype.freq_number == 1 or datatype.freq_number == 5:
                scale = (E1.freq_value / datatype.freq.freq_value) ** 2
                return scale * nav_message.BGDE1E5a
            elif datatype.freq_number == 15:
                return 0.0
            else:
                if not _warning_cache["_gal_bgd_fnav_warning"]:
                    log = get_logger(MODEL_LOG)
                    log.warning(f"Unable to get GAL BGD for datatype {datatype}. The provided FNAV message only "
                                f"contains the BGD for E1/E5a corrections")
                    _warning_cache["_gal_bgd_fnav_warning"] = True
                return 0.0
        elif nav_type == "INAV":
            if datatype.freq_number == 1 or datatype.freq_number == 7:
                scale = (E1.freq_value / datatype.freq.freq_value) ** 2
                return scale * nav_message.BGDE1E5b
            elif datatype.freq_number == 17:
                return 0.0
            elif datatype.freq_number == 5:
                scale = (E1.freq_value / datatype.freq.freq_value) ** 2
                return nav_message.BGDE1E5b - nav_message.BGDE1E5a + nav_message.BGDE1E5a * scale
            elif datatype.freq_number == 15:
                return nav_message.BGDE1E5b - nav_message.BGDE1E5a
            else:
                if not _warning_cache["_gal_bgd_inav_warning"]:
                    log = get_logger(MODEL_LOG)
                    log.warning(f"Unable to get GAL BGD for datatype {datatype}. The provided INAV message only "
                                f"contains the BGD for E1/E5a and E1/E5b corrections")
                    _warning_cache["_gal_bgd_inav_warning"] = True
                return 0.0

    else:
        raise AttributeError(f"Mismatch between navigation message constellation ({nav_message.constellation}) and "
                             f" datatype {repr(datatype)}")
    return 0.0


def compute_ggto(time_correction, week, sow):
    """
    Function to compute the GGTO (GPS-to-Galileo Time Offset)

    (GGTO) is transmitted in the broadcast navigation message and available
    in the RINEX NAV file header (GPGA entry of header line TIME SYSTEM CORR).

    It contains:
        * a0 and a1 coefficients of the polynomial
        * reference time T_ref (seconds into GNSS week)
        * reference Week number Week_ref

    So the computation of GGTO is the following:
        GGTO(t_sow) = a0 + a1 * (t_sow - T + 604800 * (week - Week_ref) )

    where (week, t_sow) is the current time in week number and seconds of week

    Args:
        time_correction(dict or None): dictionary with the `GGTO` item.
            See :py:class:`src.data_mng.gnss.navigation_data.NavigationData`
        week(int): the week number of the user epoch to compute the GGTO
        sow(float): the SoW (seconds of week) of the user epoch to compute the GGTO
    Returns:
        float: GGTO offset for the provided epoch, in [seconds]
    """
    epoch = (week, sow)

    if time_correction is None or "GGTO" not in time_correction:
        if not _warning_cache["_ggto_warning"]:
            log = get_logger(MODEL_LOG)
            log.warning(f"The provided `time_correction` argument does not contain the GGTO correction. Please check "
                        f"the RINEX Navigation file header for the GPGA entry of header line TIME SYSTEM CORR. "
                        f"Setting GGTO to 0...")
            _warning_cache["_ggto_warning"] = True
        _ggto_cache[epoch] = 0.0

    if epoch not in _ggto_cache:
        # fetch data
        a0 = time_correction["GGTO"][0]
        a1 = time_correction["GGTO"][1]
        sow_ref = time_correction["GGTO"][2]
        week_ref = time_correction["GGTO"][3]

        # apply this eq. GGTO(t_sow) = a0 + a1 * (t_sow - T + 604800 * (week - Week_ref) )
        _ggto_cache[epoch] = a0 + a1 * (sow - sow_ref + constants.SECONDS_IN_GNSS_WEEK * (week - week_ref))
    return _ggto_cache[epoch]


def compute_nav_sat_eph(nav_message, epoch):
    """
    Compute the SV position and velocity from the broadcast ephemerides (navigation message).
    Also compute the associated relativistic time correction and drift
    Implements the updating of GPS ephemerides (position) and the transformation to ECEF frame

    Reference:
        [1] NAVSTAR GPS Space Segment/Navigation User Interfaces (IS-GPS-200). May 2021. Section 20.3.3.4.3

        [2] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
            Springer Cham, 2017

    Args:
        nav_message (src.data_mng.gnss.navigation_data.NavigationPoint): instance of `NavigationPoint` used
            to compute the satellite positions
        epoch (src.data_types.date.Epoch): epoch to compute the satellite ephemerides (must be the
            transmission time epoch, that is, reception time minus transit time)
    Returns:
        tuple[numpy.ndarray,numpy.ndarray,float,float]: Position vector in [m],
            velocity vector in [m/s], relativistic clock correction in [s] and relativistic clock drift in [s/s]
    """
    GM = constants.MU_WGS84 if nav_message.constellation == "GPS" else constants.MU_GTRF

    # unpack week and seconds of week
    _, sow = epoch.gnss_time

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
    position = array([x_ECEF, y_ECEF, z_ECEF])
    velocity = array([vx_ECEF, vy_ECEF, vz_ECEF])

    # compute relativistic correction. See Equation (19.15) of [2]
    rel_correction = -2 * sqrt(GM) * sqrtA / constants.SPEED_OF_LIGHT ** 2 * eccentricity * sin(E)

    # compute time derivative of relativistic correction (drift). See Equation (19.17) of [2]
    r = linalg.norm(position)
    rel_drift = 2 * GM / constants.SPEED_OF_LIGHT ** 2 * (1.0 / A - 1.0 / r)

    return position, velocity, rel_correction, rel_drift
