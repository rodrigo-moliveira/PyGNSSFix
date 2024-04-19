"""Module with several utility functions regarding Satellite Clock Biases
"""

from datetime import timedelta
from numpy.linalg import norm

from src import constants
from src.data_types.gnss import get_data_type
from src.io.config import EnumTransmissionTime
from src.common_log import get_logger, MODEL_LOG
from src.models.frames import dcm_e_i
from .ephemeride_propagator import EphemeridePropagator, fix_gnss_week_crossovers


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


def compute_tx_time(model=None, **kwargs):
    """
    Main function to compute the transmission time of the satellite-receiver GNSS link.
    There are two alternatives, depending on the user configuration:
        * Geometric algorithm: Purely Geometric Algorithm from section 5.1.1.1 of [1]
        * Pseudorange-based algorithm: A Pseudorange-Based Algorithm from section 5.1.1.1 of [1]

    Reference:
        [1] ESA GNSS DATA PROCESSING, Volume I: Fundamentals and Algorithms, J. Sanz Subirana,
        J.M. Juan Zornoza and M. Hernández-Pajares

    Args:
        model(EnumTransmissionTime): enumeration that evaluates to `EnumTransmissionTime.GEOMETRIC` or
            `EnumTransmissionTime.PSEUDORANGE`
        kwargs(dict): depending on the model, different items are provided in `kwargs`:
            * Geometric algorithm: `r_receiver`, `t_reception`, `dt_receiver`, `nav_message`.
            * Pseudorange-based algorithm: `pseudorange_obs`, `t_reception`, `nav_message`.
    Returns:
        tuple [src.data_types.date.date.Epoch, float] : computed TX epoch, computed transit time
    """
    if model == EnumTransmissionTime.GEOMETRIC:
        return tx_time_geometric(**kwargs)
    elif model == EnumTransmissionTime.PSEUDORANGE:
        return tx_time_pseudorange(**kwargs)


def tx_time_geometric(r_receiver=None, t_reception=None, dt_receiver=None, nav_message=None, **_):
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
        t_reception (src.data_types.date.date.Epoch): time of signal reception, measured by receiver
            `t(reception)^{receiver}`, i.e., RINEX obs time tag
        dt_receiver (float): current estimate of the receiver clock bias (seconds)
        nav_message (src.data_mng.gnss.navigation_data.NavigationPoint): instance of `NavigationPoint` used
            to compute the satellite positions
    Return:
        tuple [src.data_types.basics.Epoch.Epoch, float] : computed TX epoch, computed transit time
    """
    max_iter = 5
    residual_th = 1E-8
    rho_previous = 0
    residual = 1
    N = 0

    # 1. Initialize tau
    tau = 0

    while residual > residual_th and N < max_iter:
        # 2. Get satellite coordinates
        t = t_reception + timedelta(seconds=-tau)
        r_satellite, _, _ = EphemeridePropagator.compute_nav_sat_eph(nav_message, t.gnss_time)

        # 3. Compute pseudorange (in ECEF frame associated to t_receiver epoch)
        _R = dcm_e_i(-tau)
        rho = norm(_R @ r_satellite - r_receiver)
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


def tx_time_pseudorange(pseudorange_obs=None, t_reception=None, nav_message=None, **_):
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
        t_reception (src.data_types.date.date.Epoch): time of signal reception, measured by receiver
            `t(reception)^{receiver}`, i.e., RINEX obs time tag
        nav_message (src.data_mng.gnss.navigation_data.NavigationPoint): instance of `NavigationPoint` used
            to compute the satellite positions
    Return:
        tuple [src.data_types.basics.Epoch.Epoch, float] : computed TX epoch, computed transit time
    """

    tau = pseudorange_obs.value / constants.SPEED_OF_LIGHT
    t_emission = t_reception + timedelta(seconds=-tau)  # t(emission)^{satellite} = t(reception)^{receiver} - tau
    dt_sat, _ = broadcast_clock(nav_message.af0, nav_message.af1, nav_message.af2,
                                nav_message.toc.gnss_time[1],  # to get seconds of week
                                t_emission.gnss_time[1]  # to get seconds of week
                                )

    # correct satellite clock for BGDs
    dt_sat = nav_sat_clock_correction(dt_sat, pseudorange_obs.datatype, nav_message)

    # compute transmission time, using Eq. (5.6)
    T_emission = t_emission + timedelta(seconds=-dt_sat)  # t(emission)^{GPS} = t(emission)^{satellite} - dt_sat
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
    Return:
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
        This is currently not considered in this algorithm!!!

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
            to fetch to BGDs/TGD
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
        time_correction(dict): dictionary with the `GGTO` item.
            See :py:class:`src.data_mng.gnss.navigation_data.NavigationData`
        week(int): the week number of the user epoch to compute the GGTO
        sow(float): the SoW (seconds of week) of the user epoch to compute the GGTO
    """
    epoch = (week, sow)

    if "GGTO" not in time_correction:
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
