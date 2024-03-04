from datetime import timedelta

from numpy.linalg import norm

from src.data_types.gnss.data_type import *
from src.io.config.enums import EnumTransmissionTime
from src.models.frames.frames import dcm_e_i
from src.models.gnss_obs.ephemeride_propagator import EphemeridePropagator, fix_gnss_week_crossovers
from src import constants


# utility functions related to navigation clocks

_ggto_cache = {}
L1 = get_data_type("L1", "GPS")
L2 = get_data_type("L2", "GPS")

E1 = get_data_type("E1", "GAL")
E5a = get_data_type("E5a", "GAL")
E5b = get_data_type("E5b", "GAL")


def compute_tx_time(model=None, **kwargs):
    if model == EnumTransmissionTime.GEOMETRIC:
        return tx_time_geometric(**kwargs)
    elif model == EnumTransmissionTime.PSEUDORANGE:
        return tx_time_pseudorange(**kwargs)


def tx_time_geometric(r_receiver=None, t_reception=None, dt_receiver=None, nav_message=None, **kwargs):
    """
    Implements the algorithm to compute:
        * transit time (propagation time from signal transmission to reception)
        * Transmission time in GNSS time system

    Algorithm used --- Purely Geometric Algorithm from [section 5.1.1.1 of **REF[1]**]

    Fundamental equation
        t(emission)^{receiver} = t(reception)^{receiver} - tau

    Procedure:
        1. tau = 0
        2. get satellite position r_satellite at t = t(emission)^{receiver} = t(reception)^{receiver} - tau
        3. get pseudorange rho = norm(R(-tau) @ r_satellite - r_receiver)
        4. tau = rho / c
        5. go to step 2 and iterate until convergence
    where R() is the rotation matrix from ECEF to ECI, and c is the speed of light

    Args:
        r_receiver (src.data_types.state_space.statevector.Position): position of receiver at reception time
                                                                (RINEX time tag measured by receiver clock)
        t_reception (src.data_types.basics.Epoch.Epoch): time of signal reception, measured by receiver
                                                        (receiver clock), i.e., RINEX obs time tag
        dt_receiver (float): receiver clock bias (seconds)
        nav_message (src.data_types.containers.NavigationData.NavigationPointGPS): ephemeride object to use
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


def tx_time_pseudorange(pseudorange_obs=None, t_reception=None, nav_message=None, **kwargs):
    """
    Implements the algorithm to compute:
        * transit time (propagation time from signal transmission to reception)
        * Transmission time in GNSS time system

    Algorithm used --- A Pseudorange-Based Algorithm from [section 5.1.1.1 of **REF[1]**]

    Fundamental equation
        pseudorange = c (t(reception)^{receiver} - t(emission)^{satellite})
        -> tau = t(reception)^{receiver} - t(emission)^{satellite}

    Args:
        pseudorange_obs (src.data_types.data_types.Observation.Observation): the pseudorange measured observable
        t_reception (src.data_types.basics.Epoch.Epoch): time of signal reception, measured by receiver
                                                (receiver clock), i.e., RINEX obs time tag -> t(reception)^{receiver}
        nav_message (src.data_types.containers.NavigationData.NavigationPointGPS): ephemeride object to use
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
    The polynomial defined in the following allows the user to determine the effective SV PRN code phase offset
    referenced to the phase center of the antennas (delta t sv) with respect to GNSS system time (t) at the time of
    inputs transmission.

    This estimated correction accounts for the deterministic SV clock error characteristics of bias, drift and
    aging, as well as for the SV implementation characteristics of group delay bias and mean differential group delay.
    Since these coefficients do not include corrections for relativistic effects, the user's equipment must determine
    the requisite relativistic correction.

    According to section [20.3.3.3.3.1] of **REF[3]**

    computes satellite clock bias and drift, using the parameters broadcast in the Navigation Message

    Δt sv = af0 + af1(t - toc) + af2(t - toc)^2 + Δtr

    Args:
        af0 (float) : clock model constant parameter [sec]
        af1 (float) : clock model first order parameter [sec/sec]
        af2 (float) : clock model second order parameter [sec/sec^2]
        toc (float) : time of clock (reference epoch for clock model) [seconds of week]
        tx_raw (float) : epoch when bias and drift are desired (time of signal transmission) [seconds of week]
    Return:
        tuple [float, float] : clock bias, clock drift
    """
    dt = fix_gnss_week_crossovers(tx_raw - toc)
    clock_bias = af0 + af1 * dt + af2 * dt * dt
    clock_drift = af1 + 2 * af2 * dt

    return clock_bias, clock_drift


def nav_sat_clock_correction(sat_clock, datatype, nav_message):
    bgd = get_bgd_correction(datatype, nav_message)
    return sat_clock - bgd


def get_bgd_correction(datatype, nav_message):
    # NOTE: this function is only called assuming standard SPS Mode, that is, the only BGDs available
    # are the ones from the navigation message. No extra DCB information is used
    if nav_message.constellation == "GPS" and datatype.constellation == "GPS":
        if datatype.freq_number == 1 or datatype.freq_number == 2:
            scale = (L1.freq_value / datatype.freq.freq_value) ** 2
            return scale * nav_message.TGD
        elif datatype.freq_number == 12:
            # the broadcast clock refers to the L1-L2 iono-free clock -> no correction needed
            return 0.0
        else:
            # TODO: for all other cases (PR5 or PR15) issue warning
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
                return 0.0  # TODO: issue warning
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
                return 0.0  # TODO: issue warning

    else:
        raise AttributeError(f"Mismatch between navigation message constellation ({nav_message.constellation}) and "
                             f" datatype {repr(datatype)}")
    return 0.0


def compute_ggto(time_correction, epoch):
    if epoch not in _ggto_cache:
        week, sow = epoch.gnss_time
        # fetch data
        a0 = time_correction["GGTO"][0]
        a1 = time_correction["GGTO"][1]
        sow_ref = time_correction["GGTO"][2]
        week_ref = time_correction["GGTO"][3]

        # apply this eq. GGTO(t_sow) = a0 + a1 * (t_sow - T + 604800 * (week - Week_ref) )
        _ggto_cache[epoch] = a0 + a1 * (sow - sow_ref + constants.SECONDS_IN_GNSS_WEEK * (week - week_ref))
    return _ggto_cache[epoch]
