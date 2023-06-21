from numpy.linalg import norm

from PositioningSolver.src.gnss.observation_models import ephemeride_propagator
from PositioningSolver.src.gnss.state_space.utils import matrix_ECEF2ECI
from PositioningSolver.src.math_utils.Constants import Constant

# utility functions related to navigation clocks


def compute_TX_time_geometric(r_receiver=None, t_reception=None, dt_receiver=None, nav_message=None, **kwargs):
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
        t = t_reception + (-tau)
        r_satellite, _ = ephemeride_propagator.EphemeridePropagator.compute(nav_message, t, False)

        # 3. Compute pseudorange (in ECEF frame associated to t_receiver epoch)
        _R = matrix_ECEF2ECI(-tau)
        rho = norm(_R @ r_satellite - r_receiver)
        # LOS = (_R @ r_satellite - r_receiver) / rho

        # 4. Compute new tau [s]
        tau = rho / Constant.SPEED_OF_LIGHT

        # 5. Check convergence
        residual = abs(rho - rho_previous)
        rho_previous = rho
        N += 1

    # compute transmission time, using Eq. (5.9)
    T_emission = t_reception + (-tau) + (-dt_receiver)
    return T_emission, tau


def compute_TX_time_pseudorange(pseudorange_obs=None, t_reception=None, nav_message=None, TGD=None, **kwargs):
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
        TGD (float): time group delay (TGD) to correct the satellite hardware clock delay to the appropriate frequency
        nav_message (src.data_types.containers.NavigationData.NavigationPointGPS): ephemeride object to use
    Return:
        tuple [src.data_types.basics.Epoch.Epoch, float] : computed TX epoch, computed transit time
    """

    tau = pseudorange_obs.value / Constant.SPEED_OF_LIGHT
    t_emission = t_reception + (-tau)  # t(emission)^{satellite} = t(reception)^{receiver} - tau
    dt_sat, _ = SVBroadcastCorrection(nav_message.af0, nav_message.af1, nav_message.af2, nav_message.toc, t_emission)

    # correct for TGD (already corrected for the appropriate frequency, and is 0 for IF observables)
    dt_sat = dt_sat - TGD

    # compute transmission time, using Eq. (5.6)
    T_emission = t_emission + (-dt_sat)  # t(emission)^{GPS} = t(emission)^{satellite} - dt_sat
    return T_emission, tau


def SVBroadcastCorrection(af0: float, af1: float, af2: float, toc: float, tx_raw: float) -> tuple[float, float]:
    """
    The polynomial defined in the following allows the user to determine the effective SV PRN code phase offset
    referenced to the phase center of the antennas (Δt sv) with respect to GNSS system time (t) at the time of
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

    dt = ephemeride_propagator.correct_gps_week_crossovers(tx_raw - toc)
    clock_bias = af0 + af1 * dt + af2 * dt * dt
    clock_drift = af1 + 2 * af2 * dt

    return clock_bias, clock_drift
