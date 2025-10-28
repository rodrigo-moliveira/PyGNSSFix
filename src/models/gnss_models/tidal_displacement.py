""" Module with computation of Earth Tides Displacements. Based on IERS Conventions 2010 """

from math import sin, cos
import numpy as np

from src.constants import MASS_RATIO_SUN, MASS_RATIO_MOON, SECONDS_IN_DAY, AVERAGE_DAYS_IN_YEAR, PI, DEG2RAD
from src.data_types.date import Epoch
from src.data_types.date.frames import get_gmst
from src.io.config import config_dict
from src.models.frames import latlon2dcm_e_enu
from src.spicepy_wrapper import compute_sun_pos, compute_moon_pos

# specific constants for the IERS 2010 convention
R_EARTH = 6378136.6
H2 = 0.6078
L2 = 0.0847
H3 = 0.292
L3 = 0.015

OCEAN_LOADING_ARGS = [[1.40519E-4, 2.0, -2.0, 0.0, 0.00],   # M2
                      [1.45444E-4, 0.0, 0.0, 0.0, 0.00],    # S2
                      [1.37880E-4, 2.0, -3.0, 1.0, 0.00],   # N2
                      [1.45842E-4, 2.0, 0.0, 0.0, 0.00],    # K2
                      [0.72921E-4, 1.0, 0.0, 0.0, 0.25],    # K1
                      [0.67598E-4, 1.0, -2.0, 0.0, -0.25],  # O1
                      [0.72523E-4, -1.0, 0.0, 0.0, -0.25],  # P1
                      [0.64959E-4, 1.0, -3.0, 1.0, -0.25],  # Q1
                      [0.53234E-5, 0.0, 2.0, 0.0, 0.00],    # Mf
                      [0.26392E-5, 0.0, 1.0, -1.0, 0.00],   # Mm
                      [0.03982E-5, 2.0, 0.0, 0.0, 0.00]]    # Ssa


def tide_gravity(h2, l2, r_rec_unit, r_body, MASS_RATIO, lat, long):
    """
    Computes the solar/lunar planetary contributions to the Earth tides.
    This function implements the simplified In-phase and out-of-phase (Step 1) corrections from the IERS Conventions:
        * in-phase for degree 2 and 3
        * out of phase (only radial component) for degree 2

    Reference:
        [1] IERS Conventions 2010 https://iers-conventions.obspm.fr/content/tn36.pdf (Chapter 7)

    Args:
        h2 (float)
        l2 (float)
        r_rec_unit (numpy.ndarray)
        r_body (numpy.ndarray)
        MASS_RATIO (float)
        lat (float)
        long (float)

    Returns:
        numpy.ndarray : computed correction in ENU frame (m)
    """
    r_body_norm = np.linalg.norm(r_body)
    r_body_unit = r_body / r_body_norm
    lat_body = np.arcsin(r_body[2] / np.linalg.norm(r_body))
    long_body = np.arctan2(r_body[1], r_body[0])

    dot_prod = np.dot(r_body_unit, r_rec_unit)
    K2 = MASS_RATIO * R_EARTH ** 4 / r_body_norm ** 3
    K3 = MASS_RATIO * R_EARTH ** 5 / r_body_norm ** 4

    # Step 1: in-phase for degree 2 and 3
    #   -> h2 and l2 terms for Eq. (7.5)
    h2_term = h2 * r_rec_unit * (1.5 * dot_prod ** 2 - 0.5)
    l2_term = 3 * l2 * dot_prod * (r_body_unit - dot_prod * r_rec_unit)

    #   -> h3 and l3 terms for Eq. (7.6)
    h3_term = H3 * r_rec_unit * (2.5 * dot_prod ** 3 - 1.5 * dot_prod)
    l3_term = L3 * (7.5 * dot_prod ** 2 - 1.5) * (r_body_unit - dot_prod * r_rec_unit)

    # Eq. (7.5) + Eq. (7.6)
    disp_body = K2 * (h2_term + l2_term) + K3 * (h3_term + l3_term)

    # Step 1: out of phase (only radial)
    du = 3.0 / 4.0 * 0.0025 * K2 * sin(2.0 * lat_body) * sin(2.0 * lat) * sin(long - long_body)  # Eq. (7.10a)
    du += 3.0 / 4.0 * 0.0022 * K2 * cos(lat_body) ** 2 * np.cos(lat) ** 2 * sin(2.0 * (long - long_body))  # Eq. (7.11a)
    disp_body += du * r_rec_unit

    return disp_body


def solid_disp(lat, long, r_sun, r_moon, r_rec, R_ECEF2ENU, gmst, step3=True):
    """
    Displacement by solid earth tide.

    This function implements the simplified corrections from the IERS Conventions 2010:
        * Step 1: Corrections to be computed in the time domain
        * Step 2: Corrections to be computed in the frequency domain and to be added to the results of Step 1
        * Step 3: Eliminate permanent deformation

    Reference:
        [1] IERS Conventions 2010 https://iers-conventions.obspm.fr/content/tn36.pdf (Chapter 7)

    Args:
        lat (float)
        long (float)
        r_sun (numpy.ndarray)
        r_moon (numpy.ndarray)
        r_rec (numpy.ndarray)
        R_ECEF2ENU (numpy.ndarray)
        gmst (float)
        step3 (bool)

    Returns:
        numpy.ndarray : computed correction in ECEF frame (m)
    """
    r_rec_unit = r_rec / np.linalg.norm(r_rec)

    h2 = H2 - 0.0006 * ((3 * sin(lat) ** 2 - 1) / 2)
    l2 = L2 + 0.0002 * ((3 * sin(lat) ** 2 - 1) / 2)

    # step 1: time domain
    # Sun and moon contributions
    disp_sun = tide_gravity(h2, l2, r_rec_unit, r_sun, MASS_RATIO_SUN, lat, long)
    disp_moon = tide_gravity(h2, l2, r_rec_unit, r_moon, MASS_RATIO_MOON, lat, long)
    disp = disp_sun + disp_moon

    # step 2: frequency domain
    sin2l = sin(2.0 * lat)
    du = -0.012 * sin2l * sin(gmst + long)  # Eq. (7.12a) only for K1 radial
    delta_enu = np.array([0.0, 0.0, du])
    disp = disp + R_ECEF2ENU.T @ delta_enu

    if step3:
        # step 3: eliminate permanent deformation
        P2 = (3 * np.sin(lat) ** 2 - 1) / 2
        d_radial = (-0.1206 + 0.0001 * P2) * P2  # Eq. (7.14a)
        d_north = (-0.0252 + 0.0001 * P2) * np.sin(2 * lat)  # Eq. (7.14b)
        delta_enu = np.array([0.0, d_north, d_radial])
        delta_ecef = R_ECEF2ENU.T @ delta_enu
        disp = disp - delta_ecef

    return disp


def iers_mean_pole(epoch):
    """
    This function computes the Mean Pole (xp_bar and yp_bar) from the IERS Conventions

    It implements Eq. 7.25 of the Technical Note:
        https://iers-conventions.obspm.fr/content/tn36.pdf

    Args:
        epoch(src.data_types.date.Epoch): epoch instance to compute the corrections

    Returns:
        computed xp_bar and yp_bar angles in (mas)
    """
    epoch_ut1 = epoch.change_scale("UT1").datetime
    epoch_2000 = Epoch(2000, 1, 1, 0, 0, 0, scale='UT1').datetime
    years = (epoch_ut1 - epoch_2000).total_seconds() / SECONDS_IN_DAY / AVERAGE_DAYS_IN_YEAR

    if years < 3653.0 / AVERAGE_DAYS_IN_YEAR:  # until 2010.0
        y2 = years * years
        y3 = y2 * years
        xp_bar = 55.974 + 1.8243 * years + 0.18413 * y2 + 0.007024 * y3  # in (mas)
        yp_bar = 346.346 + 1.7896 * years - 0.10729 * y2 - 0.000908 * y3
    else:  # after 2010.0
        xp_bar = 23.513 + 7.6141 * years  # (mas)
        yp_bar = 358.891 - 0.6287 * years

    return xp_bar, yp_bar


def tide_pole(epoch, lat, long):
    """
    Displacement by pole tide

    Note: for this function the correct formula is presented in the Errata of the technical note in:
        https://iers-conventions.obspm.fr/content/chapter7/icc7.pdf

    Args:
        epoch(src.data_types.date.Epoch): epoch instance to compute the correction
        lat(float)
        long(float)

    Returns:
        numpy.ndarray : computed correction in ENU frame (m)
    """
    disp_enu = [0.0, 0.0, 0.0]

    # iers mean pole (mas)
    xp_bar, yp_bar = iers_mean_pole(epoch)

    # earth rotation parameters
    eop = epoch.eop
    xp = eop.x  # in as
    yp = eop.y  # in as

    # ref Eq. (7.24)
    m1 = xp - xp_bar * 1E-3  # (as)
    m2 = -yp + yp_bar * 1E-3

    cos_long = cos(long)
    sin_long = sin(long)

    # ref Eq. (7.26) (of the Errata!!)
    # sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi)
    disp_enu[0] = 9E-3 * sin(lat) * (m1 * sin_long - m2 * cos_long)  # de = S_lambda (m)  positive east
    disp_enu[1] = -9E-3 * cos(2.0 * lat) * (m1 * cos_long + m2 * sin_long)  # dn = -S_theta (m)  positive south
    disp_enu[2] = -33E-3 * sin(2.0 * lat) * (m1 * cos_long + m2 * sin_long)  # du = S_r (m)  positive upwards

    return disp_enu


def ocean_loading(epoch, ocean_mng):
    """
    This function computes the time series of tidal displacement from the ocean tide loading.
    It requires an input file containing the ocean loading coefficients for a given station. These coefficients can be
    obtained from the ocean loading service by request from the website https://barre.oso.chalmers.se/loading/l.php

    The equations for this model can be seen in the following reference.
    Reference:
        [1]: https://apps.dtic.mil/sti/tr/pdf/ADA112278.pdf

    Args:
        epoch(src.data_types.date.Epoch): epoch instance to compute the correction
        ocean_mng(src.data_mng.gnss.ocean_loading_data.OceanLoadingData): ocean loading manager

    Returns:
        numpy.ndarray : computed correction in ENU frame (m)
    """

    # Epoch computations
    epoch_ut1 = epoch.change_scale("UT1").datetime
    epoch_ut1_00 = Epoch(epoch_ut1.year, epoch_ut1.month, epoch_ut1.day, 0, 0, 0, scale="UT1").datetime
    epoch_1975 = Epoch(1975, 1, 1, 0, 0, 0, scale="UT1").datetime

    f_day = epoch_ut1.hour * 3600.0 + epoch_ut1.minute * 60.0 + epoch_ut1.second
    d_days = (epoch_ut1_00 - epoch_1975).total_seconds() / 86400 + 1.0
    t = (27392.500528 + 1.000000035 * d_days) / 36525.0
    t2 = t * t
    t3 = t2 * t

    a = [0, 0, 0, 0, 0]
    a[0] = f_day
    a[1] = (279.69668 + 36000.768930485 * t + 3.03E-4 * t2) * DEG2RAD  # H0
    a[2] = (270.434358 + 481267.88314137 * t - 0.001133 * t2 + 1.9E-6 * t3) * DEG2RAD  # S0
    a[3] = (334.329653 + 4069.0340329577 * t - 0.010325 * t2 - 1.2E-5 * t3) * DEG2RAD  # P0
    a[4] = 2 * PI

    disp = [0.0, 0.0, 0.0]  # radial, west, south
    d_enu = [0.0, 0.0, 0.0]  # east, north, up

    for i in range(11):
        ang = 0.0
        for j in range(5):
            ang += a[j] * OCEAN_LOADING_ARGS[i][j]
        disp[0] += ocean_mng.radial_amplitude[i] * cos(ang - ocean_mng.radial_phase[i] * DEG2RAD)
        disp[1] += ocean_mng.west_amplitude[i] * cos(ang - ocean_mng.west_phase[i] * DEG2RAD)
        disp[2] += ocean_mng.south_amplitude[i] * cos(ang - ocean_mng.south_phase[i] * DEG2RAD)

    d_enu[0] = -disp[1]  # east (=-west)
    d_enu[1] = -disp[2]  # north(=-south)
    d_enu[2] = disp[0]  # up (=+radial)
    return np.array(d_enu)


def compute_displacement(epoch, r_rec, ocean_loading_mng=None, gmst_model='IAU82'):
    """
    Main function to compute displacement corrections by earth tides. These corrections are applied to the station
    positions, or directly as a correction in the PR / CP measurements (when multiplied by the LoS for each GNSS SV
    link)

    The tidal model is based on the IERS Conventions 2010 technical note:
        * https://iers-conventions.obspm.fr/conventions_material.php
        * Technical Note: https://iers-conventions.obspm.fr/content/tn36.pdf

    The implemented tides are:
        * Solid Earth Tide (with elimination of permanent deformation)
        * Ocean Tide Loading (Ocean loading file is required)
        * Pole Tide

    Args:
        epoch(src.data_types.date.Epoch): epoch instance to compute the corrections
        r_rec(numpy.ndarray): receiver position array in ECEF (m)
        ocean_loading_mng(src.data_mng.gnss.ocean_loading_data.OceanLoadingData): ocean loading manager
        gmst_model(str): Model to compute Greenwich Mean Sidereal Time (IAU82, IAU00 or IAU06)

    Returns:
        numpy.ndarray : computed tidal displacement in ECEF frame
    """
    # Main function
    disp_ecef = np.array([0.0, 0.0, 0.0])  # site displacement vector in ECEF

    enable_solid_tides = config_dict.get("model", "earth_deformation_effects", "enable_solid_tides")
    enable_ocean_loading = config_dict.get("model", "earth_deformation_effects", "enable_ocean_loading")
    enable_pode_tides = config_dict.get("model", "earth_deformation_effects", "enable_pode_tides")

    # pre-computations required
    lat = np.arcsin(r_rec[2] / np.linalg.norm(r_rec))
    long = np.arctan2(r_rec[1], r_rec[0])
    R_ECEF2ENU = latlon2dcm_e_enu(lat, long)

    # Solid Displacement
    if enable_solid_tides:
        gmst = get_gmst(epoch, gmst_model)
        r_sun = compute_sun_pos(epoch)
        r_moon = compute_moon_pos(epoch)
        disp_solid = solid_disp(lat, long, r_sun, r_moon, r_rec, R_ECEF2ENU, gmst, step3=True)
        disp_ecef += disp_solid

    # ocean loading
    if enable_ocean_loading and ocean_loading_mng is not None:
        disp_ocn = ocean_loading(epoch, ocean_loading_mng)  # in ENU
        disp_ecef += R_ECEF2ENU.T @ disp_ocn

    # displacement by pole tide
    if enable_pode_tides:
        disp_pole = tide_pole(epoch, lat, long)  # in ENU
        disp_ecef += R_ECEF2ENU.T @ disp_pole

    return disp_ecef
