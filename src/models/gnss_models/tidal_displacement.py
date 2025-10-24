# see https://github.com/tomojitakasu/RTKLIB/blob/71db0ffa0d9735697c6adfd06fdf766d0e5ce807/src/ppp.c#L184
# see also file:///C:/Users/rooo/Downloads/iersconventions_v1_3_0/tn36.pdf
# see also https://www.rtklib.com/prog/manual_2.4.2.pdf

# TODO: add docstrings

from math import sin, cos
import numpy as np

from src.constants import MASS_RATIO_SUN, MASS_RATIO_MOON, SECONDS_IN_DAY, AVERAGE_DAYS_IN_YEAR
from src.data_types.date import Epoch
from src.data_types.date.frames import get_gmst
from src.models.frames import latlon2dcm_e_enu

R_EARTH = 6378136.6  # specific value for the IERS 2010 convention
H2 = 0.6078
L2 = 0.0847
H3 = 0.292
L3 = 0.015


def tide_gravity(h2, l2, r_rec_unit, r_body, MASS_RATIO, lat, long):
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
    du += 3.0 / 4.0 * 0.0022 * K2 * cos(lat_body)**2 * np.cos(lat)**2 * sin(2.0 * (long - long_body))  # Eq. (7.11a)
    disp_body += du * r_rec_unit

    return disp_body


def solid_disp(lat, long, r_sun, r_moon, r_rec, R_ECEF2ENU, gmst, step3):
    r_rec_unit = r_rec / np.linalg.norm(r_rec)

    h2 = H2 - 0.0006 * ((3 * sin(lat)**2 - 1) / 2)
    l2 = L2 + 0.0002 * ((3 * sin(lat)**2 - 1) / 2)

    # step 1: time domain
    # Sun and moon contributions
    disp_sun = tide_gravity(h2, l2, r_rec_unit, r_sun, MASS_RATIO_SUN, lat, long)
    disp_moon = tide_gravity(h2, l2, r_rec_unit, r_moon, MASS_RATIO_MOON, lat, long)
    disp = disp_sun + disp_moon

    # step 2: frequency domain
    sin2l = sin(2.0 * lat)
    du = -0.012 * sin2l * sin(gmst + long)  # Eq. (7.12a) only for K1 radial
    delta_enu = np.array([0, 0, du])
    disp = disp + R_ECEF2ENU.T @ delta_enu

    if step3:
        # step 3: eliminate permanent deformation
        P2 = (3 * np.sin(lat) ** 2 - 1) / 2
        d_radial = (-0.1206 + 0.0001 * P2) * P2  # Eq. (7.14a)
        d_north = (-0.0252 + 0.0001 * P2) * np.sin(2*lat)  # Eq. (7.14b)
        delta_enu = np.array([0, d_north, d_radial])
        delta_ecef = R_ECEF2ENU.T @ delta_enu
        disp = disp - delta_ecef

    return disp


def iers_mean_pole(epoch):
    epoch_ut1 = epoch.change_scale("UT1")
    epoch_2000 = Epoch(2000, 1, 1, 0, 0, 0, scale='UT1')
    years = (epoch_ut1 - epoch_2000).total_seconds() / SECONDS_IN_DAY / AVERAGE_DAYS_IN_YEAR

    if years < 3653.0/AVERAGE_DAYS_IN_YEAR:  # until 2010.0
        y2 = years * years
        y3 = y2 * years
        xp_bar = 55.974 + 1.8243 * years + 0.18413 * y2 + 0.007024 * y3  # in (mas)
        yp_bar = 346.346 + 1.7896 * years - 0.10729 * y2 - 0.000908 * y3
    else:  # after 2010.0
        xp_bar = 23.513 + 7.6141 * years  # (mas)
        yp_bar = 358.891 - 0.6287 * years

    return xp_bar, yp_bar


def tide_pole(epoch, lat, long):
    disp_enu = [0, 0, 0]

    # iers mean pole (mas)
    xp_bar, yp_bar = iers_mean_pole(epoch)

    # earth rotation parameters
    eop = epoch.eop
    xp = eop.x  # in as
    yp = eop.y  # in as

    # ref [7] eq.7.24
    m1 = xp - xp_bar * 1E-3  # (as)
    m2 = -yp + yp_bar * 1E-3

    cos_long = cos(long)
    sin_long = sin(long)

    # ref [7] eq.7.26 (correct formula is from the Errata
    # https://iers-conventions.obspm.fr/content/chapter7/icc7.pdf
    # sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi)
    disp_enu[0] = 9E-3 * sin(lat) * (m1 * sin_long - m2 * cos_long)  # de = S_lambda (m)  positive east
    disp_enu[1] = -9E-3 * cos(2.0 * lat) * (m1 * cos_long + m2 * sin_long)  # dn = -S_theta (m)  positive south
    disp_enu[2] = -33E-3 * sin(2.0 * lat) * (m1 * cos_long + m2 * sin_long)  # du = S_r (m)  positive upwards

    return disp_enu


def ocean_loading():
    """This routine computes time series of tidal displacements
    from an input file containing the ocean loading coefficients for a given station. These coefficients can be obtained from the ocean loading service by
    request from the website https://barre.oso.chalmers.se/loading/l.php"""
    pass


def compute_displacement(epoch, r_sun, r_moon, r_rec, step3=True, gmst_model='IAU82'):
    # site displacement vector in ECEF
    lat = np.arcsin(r_rec[2] / np.linalg.norm(r_rec))
    long = np.arctan2(r_rec[1], r_rec[0])

    # pre-computations required
    R_ECEF2ENU = latlon2dcm_e_enu(lat, long)
    gmst = get_gmst(epoch, gmst_model)

    # Solid Displacement
    disp = solid_disp(lat, long, r_sun, r_moon, r_rec, R_ECEF2ENU, gmst, step3)

    # displacement by pole tide
    disp_enu = tide_pole(epoch, lat, long)
    disp_pole = R_ECEF2ENU.T @ disp_enu
    disp = disp + disp_pole

    return disp
