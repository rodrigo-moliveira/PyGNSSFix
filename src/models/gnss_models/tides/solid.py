# see https://github.com/tomojitakasu/RTKLIB/blob/71db0ffa0d9735697c6adfd06fdf766d0e5ce807/src/ppp.c#L184
# see also file:///C:/Users/rooo/Downloads/iersconventions_v1_3_0/tn36.pdf
# see also https://www.rtklib.com/prog/manual_2.4.2.pdf

from math import sin
import numpy as np

from src.constants import MASS_RATIO_SUN, MASS_RATIO_MOON
from src.models.frames import latlon2dcm_e_enu

R_EARTH = 6378136.6  # specific value for the IERS 2010 convention
H2 = 0.6078
L2 = 0.0847
H3 = 0.292
L3 = 0.015


def tide_gravity(h2, l2, r_rec_unit, r_body, MASS_RATIO):
    # order 2 and 3
    r_body_norm = np.linalg.norm(r_body)
    r_body_unit = r_body / r_body_norm
    dot_prod = np.dot(r_body_unit, r_rec_unit)
    h2_term = h2 * r_rec_unit * (1.5 * dot_prod ** 2 - 0.5)
    l2_term = 3 * l2 * dot_prod * (r_body_unit - dot_prod * r_rec_unit)
    h3_term = H3 * r_rec_unit * (2.5 * dot_prod ** 3 - 1.5 * dot_prod)
    l3_term = L3 * (7.5 * dot_prod ** 2 - 1.5) * (r_body_unit - dot_prod * r_rec_unit)
    disp_body = MASS_RATIO * R_EARTH ** 4 / r_body_norm ** 3 * (h2_term + l2_term) + \
                MASS_RATIO * R_EARTH ** 5 / r_body_norm ** 4 * (h3_term + l3_term)
    return disp_body


def solid_disp(r_sun, r_moon, r_rec):
    r_rec_unit = r_rec / np.linalg.norm(r_rec)
    lat = np.arcsin(r_rec_unit[2])
    long = np.arctan2(r_rec[1], r_rec[0])

    h2 = H2 - 0.0006 * ((3 * sin(lat)**2 - 1) / 2)
    l2 = L2 + 0.0002 * ((3 * sin(lat)**2 - 1) / 2)

    # step 1: time domain
    # Sun and moon contributions
    disp_sun = tide_gravity(h2, l2, r_rec_unit, r_sun, MASS_RATIO_SUN)
    disp_moon = tide_gravity(h2, l2, r_rec_unit, r_moon, MASS_RATIO_MOON)
    disp = disp_sun + disp_moon

    # step 2: frequency domain
    # not yet done...

    # step 3: eliminate permanent deformation
    P2 = (3 * np.sin(lat) ** 2 - 1) / 2
    d_radial = (-0.1206 + 0.0001 * P2) * P2
    d_north = (-0.0252 + 0.0001 * P2) * np.sin(2*lat)

    delta_radial = d_radial * r_rec_unit
    delta_enu = np.array([0, d_north, 0])
    R_ECEF2ENU = latlon2dcm_e_enu(lat, long)
    delta_ecef = delta_radial + R_ECEF2ENU.T @ delta_enu

    return disp - delta_ecef

def compute_displacement(los, r_sun, r_moon, r_rec):
    # site displacement vector in ECEF
    disp = solid_disp(r_sun, r_moon, r_rec)
    dr = np.dot(los, disp)
    return dr



def solid_testing():
    moon = np.array([-179996231.920342, -312468450.131567, -169288918.592160])
    sun = np.array([137859926952.015, 54228127881.4350, 23509422341.6960])
    rec = np.array([4075578.385, 931852.890, 4801570.154])
    disp = solid_disp(sun, moon, rec)
    iers_disp = [0.7700420357108125891E-01, 0.6304056321824967613E-01, 0.5516568152597246810E-01]
    print("disp", disp)
    print("iers_disp", iers_disp)


solid_testing()
def rtk1():
    # test 1
    moon = np.array([-179996231.920342, -312468450.131567, -169288918.592160])
    sun = np.array([137859926952.015, 54228127881.4350, 23509422341.6960])
    rec = np.array([4075578.385, 931852.890, 4801570.154])

    r_rec_unit = rec / np.linalg.norm(rec)
    lat = np.arcsin(r_rec_unit[2])
    long = np.arctan2(rec[1], rec[0])
    R_ECEF2ENU = latlon2dcm_e_enu(lat, long)

    disp = solid_disp(sun, moon, rec)
    print("disp", disp)

    # eliminate permanent deformation
    eu = np.array([R_ECEF2ENU[0, 2], R_ECEF2ENU[1, 2], R_ECEF2ENU[2, 2]])
    sinl = np.sin(lat)
    sin2l = np.sin(2.0 * lat)
    du = 0.1196 * (1.5 * sinl ** 2 - 0.5)
    dn = 0.0247 * sin2l
    delta1 = [0, 0, 0]
    delta1[0] = du * R_ECEF2ENU[2,0] + dn * R_ECEF2ENU[1,0]
    delta1[1] = du * R_ECEF2ENU[2,1] + dn * R_ECEF2ENU[1,1]
    delta1[2] = du * R_ECEF2ENU[2,2] + dn * R_ECEF2ENU[1,2]




    P2 = (3 * sinl**2 - 1) / 2
    du2 = (-0.1206 + 0.0001 * P2) * P2
    dn2 = (-0.0252 + 0.0001 * P2) * sin2l

    delta_r = du2 * r_rec_unit
    delta_enu = np.array([0, dn2, 0])
    delta_ecef = delta_r + R_ECEF2ENU.T @ delta_enu
    print("delta1", delta1)
    print("delta_ecef", delta_ecef[0], delta_ecef[1], delta_ecef[2])








########## TEST LATER ##########
# --- Constants (fill with your actual values) ---
RE_WGS84 = 6378137.0       # [m] Earth radius (WGS84)
GME = 3.986004415e14       # [m^3/s^2] Earth GM
GMS = 1.32712442099e20     # [m^3/s^2] Sun GM
GMM = 4.902801e12          # [m^3/s^2] Moon GM

def norm(v):
    return np.linalg.norm(v)

def dot(a, b):
    return float(np.dot(a, b))

def tide_pl(eu, rp, GMp, pos):
    """
    Planetary tide contribution (Sun or Moon).
    eu : unit vector of local Up (3,)
    rp : position vector of planet (Sun/Moon) in ECEF? (3,)
    GMp: GM of planet
    pos: [lat, lon, height?] geodetic latitude, longitude (rad)
    returns dr (3,)
    """
    H3, L3 = 0.292, 0.015
    r = norm(rp)
    if r <= 0.0:
        return np.zeros(3)

    ep = rp / r

    K2 = GMp / GME * (RE_WGS84**4) / (r**3)
    K3 = K2 * RE_WGS84 / r

    latp = np.arcsin(ep[2])
    lonp = np.arctan2(ep[1], ep[0])

    cosp = np.cos(latp)
    sinl = np.sin(pos[0])
    cosl = np.cos(pos[0])

    # step1 in phase (degree 2)
    p = (3.0 * sinl**2 - 1.0) / 2.0
    H2 = 0.6078 - 0.0006 * p
    L2 = 0.0847 + 0.0002 * p
    a = dot(ep, eu)
    dp = K2 * 3.0 * L2 * a
    du = K2 * (H2 * (1.5 * a * a - 0.5) - 3.0 * L2 * a * a)

    # step1 in phase (degree 3)
    dp += K3 * L3 * (7.5 * a * a - 1.5)
    du += K3 * (H3 * (2.5 * a * a * a - 1.5 * a) - L3 * (7.5 * a * a - 1.5) * a)

    # step1 out-of-phase (only radial)
    du += 3.0/4.0 * 0.0025 * K2 * np.sin(2.0*latp) * np.sin(2.0*pos[0]) * np.sin(pos[1]-lonp)
    du += 3.0/4.0 * 0.0022 * K2 * cosp**2 * cosl**2 * np.sin(2.0*(pos[1]-lonp))

    dr = dp * ep + du * eu
    return dr


def tide_solid(rsun, rmoon, pos, E, gmst, opt):
    """
    Solid Earth tide displacement.
    rsun, rmoon : position vectors of Sun and Moon
    pos         : [lat, lon] in radians
    E           : transformation matrix (flattened or 3x3)
    gmst        : Greenwich Mean Sidereal Time [rad]
    opt         : option flags
    returns dr (3,)
    """
    E = np.array(E).reshape(3,3)
    eu = np.array([E[0,2], E[1,2], E[2,2]])

    dr1 = tide_pl(eu, rsun, GMS, pos)
    dr2 = tide_pl(eu, rmoon, GMM, pos)

    # step2: frequency domain, only K1 radial
    sin2l = np.sin(2.0 * pos[0])
    du = -0.012 * sin2l * np.sin(gmst + pos[1])

    dr = dr1 + dr2 + du * eu

    # eliminate permanent deformation
    if opt & 8:
        sinl = np.sin(pos[0])
        du = 0.1196 * (1.5 * sinl**2 - 0.5)
        dn = 0.0247 * sin2l
        dr += du * eu + dn * E[:,1]

    return dr