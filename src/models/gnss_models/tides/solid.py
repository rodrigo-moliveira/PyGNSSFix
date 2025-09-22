# see https://github.com/tomojitakasu/RTKLIB/blob/71db0ffa0d9735697c6adfd06fdf766d0e5ce807/src/ppp.c#L184
# see also file:///C:/Users/rooo/Downloads/iersconventions_v1_3_0/tn36.pdf
# see also https://www.rtklib.com/prog/manual_2.4.2.pdf

from math import sin
import numpy as np

from src.constants import MASS_RATIO_SUN, MASS_RATIO_MOON

R_EARTH = 6378136.6  # specific value for the IERS 2010 convention
H2 = 0.6078
L2 = 0.0847
def solid_disp(r_sun, r_moon, r_rec, geoc_lat):
    disp = 0.0

    r_rec_unit = r_rec / np.linalg.norm(r_rec)
    h2 = H2 - 0.0006 * ((3 * sin(geoc_lat)**2 - 1) / 2)
    l2 = L2 + 0.0002 * ((3 * sin(geoc_lat)**2 - 1) / 2)

    # Sun contribution (order 2 only)
    r_sun_norm = np.linalg.norm(r_sun)
    r_sun_unit = r_sun / r_sun_norm
    dot_prod = np.dot(r_sun_unit, r_rec_unit)
    h2_term = h2 * r_rec_unit * (1.5 * dot_prod**2 - 0.5)
    l2_term = 3 * l2 * dot_prod * (r_sun_unit - dot_prod * r_rec_unit)
    disp_sun = MASS_RATIO_SUN * R_EARTH**4 / r_sun_norm**3 * (h2_term + l2_term)

    # Moon contribution (order 2 and 3)
    r_moon_norm = np.linalg.norm(r_moon)
    r_moon_unit = r_moon / r_moon_norm
    dot_prod = np.dot(r_moon_unit, r_rec_unit)
    h2_term = h2 * r_rec_unit * (1.5 * dot_prod ** 2 - 0.5)
    l2_term = 3 * l2 * dot_prod * (r_moon_unit - dot_prod * r_rec_unit)
    disp_sun = MASS_RATIO_MOON * R_EARTH ** 4 / r_moon_norm ** 3 * (h2_term + l2_term)
