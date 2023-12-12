import numpy as np
from tropo_saastamoinen import tropo_saastamoinen
from gmf_GMF_function import gmf
from gpt3_5_fast_readGrid import gpt3_5_fast_readGrid
from gpt3_5_fast import  gpt3_5_fast

print("### TROPO TEST ###")

height = 0
lat = 0.9
long = 1.0
doy = 1
mjd = 59045 # 1st of january 2023
el = 60 * np.pi/180.0
zenith_angle = np.pi/2.0 - el

# SAASTAMOINEN MODEL
t1 = tropo_saastamoinen(height, lat, doy, el)

# GPT3 Model with VMF1 mapping functions
grid = gpt3_5_fast_readGrid()  # GPT grid file
gpt3_data = gpt3_5_fast(mjd=mjd, lat=[lat,], lon=[long,], h_ell=[height,], it=1, grid=grid)  # GPT data
# gpt3_data = [p, T, dT, Tm, e, ah, aw, la, undu, Gn_h, Ge_h, Gn_w, Ge_w]

# VMF1 Map
# ah_vec, aw_vec = read_vmf1_grid(mjd=None, ell=None)
# vmf1_map = vmf1_ht(ah=None, aw=None, dmjd=None, dlat=None, ht=None, zd=None)

# GMF Map
gmf_map = gmf(dmjd=mjd, dlat=lat, dlon=long, dhgt=height, zd=zenith_angle)


# k2_ and k3 empirically determined refractivity constants
# Rd specific gas 287.0464; gm 9.80665
# e: water vapor pressure, Tm: mean temperature weighted with water vapor pressure, lambda: water vator decrease factor (from GPT3 model)

tropo = map_hydro * ZHD + map_wet * ZWD

# TODO: see what to do with the east/north gradients..