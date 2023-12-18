import numpy as np

from src.models.gnss_obs.troposphere.tropo_gpt import GPTModel
from tropo_saastamoinen import tropo_saastamoinen
from src.data_types.date.date import Epoch
from gmf_GMF_function import gmf

print("### TROPO TEST ###")

height = 180
lat = 0.9
long = 1.0
doy = 1
mjd = 59045 # 1st of january 2023
el = 60 * np.pi/180.0
zenith_angle = np.pi/2.0 - el

# SAASTAMOINEN MODEL
t1 = tropo_saastamoinen(height, lat, doy, el)

# GPT3 Model with VMF1 mapping functions
model = GPTModel()
t2 = model.compute(mjd, lat, long, height, zenith_angle, it=0)


print(t1, t2)

# VMF1 Map
# ah_vec, aw_vec = read_vmf1_grid(mjd=None, ell=None)
# vmf1_map = vmf1_ht(ah=None, aw=None, dmjd=None, dlat=None, ht=None, zd=None)

# GMF Map
#gmf_map = gmf(dmjd=mjd, dlat=lat, dlon=long, dhgt=height, zd=zenith_angle)


# TODO: see what to do with the east/north gradients..