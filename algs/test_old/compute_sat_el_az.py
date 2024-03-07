import numpy as np

from src.models.frames import cartesian2geodetic
from src.models.frames.frames import ecef2enu, enu2azel

rec_pos = [3929894.266,
        309272.3398,
        4997323.8144]

pos_sat_1 = np.array([17221.403071, -13865.306516,  19674.221505]) * 1000
pos_sat_2 = np.array([18760.239983, -19710.173016 , 11655.159833]) * 1000


#pos_sat_1 = np.array([18361.847030, -14036.730292,  18484.970725]) * 1000
#pos_sat_2 = np.array([-23969.480029, -16642.417861,  -4964.261518]) * 1000

#pos_sat_1 = np.array([18901.243130, -14139.960018 , 17851.553740]) * 1000
#pos_sat_2 = np.array([-24916.772892, -15184.084759 , -4974.259012]) * 1000


lat, long, h = cartesian2geodetic(rec_pos[0], rec_pos[1], rec_pos[2])
enu_coord = ecef2enu(pos_sat_1[0], pos_sat_1[1], pos_sat_1[2], lat, long, h)
az, el = enu2azel(*enu_coord)
print("sat E10", az, el)

enu_coord = ecef2enu(pos_sat_2[0], pos_sat_2[1], pos_sat_2[2], lat, long, h)
az, el = enu2azel(*enu_coord)
print("sat E04", az, el)