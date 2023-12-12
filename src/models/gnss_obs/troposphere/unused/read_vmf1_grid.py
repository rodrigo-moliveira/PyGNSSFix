"""# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import pandas as pd
import os
from scipy.interpolate import griddata


def read_vmf1_grid(mjd=None, ell=None):
    # =============================================================================
    #     Function read_vmf1_grid searches for grid file close in time and computes
    #     mapping coefficients ah, aw for selected station

    #     Input parameters:
    #     -------------------
    #     mjd   ... Modified Julian Date
    #     ell   ... ellipsoidal station coordinates lat [degrees], lon [degrees] ,hell [m] 'Should be in list form e.g 'ell=[42.72,56.13,570]''
    #     File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
    # =============================================================================

    # Closest 6 hour period
    mjd = np.array(mjd)
    mjd_6h = np.round(mjd * 4) / 4

    # List of relevant 6 hour intervals
    mjd_6h_list = np.unique(mjd_6h)

    # Read vmf1 grid of period before and after current date
    # (needed for temporal interpolation)

    if len(mjd_6h_list) == 1:
        mjd_6h_list = np.array([mjd_6h_list - 0.25, mjd_6h_list, mjd_6h_list + 0.25])

    # Initialise variable
    ah_ell = np.zeros((len(mjd_6h_list), 1))
    aw_ell = np.zeros((len(mjd_6h_list), 1))
    for i in range(0, len(mjd_6h_list)):
        # Convert mjd to date

        yr, mn, dy = jd_to_date(np.rint(np.amin(mjd_6h_list[i])) + 2400000.5)
        # Hour of the day
        hh = (mjd_6h_list[i] - np.rint(mjd_6h_list[i])) * 24

        # Open grid file with mapping coefficients
        pwd = os.getcwd()
        directory_vmf = pwd + "\\Mapping_Fcn\\vmf1\\VMFG_" + str(int(yr)) + str(int(mn)).zfill(2) + str(int(dy)).zfill(
            2) + '.H' + str(int(hh)).zfill(2)
        file = open(directory_vmf)

        ###############################
        pos = 0
        line0 = '! Comment:'
        cur_line = file.readline()
        Columns = ['lat', 'lon', 'ah', 'aw', 'zhd', 'zwd']
        while not cur_line.startswith(line0):
            pos = pos + 1
            cur_line = file.readline()

        file.close()
        file = open(directory_vmf)
        dat = pd.read_csv(file, skiprows=[i for i in range(0, pos)], sep='\s+', skipfooter=0, skip_blank_lines=True,
                          header=0, index_col=False, names=Columns)
        file.close()

        # Restore grid data
        lat = np.array(dat.iloc[:, 0])

        lon = np.array(dat.iloc[:, 1])

        ah = np.array(dat.iloc[:, 2])

        aw = np.array(dat.iloc[:, 3])
        # Interpolate data to station coordinates
        ah_ell[i] = griddata((lat, lon), ah, (ell[0], ell[1]))

        aw_ell[i] = griddata((lat, lon), aw, (ell[0], ell[1]))

    # Interpolate and store final coefficients to return
    ah_vec = np.interp(np.transpose(mjd), mjd_6h_list, ah_ell)
    aw_vec = np.interp(np.transpose(mjd), mjd_6h_list, aw_ell)

    return ah_vec, aw_vec"""