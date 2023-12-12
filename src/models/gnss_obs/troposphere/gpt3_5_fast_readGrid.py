# troposphere-5 Tropospheric Model
# (c) Department of Geodesy and Geoinformation, Vienna University of Technology, 2017
# Source: https://vmf.geo.tuwien.ac.at/codes/Python_Tools_Adavi/
# This troposphere model has been adapted from the reference above
# Disclaimer: All credits for the original code go to the Department of Geodesy and Geoinformation,
#              Vienna University of Technology. This code has been adapted for use in this project

import numpy as np
import pandas as pd


def gpt3_5_fast_readGrid(filename='gpt3_5.grd'):
    # =============================================================================
    #     gpt3_5_fast_readGrid.py
    #     This routine reads the grid 'gpt3_5.grd' and overgives the respective
    #     parameters in a cell. It must be placed ahead of the for loop, which runs
    #    through all observations in order to save a huge time amount. The related
    #    function "gpt3_5_fast.py" must then replace the default
    #   "gpt3_5.m" at the respective location in the text.

    #    File created by Daniel Landskron, 2017-04-12

    #    File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
    # =============================================================================

    # read gridfile

    file = open(filename)
    C = pd.read_csv(file, skiprows=[i for i in range(0, 1)], sep='\s+', skipfooter=0, skip_blank_lines=True,
                    header=None, index_col=False)
    file.close()

    p_grid = C.loc[:, np.arange(2, 7)]

    T_grid = C.loc[:, np.arange(7, 12)]

    Q_grid = C.loc[:, np.arange(12, 17)] / 1000

    dT_grid = C.loc[:, np.arange(17, 22)] / 1000

    u_grid = C.loc[:, 22]

    Hs_grid = C.loc[:, 23]

    ah_grid = C.loc[:, np.arange(24, 29)] / 1000

    aw_grid = C.loc[:, np.arange(29, 34)] / 1000

    la_grid = C.loc[:, np.arange(34, 39)]

    Tm_grid = C.loc[:, np.arange(39, 44)]

    Gn_h_grid = C.loc[:, np.arange(44, 49)] / 100000

    Ge_h_grid = C.loc[:, np.arange(49, 54)] / 100000

    Gn_w_grid = C.loc[:, np.arange(54, 59)] / 100000

    Ge_w_grid = C.loc[:, np.arange(59, 64)] / 100000

    # combine all data to on cell grid
    gpt3_grid = {}
    gpt3_grid[0] = p_grid
    gpt3_grid[1] = T_grid
    gpt3_grid[2] = Q_grid
    gpt3_grid[3] = dT_grid
    gpt3_grid[4] = u_grid
    gpt3_grid[5] = Hs_grid
    gpt3_grid[6] = ah_grid
    gpt3_grid[7] = aw_grid
    gpt3_grid[8] = la_grid
    gpt3_grid[9] = Tm_grid
    gpt3_grid[10] = Gn_h_grid
    gpt3_grid[11] = Ge_h_grid
    gpt3_grid[12] = Gn_w_grid
    gpt3_grid[13] = Ge_w_grid

    return gpt3_grid