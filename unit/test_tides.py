import numpy as np

from src.models.gnss_models.tidal_displacement import compute_displacement
from src.data_types.date import Epoch


def solid_testing():
    epoch = Epoch(2010, 1, 1, 0, 1, 0, scale="UTC")
    moon = np.array([-179996231.920342, -312468450.131567, -169288918.592160])
    sun = np.array([137859926952.015, 54228127881.4350, 23509422341.6960])
    rec = np.array([4075578.385, 931852.890, 4801570.154])

    disp_gmst82 = compute_displacement(epoch, sun, moon, rec, step3=False, gmst_model='IAU82')
    disp_gmst00 = compute_displacement(epoch, sun, moon, rec, step3=False, gmst_model='IAU00')
    disp_gmst06 = compute_displacement(epoch, sun, moon, rec, step3=False, gmst_model='IAU06')
    iers_disp = [0.7700420357108125891E-01, 0.6304056321824967613E-01, 0.5516568152597246810E-01]
    print("solid displacement (IAU 82)= ", disp_gmst82)
    print("solid displacement (IAU 00)= ", disp_gmst00)
    print("solid displacement (IAU 06)= ", disp_gmst06)
    print("IERS solid displacement = ", iers_disp)
    print("diff (IAU 82) = ", disp_gmst82 - iers_disp)
    print("diff (IAU 00) = ", disp_gmst00 - iers_disp)
    print("diff (IAU 06) = ", disp_gmst06 - iers_disp)


if __name__ == "__main__":
    solid_testing()
