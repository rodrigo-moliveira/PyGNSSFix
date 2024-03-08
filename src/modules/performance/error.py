import numpy as np

from src.data_mng.csv_data import CSVData
from src.models.frames import cartesian2geodetic, latlon2dcm_e_enu, latlon2dcm_e_ned


def error_vector(x, true_pos, rot_matrix=np.eye(3)):
    return rot_matrix @ (true_pos - x)


def compute_error_static(est_vec: CSVData, true_vec, true_pos=None, local="ENU"):
    if true_pos is None:
        true_pos = true_vec
    [lat, lon, _] = cartesian2geodetic(*true_pos)

    if local == "ENU":
        rot_e_local = latlon2dcm_e_enu(lat, lon)  # rotation matrix from ECEF to ENU
    elif local == "NED":
        rot_e_local = latlon2dcm_e_ned(lat, lon)  # rotation matrix from ECEF to NED
    else:
        raise ValueError(f"argument local must either be 'ENU' or 'NED'")

    error_matrix_ecef = np.apply_along_axis(error_vector, 1, est_vec.data, true_vec)
    error_matrix_local = np.apply_along_axis(error_vector, 1, est_vec.data, true_vec, rot_e_local)

    return error_matrix_ecef, error_matrix_local


def compute_rms_error(error_matrix):
    accum_x = accum_y = accum_z = accum_2d = accum_3d = 0

    for error in error_matrix:
        accum_x += error[0] * error[0]
        accum_y += error[1] * error[1]
        accum_z += error[2] * error[2]

        accum_2d += error[0] * error[0] + error[1] * error[1]
        accum_3d += error[0] * error[0] + error[1] * error[1] + error[2] * error[2]

    # 1D RMS stats
    rms_x = np.sqrt(1 / len(error_matrix) * accum_x)
    rms_y = np.sqrt(1 / len(error_matrix) * accum_y)
    rms_z = np.sqrt(1 / len(error_matrix) * accum_z)

    # 2D RMS stats
    rms_2d = np.sqrt(1 / len(error_matrix) * accum_2d)

    # 3D RMS stats
    rms_3d = np.sqrt(1 / len(error_matrix) * accum_3d)

    return {"x": rms_x,
            "y": rms_y,
            "z": rms_z,
            "2D": rms_2d,
            "3D": rms_3d
            }
