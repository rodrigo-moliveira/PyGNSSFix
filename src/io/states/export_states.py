import numpy as np

from src import constants
from src.algorithms.gnss.estimators.state_space import GnssStateSpace


def get_file_header(exportable):
    if exportable == "position":
        return "GPS_Week,GPS_TOW,X_ECEF[m],Y_ECEF[m],Z_ECEF[m],cov_XX[m^2],cov_YY[m^2],cov_ZZ[m^2],cov_XY[m^2]," \
               "cov_XZ[m^2],cov_YZ[m^2]"
    elif exportable == "clock_bias":
        return "GPS_Week,GPS_TOW,clock_bias[s],cov[s^2]"
    elif exportable == "iono":
        return "GPS_Week,GPS_TOW,sat,iono[m],cov[m^2]"
    elif exportable == "epoch":
        return "GPS_Week,GPS_TOW,Epoch_timetag"
    elif exportable == "prefit_residuals":
        return "GPS_Week,GPS_TOW,sat,prefit_residuals_i[m^2]"
    elif exportable == "postfit_residuals":
        return "GPS_Week,GPS_TOW,sat,postfit_residuals_i[m^2]"
    elif exportable == "satellite_azel":
        return "GPS_Week,GPS_TOW,sat,azimuth[deg],elevation[deg]"
    elif exportable == "dop_ecef":
        return "GPS_Week,GPS_TOW,DOP_x,DOP_y,DOP_z,DOP_t,DOP_geometry,DOP_position"
    elif exportable == "dop_local":
        return "GPS_Week,GPS_TOW,DOP_East,DOP_North,DOP_Up,DOP_Horizontal"
    else:
        raise ValueError(f"Undefined header due to unknown exportable '{exportable}'")


def export_to_file(gnss_state: GnssStateSpace, exportable):
    if exportable == "position":
        cov = gnss_state._info.get("cov", None)
        if cov is not None:
            cov_xx = cov[0, 0]
            cov_yy = cov[1, 1]
            cov_zz = cov[2, 2]
            cov_xy = cov[0, 1]
            cov_xz = cov[0, 2]
            cov_yz = cov[1, 2]
            return f"{gnss_state.position[0]},{gnss_state.position[1]},{gnss_state.position[2]},{cov_xx},{cov_yy}," \
                   f"{cov_zz},{cov_xy},{cov_xz},{cov_yz}"
        else:
            return f"{gnss_state.position[0]},{gnss_state.position[1]},{gnss_state.position[2]}"

    if exportable == "clock_bias":
        cov = gnss_state._info.get("cov", None)
        if cov is not None:
            cov_t = cov[3, 3] / (constants.SPEED_OF_LIGHT ** 2)  # in seconds^2
            return f"{gnss_state.clock_bias},{cov_t}"
        else:
            return f"{gnss_state.clock_bias}"

    if exportable == "iono":
        cov = gnss_state._info.get("cov", None)
        sat_list = gnss_state.get_additional_info("geometry").get_satellites()
        data = []
        iono_list = gnss_state.iono

        for sat, iono, cov_iono in zip(sat_list, iono_list, np.diag(cov)[4:]):
            data.append(f"{sat},{iono},{cov_iono}")

        return data

    elif exportable == "prefit_residuals" or exportable == "postfit_residuals":
        residuals = gnss_state._info[exportable]
        sat_list = gnss_state.get_additional_info("geometry").get_satellites()
        data = []
        if len(sat_list) > 0:
            for sat, res in zip(sat_list, residuals):
                data.append(f"{sat},{res}")
        return data

    elif exportable == "satellite_azel":
        geometry = gnss_state.get_additional_info("geometry")
        sat_list = geometry.get_satellites()
        data = []
        if len(sat_list) > 0:
            for sat in sat_list:
                el = geometry.get("el", sat) * constants.RAD2DEG
                az = geometry.get("az", sat) * constants.RAD2DEG
                data.append(f"{sat},{az},{el}")
        return data

    elif exportable == "dop_ecef":
        dop = gnss_state._info["dop_ecef"]
        return f"{dop['x_ecef']},{dop['y_ecef']},{dop['z_ecef']},{dop['time']},{dop['geometry']},{dop['position']}"

    elif exportable == "dop_local":
        dop = gnss_state._info["dop_local"]
        return f"{dop['east']},{dop['north']},{dop['up']},{dop['horizontal']}"

    elif exportable == "epoch":
        return f"{str(gnss_state.epoch)}"

    else:
        raise ValueError(f"Undefined data due to unknown exportable '{exportable}'")
