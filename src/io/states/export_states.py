from src import constants
from src.algorithms.gnss.estimators.state_space import GnssStateSpace


def get_file_header(exportable, state):
    epoch_system = state.epoch.scale
    if exportable == "position":
        return f"Week_Number({epoch_system}),Time_of_Week[s],X_ECEF[m],Y_ECEF[m],Z_ECEF[m]," \
               f"cov_XX[m^2],cov_YY[m^2],cov_ZZ[m^2],cov_XY[m^2],cov_XZ[m^2],cov_YZ[m^2]"
    elif exportable == "clock_bias":
        master = state.get_additional_info("clock_master")
        return f"Week_Number({epoch_system}),Time_of_Week[s],clock_bias(master={master})[s],cov[s^2]"
    elif exportable == "iono":
        return f"Week_Number({epoch_system}),Time_of_Week[s],sat,iono[m],cov[m^2]"
    elif exportable == "tropo_wet":
        return f"Week_Number({epoch_system}),Time_of_Week[s],tropo_wet[m],cov[m^2]"
    elif exportable == "time":
        return f"Week_Number({epoch_system}),Time_of_Week[s],Epoch_timetag"
    elif exportable == "prefit_residuals":
        return f"Week_Number({epoch_system}),Time_of_Week[s],constellation,sat,data_type,residual[m]"
    elif exportable == "postfit_residuals":
        return f"Week_Number({epoch_system}),Time_of_Week[s],constellation,sat,data_type,residual[m]"
    elif exportable == "satellite_azel":
        return f"Week_Number({epoch_system}),Time_of_Week[s],sat,azimuth[deg],elevation[deg]"
    elif exportable == "dop_ecef":
        return f"Week_Number({epoch_system}),Time_of_Week[s],DOP_x,DOP_y,DOP_z,DOP_t,DOP_geometry,DOP_position"
    elif exportable == "dop_local":
        return f"Week_Number({epoch_system}),Time_of_Week[s],DOP_East,DOP_North,DOP_Up,DOP_Horizontal"
    elif exportable == "isb":
        master = state.get_additional_info("clock_master")
        slave = state.get_additional_info("clock_slave")
        return f"Week_Number({epoch_system}),Time_of_Week[s],ISB(master={master},slave={slave})[s],cov[s^2]"
    else:
        raise ValueError(f"Undefined header due to unknown exportable '{exportable}'")


def export_to_file(gnss_state: GnssStateSpace, exportable):
    if exportable == "position":
        cov = gnss_state.cov_position
        cov_xx = cov[0, 0]
        cov_yy = cov[1, 1]
        cov_zz = cov[2, 2]
        cov_xy = cov[0, 1]
        cov_xz = cov[0, 2]
        cov_yz = cov[1, 2]
        return f"{gnss_state.position[0]},{gnss_state.position[1]},{gnss_state.position[2]},{cov_xx},{cov_yy}," \
               f"{cov_zz},{cov_xy},{cov_xz},{cov_yz}"

    if exportable == "clock_bias":
        return f"{gnss_state.clock_bias},{gnss_state.cov_clock_bias}"

    if exportable == "tropo_wet":
        return f"{gnss_state.tropo_wet},{gnss_state.cov_tropo_wet}"

    if exportable == "isb":
        return f"{gnss_state.isb},{gnss_state.cov_isb}"

    if exportable == "iono":
        cov = gnss_state.cov_iono

        data = []
        for sat, iono_data in gnss_state.iono.items():
            data.append(f"{sat},{iono_data},{cov[sat]}")
        return data

    elif exportable == "prefit_residuals" or exportable == "postfit_residuals":
        residuals = gnss_state.get_additional_info(exportable)
        data = []

        for constellation, dict1 in residuals.items():
            for sat, dict2 in dict1.items():
                for observable, residual in dict2.items():
                    data.append(f"{constellation},{sat},{observable.data_type},{residual}")
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

    elif exportable == "time":
        return f"{str(gnss_state.epoch)}"

    else:
        raise ValueError(f"Undefined data due to unknown exportable '{exportable}'")
