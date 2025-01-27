""" Module with functions to write the State Space vectors to output files. """

from src import constants

__all__ = ["get_file_header", "export_to_file"]


def get_file_header(state_variable, state):
    """
    Get the header line of the provided `state_variable`.

    Args:
        state_variable(str): name of the state variable
        state(src.data_mng.gnss.state_space.GnssStateSpace): instance of the `GnssStateSpace`
            (to get fetch some metadata info about the state)
    Returns:
        str: header line string (CSV format)
    Raises:
        ValueError: the exception is raised if the provided `state_variable` is unavailable or not known.
    """
    epoch_system = state.epoch.scale
    if state_variable == "position":
        return f"Week_Number({epoch_system}),Time_of_Week[s],X_ECEF[m],Y_ECEF[m],Z_ECEF[m]," \
               f"cov_XX[m^2],cov_YY[m^2],cov_ZZ[m^2],cov_XY[m^2],cov_XZ[m^2],cov_YZ[m^2]"
    elif state_variable == "velocity":
        return f"Week_Number({epoch_system}),Time_of_Week[s],VX_ECEF[m/s],VY_ECEF[m/s],VZ_ECEF[m/s]," \
               f"cov_XX[m^2/s^2],cov_YY[m^2/s^2],cov_ZZ[m^2/s^2],cov_XY[m^2/s^2],cov_XZ[m^2/s^2],cov_YZ[m^2/s^2]"
    elif state_variable == "clock_bias":
        master = state.get_additional_info("clock_master")
        return f"Week_Number({epoch_system}),Time_of_Week[s],clock_bias(master={master})[s],cov[s^2]"
    elif state_variable == "clock_bias_rate":
        return f"Week_Number({epoch_system}),Time_of_Week[s],constellation,clock_bias_rate[s/s],cov[(s/s)^2]"
    elif state_variable == "iono":
        return f"Week_Number({epoch_system}),Time_of_Week[s],sat,iono[m],cov[m^2]"
    elif state_variable == "tropo_wet":
        return f"Week_Number({epoch_system}),Time_of_Week[s],tropo_wet[m],cov[m^2]"
    elif state_variable == "time":
        return f"Week_Number({epoch_system}),Time_of_Week[s],Epoch_timetag"
    elif state_variable == "pr_prefit_residuals":
        return f"Week_Number({epoch_system}),Time_of_Week[s],constellation,sat,data_type,residual[m]"
    elif state_variable == "pr_postfit_residuals":
        return f"Week_Number({epoch_system}),Time_of_Week[s],constellation,sat,data_type,residual[m]"
    elif state_variable == "pr_rate_prefit_residuals":
        return f"Week_Number({epoch_system}),Time_of_Week[s],constellation,sat,data_type,residual[m/s]"
    elif state_variable == "pr_rate_postfit_residuals":
        return f"Week_Number({epoch_system}),Time_of_Week[s],constellation,sat,data_type,residual[m/s]"
    elif state_variable == "satellite_azel":
        return f"Week_Number({epoch_system}),Time_of_Week[s],sat,azimuth[deg],elevation[deg]"
    elif state_variable == "dop_ecef":
        return f"Week_Number({epoch_system}),Time_of_Week[s],DOP_x,DOP_y,DOP_z,DOP_t,DOP_geometry,DOP_position"
    elif state_variable == "dop_local":
        return f"Week_Number({epoch_system}),Time_of_Week[s],DOP_East,DOP_North,DOP_Up,DOP_Horizontal"
    elif state_variable == "isb":
        master = state.get_additional_info("clock_master")
        slave = state.get_additional_info("clock_slave")
        return f"Week_Number({epoch_system}),Time_of_Week[s],ISB(master={master} slave={slave})[s],cov[s^2]"
    else:
        raise ValueError(f"Undefined header due to unknown state_variable '{state_variable}'")


def export_to_file(state_variable, state):
    """
    Export the provided state to a string (to be written in the output file, in CSV format).

    Args:
        state_variable(str): name of the state variable
        state(src.data_mng.gnss.state_space.GnssStateSpace): instance of the `GnssStateSpace`
    Returns:
        str: state variable converted to a string (along with covariance) in CSV format
    Raises:
        ValueError: the exception is raised if the provided `state_variable` is unavailable or not known.
    """
    if state_variable == "position":
        cov = state.cov_position
        cov_xx = cov[0, 0]
        cov_yy = cov[1, 1]
        cov_zz = cov[2, 2]
        cov_xy = cov[0, 1]
        cov_xz = cov[0, 2]
        cov_yz = cov[1, 2]
        return f"{state.position[0]},{state.position[1]},{state.position[2]},{cov_xx},{cov_yy}," \
               f"{cov_zz},{cov_xy},{cov_xz},{cov_yz}"

    if state_variable == "velocity":
        cov = state.cov_velocity
        cov_xx = cov[0, 0]
        cov_yy = cov[1, 1]
        cov_zz = cov[2, 2]
        cov_xy = cov[0, 1]
        cov_xz = cov[0, 2]
        cov_yz = cov[1, 2]
        return f"{state.velocity[0]},{state.velocity[1]},{state.velocity[2]},{cov_xx},{cov_yy}," \
               f"{cov_zz},{cov_xy},{cov_xz},{cov_yz}"

    if state_variable == "clock_bias":
        return f"{state.clock_bias/constants.SPEED_OF_LIGHT},{state.cov_clock_bias/constants.SPEED_OF_LIGHT**2}"

    if state_variable == "clock_bias_rate":
        data = []
        for constellation, clock_rate in state.clock_bias_rate.items():
            try:
                data.append(f"{constellation},{clock_rate/constants.SPEED_OF_LIGHT},"
                            f"{state.cov_clock_bias_rate[constellation]/constants.SPEED_OF_LIGHT**2}")
            except KeyError:
                pass
        return data

    if state_variable == "tropo_wet":
        return f"{state.tropo_wet},{state.cov_tropo_wet}"

    if state_variable == "isb":
        return f"{state.isb/constants.SPEED_OF_LIGHT},{state.cov_isb/constants.SPEED_OF_LIGHT**2}"

    if state_variable == "iono":
        cov = state.cov_iono

        data = []
        for sat, iono_data in state.iono.items():
            try:
                data.append(f"{sat},{iono_data},{cov[sat]}")
            except KeyError:
                pass
        return data

    elif state_variable in ["pr_prefit_residuals", "pr_postfit_residuals", "pr_rate_prefit_residuals",
                            "pr_rate_postfit_residuals"]:
        residuals = state.get_additional_info(state_variable)
        data = []

        for constellation, dict1 in residuals.items():
            for sat, dict2 in dict1.items():
                for observable, residual in dict2.items():
                    data.append(f"{constellation},{sat},{observable.data_type},{residual}")
        return data

    elif state_variable == "satellite_azel":
        geometry = state.get_additional_info("geometry")
        sat_list = geometry.get_satellites()
        data = []
        if len(sat_list) > 0:
            for sat in sat_list:
                el = geometry.get("el", sat) * constants.RAD2DEG
                az = geometry.get("az", sat) * constants.RAD2DEG
                data.append(f"{sat},{az},{el}")
        return data

    elif state_variable == "dop_ecef":
        dop = state._info["dop_ecef"]
        return f"{dop['x_ecef']},{dop['y_ecef']},{dop['z_ecef']},{dop['time']},{dop['geometry']},{dop['position']}"

    elif state_variable == "dop_local":
        dop = state._info["dop_local"]
        return f"{dop['east']},{dop['north']},{dop['up']},{dop['horizontal']}"

    elif state_variable == "time":
        return f"{str(state.epoch)}"

    else:
        raise ValueError(f"Undefined data due to unknown state_variable '{state_variable}'")
