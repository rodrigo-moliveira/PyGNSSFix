from src.data_mng.csv_data_mng import GnssRunStorageManager
from src.data_mng.timeseries import TimeSeries


def rotate_dops(dop_ecef):
    return None

class PerformanceEvaluation:
    @staticmethod
    def process(data_manager: GnssRunStorageManager):


        # compute RMS
        RMS_ECEF = TimeSeries()
        RMS_ENU = TimeSeries()
        RMS_ECEF.compute_errors(receiver_pos, true_position, "static", "ECEF")
        RMS_ENU.compute_errors(receiver_pos, true_position, "static", "ENU")

        # 2- save to files
        GNSSQualityManager._write_trace_data(trace_path, sat_info, prefit_residuals, postfit_residuals, estimated_iono)
        GNSSQualityManager._write_outputs(output_path, receiver_pos, receiver_bias, DOPs, RMS_ECEF, RMS_ENU)

        log.info("########## End of module 'PVT Quality Check' ... ###########\n")

        # 3- plots
        if plot:
            plot_1D_TimeSeries(receiver_bias, x_label="Time", y_label="clock [s]", title="Receiver Clock Bias")

            plot_3D_trajectory(list(receiver_pos.values()), true_position=true_position.copy(form="cartesian"),
                               x_label="X ECEF [m]", y_label="Y ECEF [m]", z_label="Z ECEF [m]",
                               title="Estimated VS True position")

            GNSSQualityManager.plot_errors(RMS_ECEF, "RMS Estimation Error - ECEF", "X", "Y", "Z", norm=True)
            GNSSQualityManager.plot_errors(RMS_ENU, "RMS Estimation Error - ENU", "East", "North", "Up", norm=True)
            GNSSQualityManager.plot_dops(DOPs)

            plot_skyplot(sat_info)

            plot_satellite_availability(sat_info, x_label="Time", y_label="Number of Sats",
                                        title="Satellite Availability")

            _, x, y, _ = RMS_ENU.export2time_data()
            plot_2D_trajectory(x, y, x_label="East [m]", y_label="North [m]", title="Horizontal Position Error")

            if not estimated_iono.is_empty():
                GNSSQualityManager.plot_iono(estimated_iono)

            show_all()

    @staticmethod
    def _write_trace_data(trace_path, sat_info, prefit_residuals, postfit_residuals, estimated_iono):
        f = open(trace_path + "/Satellite_Info.txt", "w")
        f.write(str(sat_info))
        f.close()

        f = open(trace_path + "/Prefit_Residuals.txt", "w")
        f.write(str(prefit_residuals))
        f.close()

        f = open(trace_path + "/Postfit_Residuals.txt", "w")
        f.write(str(postfit_residuals))
        f.close()

        if not estimated_iono.is_empty():
            f = open(trace_path + "/EstimatedIono.txt", "w")
            f.write(str(estimated_iono))
            f.close()

    @staticmethod
    def _write_outputs(output_path, receiver_pos, receiver_bias, DOPs, RMS_ECEF, RMS_ENU):
        f_PT = open(output_path + "/PositionTime.txt", "w")
        f_DOP = open(output_path + "/DOPs.txt", "w")
        f_DOP_ECEF = open(output_path + "/DOPs_ECEF.txt", "w")
        f_DOP_ENU = open(output_path + "/DOPs_ENU.txt", "w")
        f_RMS_ECEF = open(output_path + "/RMS_ECEF.txt", "w")
        f_RMS_ENU = open(output_path + "/RMS_ENU.txt", "w")
        f_RMS_stats = open(output_path + "/Stats.txt", "w")

        # write headers
        f_PT.write(f"Time,X[m],Y[m],Z[m],lat[deg],long[deg],height[m],clock_dt[s]\n")
        f_DOP.write(f"Time,geometry_DOP,position_DOP,time_DOP,horizontal_DOP\n")
        f_DOP_ECEF.write(f"Time,x_DOP,y_DOP,z_DOP\n")
        f_DOP_ENU.write(f"Time,east_DOP,north_DOP,up_DOP\n")
        f_RMS_ECEF.write(f"Time,x_RMS[m],y_RMS[m],z_RMS[m]\n")
        f_RMS_ENU.write(f"Time,east_RMS[m],north_RMS[m],up_RMS[m]\n")

        for epoch, position in receiver_pos.items():
            # get data for epoch
            bias = receiver_bias.get_data_for_epoch(epoch)
            dop = DOPs.get_data_for_epoch(epoch)
            rms_ecef = RMS_ECEF.get_data_for_epoch(epoch)
            rms_enu = RMS_ENU.get_data_for_epoch(epoch)

            position_geodetic = position.copy()
            position_geodetic.form = "geodetic"

            # save to files
            f_PT.write(f"{epoch.to_time_stamp()},{position[0]},{position[1]},{position[2]},"
                       f"{position_geodetic[0] * Constant.RAD2DEG},{position_geodetic[1] * Constant.RAD2DEG},"
                       f"{position_geodetic[2]},{bias}\n")
            f_DOP.write(f"{epoch.to_time_stamp()},{dop.geometry},{dop.position},{dop.time},{dop.horizontal}\n")
            f_DOP_ECEF.write(f"{epoch.to_time_stamp()},{dop.x_ecef},{dop.y_ecef},{dop.z_ecef}\n")
            f_DOP_ENU.write(f"{epoch.to_time_stamp()},{dop.east},{dop.north},{dop.up}\n")
            f_RMS_ECEF.write(f"{epoch.to_time_stamp()},{rms_ecef.x},{rms_ecef.y},{rms_ecef.z}\n")
            f_RMS_ENU.write(f"{epoch.to_time_stamp()},{rms_enu.x},{rms_enu.y},{rms_enu.z}\n")

        # Overall Estimation Stats
        stats = "Root Mean Square Error:\n" \
                "\tECEF frame\n" \
                f"\t\tx = {RMS_ECEF.stats['x']} [m]\n" \
                f"\t\ty = {RMS_ECEF.stats['y']} [m]\n" \
                f"\t\tz = {RMS_ECEF.stats['z']} [m]\n" \
                "\n\n\tENU frame\n" \
                f"\t\teast = {RMS_ENU.stats['x']} [m]\n" \
                f"\t\tnorth = {RMS_ENU.stats['y']} [m]\n" \
                f"\t\tup = {RMS_ENU.stats['z']} [m]\n" \
                "\n\n\tOverall\n" \
                f"\t\tHorizontal Error = {RMS_ENU.stats['2D']} [m]\n" \
                f"\t\t3D Position Error = {RMS_ENU.stats['3D']} [m]\n"

        f_RMS_stats.write(stats)
        print(stats)

        f_PT.close()
        f_DOP.close()
        f_DOP_ECEF.close()
        f_DOP_ENU.close()
        f_RMS_ECEF.close()
        f_RMS_ENU.close()
        f_RMS_stats.close()

    @staticmethod
    def plot_outputs():
        pass

    @staticmethod
    def plot_errors(rms, title, x_lab, y_lab, z_lab, norm=False):
        ax = None

        if norm:
            t, x, y, z, _norm = rms.export2time_data(norm=True)
            if norm:
                ax = plot_1D(t, _norm, label="3D")
        else:
            t, x, y, z = rms.export2time_data(norm=False)

        ax = plot_1D(t, x, ax=ax, label=x_lab)
        plot_1D(t, y, ax=ax, label=y_lab)
        plot_1D(t, z, ax=ax, x_label="Time", y_label="Position error [m]", label=z_lab,
                title=title, set_legend=True)

    @staticmethod
    def plot_dops(rms):
        t, _, geometry, position, time, horizontal, x, y, z, east, north, up = rms.export2time_data()

        ax = plot_1D(t, geometry, label="geometry")
        plot_1D(t, position, ax=ax, label="position")
        plot_1D(t, horizontal, ax=ax, x_label="Time", y_label="DOPs [m]", label="horizontal",
                title="Dilution of Precision", set_legend=True)

        ax = plot_1D(t, x, label="x")
        plot_1D(t, y, ax=ax, label="y")
        plot_1D(t, z, ax=ax, x_label="Time", y_label="DOPs [m]", label="z",
                title="Dilution of Precision in ECEF", set_legend=True)

        ax = plot_1D(t, east, label="east")
        plot_1D(t, north, ax=ax, label="north")
        plot_1D(t, up, ax=ax, x_label="Time", y_label="DOPs [m]", label="up",
                title="Dilution of Precision in ENU", set_legend=True)

    @staticmethod
    def plot_iono(iono_tm):
        ax = None
        sat_data = {}
        epochs = iono_tm.get_all_epochs()

        for epoch in epochs:
            data = iono_tm.get_data_for_epoch(epoch)
            for sat, iono in data:
                if sat not in sat_data:
                    sat_data[sat] = []
                sat_data[sat].append([epoch, iono])

        for sat, data in sat_data.items():
            epoch = [x[0] for x in data]
            iono = [x[1] for x in data]

            ax = plot_1D(epoch, iono, ax=ax, label=sat, x_label="Time", y_label="Ionosphere [m]",
                         title="Estimated Ionosphere", set_legend=True)
