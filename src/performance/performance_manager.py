import os

from src.data_mng.csv_data_mng import GnssRunStorageManager
from .error import *
from .plot_manager import *


class PerformanceEvaluation:
    @staticmethod
    def process(output_dir, data_manager: GnssRunStorageManager, config):

        # create post-processing folder, if not exists
        post_proc_dir = output_dir / 'post_proc'
        if not os.path.exists(post_proc_dir):
            print(f"Creating directory to store computed errors {post_proc_dir}...")
            os.makedirs(post_proc_dir)

        # create plot folder, if not exists
        plot_dir = output_dir / 'plot'
        if not os.path.exists(plot_dir):
            print(f"Creating directory to store plots {plot_dir}...")
            os.makedirs(plot_dir)

        PerformanceEvaluation._process_errors(post_proc_dir, config, data_manager)

        # 3- plots
        print("Plotting simulation data...")
        PerformanceEvaluation._plot()

    @staticmethod
    def _process_errors(post_proc_dir, config, data_manager):
        # fetch estimated data
        time = data_manager.get_data("time")
        position = data_manager.get_data("position")

        # compute Root Mean Square Error (RMSE) in ECEF and ENU frames
        print("Computing estimation error and root mean square error...")
        is_static = config["performance_evaluation"]["static"]
        if is_static:
            true_pos_dict = config["performance_evaluation"]["true_static_position"]
            true_pos = [true_pos_dict["x_ecef"], true_pos_dict["y_ecef"], true_pos_dict["z_ecef"]]
            error_ecef, error_enu = compute_error_static(position, true_pos, local="ENU")
        else:
            raise NotImplementedError(f"Dynamic error computation not yet implemented...")

        rms_ecef = compute_rms_error(error_ecef)
        rms_enu = compute_rms_error(error_enu)

        # 2- save computed errors to files
        print("Saving computed errors to files...")
        PerformanceEvaluation._write_outputs(post_proc_dir, time, error_ecef, error_enu, rms_ecef, rms_enu)

    @staticmethod
    def _write_outputs(output_dir, time, error_ecef, error_enu, rms_ecef, rms_enu):
        f_rms_stats = open(output_dir / "EstimationStats.txt", "w")

        error_ecef_ = np.concatenate((time.data, error_ecef), axis=1)
        error_enu_ = np.concatenate((time.data, error_enu), axis=1)

        # write headers
        np.savetxt(output_dir / "error_ecef.txt", error_ecef_, delimiter=',',
                   header='GPS_Week,GPS_TOW,x_error[m],y_error[m],z_error[m]', comments="")
        np.savetxt(output_dir / "error_enu.txt", error_enu_, delimiter=',',
                   header='GPS_Week,GPS_TOW,east_error[m],north_error[m],up_error[m]', comments="")

        # Overall Estimation Stats
        stats = "Root Mean Square Error:\n" \
                "\tECEF frame\n" \
                f"\t\tx = {rms_ecef['x']} [m]\n" \
                f"\t\ty = {rms_ecef['y']} [m]\n" \
                f"\t\tz = {rms_ecef['z']} [m]\n" \
                "\n\n\tENU frame\n" \
                f"\t\teast = {rms_enu['x']} [m]\n" \
                f"\t\tnorth = {rms_enu['y']} [m]\n" \
                f"\t\tup = {rms_enu['z']} [m]\n" \
                "\n\n\tOverall\n" \
                f"\t\tHorizontal Error = {rms_enu['2D']} [m]\n" \
                f"\t\t3D Position Error = {rms_enu['3D']} [m]\n"

        f_rms_stats.write(stats)
        print(stats)
        f_rms_stats.close()

    @staticmethod
    def _plot():
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
