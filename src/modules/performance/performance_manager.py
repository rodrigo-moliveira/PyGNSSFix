import os

from .error import *
from .plots import *
import pandas as pd


class PerformanceManager:

    def __init__(self, data_manager, config):
        # TODO: add logger..
        self.data_manager = data_manager
        self.config = config

        self.true_pos = None
        self.error_ecef = None
        self.error_enu = None

    def process(self, output_dir):

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

        print("Running post-processing scripts ")
        print("1 - Processing estimation errors...")
        self._process_errors(post_proc_dir)

        print("2 - Running residual analysis")
        # TODO...

        # plots
        print("Plotting simulation data...")
        self._plot()

    def _process_errors(self, post_proc_dir):
        # fetch estimated data
        time = self.data_manager.get_data("time")
        position = self.data_manager.get_data("position")

        # compute Root Mean Square Error (RMSE) in ECEF and ENU frames
        print("Computing estimation error and root mean square error...")
        is_static = self.config.get("performance_evaluation", "static")
        if is_static:
            self.true_pos = self.config.get("performance_evaluation", "true_static_position")
            self.error_ecef, self.error_enu = compute_error_static(position, self.true_pos, local="ENU")
        else:
            raise NotImplementedError(f"Dynamic error computation not yet implemented...")

        rms_ecef = compute_rms_error(self.error_ecef)
        rms_enu = compute_rms_error(self.error_enu)

        # 2- save computed errors to files
        print("Saving computed errors to files...")
        self._write_outputs(post_proc_dir, time, rms_ecef, rms_enu)

    def _write_outputs(self, output_dir, time, rms_ecef, rms_enu):
        f_rms_stats = open(output_dir / "EstimationStats.txt", "w")

        # Create DataFrame for error matrix
        ecef_df = pd.DataFrame(self.error_ecef, columns=['x_error[m]', 'y_error[m]', 'z_error[m]'])
        error_ecef_df = pd.concat([time.data, ecef_df], axis=1)
        enu_df = pd.DataFrame(self.error_enu, columns=['east_error[m]', 'north_error[m]', 'up_error[m]'])
        error_enu_df = pd.concat([time.data, enu_df], axis=1)

        # Save the error dataframes
        error_ecef_df.to_csv(output_dir / "error_ecef.txt", index=False)
        error_enu_df.to_csv(output_dir / "error_enu.txt", index=False)

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

    def _plot(self):
        # fetch estimated data
        time = format_time_for_plot(self.data_manager.get_data("time").data)
        position = self.data_manager.get_data("position")

        dop_ecef = self.data_manager.get_data("dop_ecef")
        dop_enu = self.data_manager.get_data("dop_local")
        azel = self.data_manager.get_data("satellite_azel")
        residuals = self.data_manager.get_data("postfit_residuals")

        self._plot_clock_bias(time)

        # plot_3D_trajectory(position.data, true_position=self.true_pos,
        #                    x_label="X ECEF [m]", y_label="Y ECEF [m]", z_label="Z ECEF [m]",
        #                    title="Estimated VS True position")
        #
        # self._plot_errors(time, self.error_ecef, "ECEF", "X", "Y", "Z")
        # self._plot_errors(time, self.error_enu, "ENU", "East", "North", "Up")
        # self._plot_dops(time, dop_ecef, dop_enu)
        #
        # plot_satellite_availability(residuals.data, x_label="Time", y_label="Number of Sats",
        #                             title="Satellite Availability")
        # plot_skyplot(azel.data)
        #
        # plot_2D_trajectory(self.error_enu,
        #                    x_label="East [m]", y_label="North [m]", title="Horizontal Position Error")

        """
        if not estimated_iono.is_empty():
            GNSSQualityManager.plot_iono(estimated_iono)
        """
        show_all()

    def _plot_clock_bias(self, time):
        clock_bias = self.data_manager.get_data("clock_bias")
        data_matrix = clock_bias.to_data_array()

        clock_series = data_matrix[:, 0]
        sigma_series = np.sqrt(data_matrix[:, 1])

        ax = plot_1D(time, clock_series, x_label="Time", label='clock bias', color='blue')
        ax = plot_1D(time, clock_series+sigma_series, ax=ax, x_label="Time", label="+/- sigma [s]", linestyle='--',
                     linewidth=1.0, color='lightblue')
        ax = plot_1D(time, clock_series-sigma_series, ax=ax, x_label="Time", linestyle='--',
                     linewidth=1.0, title=clock_bias.title, set_legend=True, y_label="clock bias [s]",
                     color='lightblue')
        ax.fill_between(time, clock_series + sigma_series, clock_series - sigma_series, color='lightblue', alpha=0.3)
        # Plot +/- sigma now...
        # ax = plot_1D(time, clock_bias., x_label="Time", y_label="clock [s]", title="Receiver Clock Bias")
        # ax = plot_1D(time, clock_bias. x_label="Time", y_label="clock [s]", title="Receiver Clock Bias")

    def _plot_errors(self, time, error, title, x_label, y_label, z_label):
        ax = None
        x_error = error[:, 0]
        y_error = error[:, 1]
        z_error = error[:, 2]
        norm = np.linalg.norm(error, axis=1)

        ax = plot_1D(time, norm, ax=ax, label="Norm")
        plot_1D(time, x_error, ax=ax, label=x_label)
        plot_1D(time, y_error, ax=ax, label=y_label)
        plot_1D(time, z_error, ax=ax, x_label="Time", y_label="Position error [m]", label=z_label,
                title=f"Estimation Error in {title}", set_legend=True)

    def _plot_dops(self, time, dop_ecef, dop_enu):

        x_dop = [dop[0] for dop in dop_ecef.data]
        y_dop = [dop[1] for dop in dop_ecef.data]
        z_dop = [dop[2] for dop in dop_ecef.data]
        t_dop = [dop[3] for dop in dop_ecef.data]
        geometry_dop = [dop[4] for dop in dop_ecef.data]
        position_dop = [dop[5] for dop in dop_ecef.data]

        east_dop = [dop[0] for dop in dop_enu.data]
        north_dop = [dop[1] for dop in dop_enu.data]
        up_dop = [dop[2] for dop in dop_enu.data]
        horizontal_dop = [dop[3] for dop in dop_enu.data]

        ax = plot_1D(time, geometry_dop, label="geometry")
        plot_1D(time, position_dop, ax=ax, label="position")
        plot_1D(time, horizontal_dop, ax=ax, x_label="Time", y_label="DOPs [m]", label="horizontal",
                title="Dilution of Precision", set_legend=True)

        ax = plot_1D(time, x_dop, label="x")
        plot_1D(time, y_dop, ax=ax, label="y")
        plot_1D(time, t_dop, ax=ax, label="t")
        plot_1D(time, z_dop, ax=ax, x_label="Time", y_label="DOPs [m]", label="z",
                title="Dilution of Precision in ECEF", set_legend=True)

        ax = plot_1D(time, east_dop, label="east")
        plot_1D(time, north_dop, ax=ax, label="north")
        plot_1D(time, up_dop, ax=ax, x_label="Time", y_label="DOPs [m]", label="up",
                title="Dilution of Precision in ENU", set_legend=True)


    def plot_iono(self, iono_tm):
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
