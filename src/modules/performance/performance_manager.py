import os

from . import error
from .plots import *
import pandas as pd
from src.data_types.date import Epoch
import allantools as at

# TODO List
#   * todos os datamanager.get_data() têm que estar protegidos por um try catch

import re

from ...io.rinex_parser.utils import RINEX_SATELLITE_SYSTEM, RINEX_OBS_TYPES_UNITS


# TODO: tirar isto daqui
def extract_constellations(input_string):
    # Define the regex pattern to match the master and slave constellations
    pattern = r'ISB\(master=(\w+) slave=(\w+)\)'

    # Search for the pattern in the input string
    match = re.search(pattern, input_string)

    # If a match is found, extract the master and slave values
    if match:
        master = match.group(1)
        slave = match.group(2)
        return master, slave
    else:
        return None, None


class PerformanceManager:
    """
    TODO: add docs
    """

    def __init__(self, data_manager, config, log):
        """TODO: add docstring"""
        self.data_manager = data_manager
        self.config = config
        self.log = log

        # Declaring variables to store error computations
        # True static position and velocity
        self.true_pos = None
        self.true_vel = None

        # Position and Velocity error dicts (estimation error + covariance)
        self.pos_error = {
            "error_ecef": None,
            "error_enu": None,
            "cov_ecef": None,
            "cov_enu": None
        }
        self.vel_error = {
            "error_ecef": None,
            "error_enu": None,
            "cov_ecef": None,
            "cov_enu": None
        }

    def process(self, output_dir):
        """ Main function of the performance manager. 3 tasks are performed:
            * Processing estimation errors
            * Processing residual analysis
            * Plotting data
         """

        # create post-processing folder, if not yet created
        post_proc_dir = output_dir / 'post_proc'
        if not os.path.exists(post_proc_dir):
            self.log.info(f"Creating directory to store computed errors {post_proc_dir}...")
            os.makedirs(post_proc_dir)

        # create plot folder, if not yet created
        plot_dir = output_dir / 'plot'
        if not os.path.exists(plot_dir):
            self.log.info(f"Creating directory to store plots {plot_dir}...")
            os.makedirs(plot_dir)

        if self.config.get("performance_evaluation", "run_error_computation"):
            self.log.info("Running Task: Processing estimation errors...")
            try:
                self._process_errors(post_proc_dir)
            except Exception as e:
                self.log.error(f"Error when processing estimation errors due to: {e}")

        if self.config.get("performance_evaluation", "run_residual_analysis"):
            self.log.info("Running Task: Processing residual analysis...")
            try:
                self._residual_analysis(post_proc_dir)
            except Exception as e:
                self.log.error(f"Error when running residual analysis due to: {e}")

        # plots
        if self.config.get("performance_evaluation", "run_plot_manager"):
            self.log.info("Running Task: Plotting data...")
            try:
                self._plot()
            except Exception as e:
                self.log.error(f"Error when plotting data due to: {e}")

    def _process_errors(self, post_proc_dir):
        """ Computes the estimation error and RMS statistics for position and velocity in ECEF and local ENU frames """
        # fetch estimated data
        position = self.data_manager.get_data("position")
        velocity = self.data_manager.get_data("velocity")

        # compute Root Mean Square Error (RMSE) in ECEF and ENU frames
        self.log.info("Computing estimation error and root mean square error...")
        if self.config.get("performance_evaluation", "error_configs", "static"):
            self.true_pos = self.config.get("performance_evaluation", "error_configs", "true_static_position")
            self.true_vel = [0, 0, 0]  # in static conditions, ECEF velocity is 0

            # compute position errors
            self.pos_error["error_ecef"], self.pos_error["error_enu"], self.pos_error["cov_ecef"], \
                self.pos_error["cov_enu"] = error.compute_error_static(position, self.true_pos, self.true_pos,
                                                                       local="ENU")

            # compute velocity errors
            self.vel_error["error_ecef"], self.vel_error["error_enu"], self.vel_error["cov_ecef"], \
                self.vel_error["cov_enu"] = error.compute_error_static(velocity, self.true_vel, self.true_pos,
                                                                       local="ENU")

        else:
            raise NotImplementedError(f"Dynamic error computation not yet implemented...")

        rms_pos_ecef = error.compute_rms_error(self.pos_error["error_ecef"])
        rms_pos_enu = error.compute_rms_error(self.pos_error["error_enu"])
        rms_vel_ecef = error.compute_rms_error(self.vel_error["error_ecef"])
        rms_vel_enu = error.compute_rms_error(self.vel_error["error_enu"])

        # save computed errors to files
        self._write_outputs(post_proc_dir, rms_pos_ecef, rms_pos_enu, rms_vel_ecef, rms_vel_enu)

    def _write_outputs(self, output_dir, rms_pos_ecef, rms_pos_enu, rms_vel_ecef, rms_vel_enu):
        """ Saves the error statistics to output files  """
        self.log.info("Saving computed errors to files...")

        time = self.data_manager.get_data("position").to_time_array()

        # Create DataFrame for error matrix
        pos_ecef_df = pd.DataFrame(self.pos_error["error_ecef"], columns=['Error x [m]', 'Error y [m]', 'Error z [m]'])
        vel_ecef_df = pd.DataFrame(self.vel_error["error_ecef"],
                                   columns=['Error vx [m/s]', 'Error vy [m/s]', 'Error vz [m/s]'])
        error_ecef_df = pd.concat([time, pos_ecef_df, vel_ecef_df], axis=1)
        pos_enu_df = pd.DataFrame(self.pos_error["error_enu"],
                                  columns=['Error East [m]', 'Error North [m]', 'Error Up [m]'])
        vel_enu_df = pd.DataFrame(self.vel_error["error_enu"],
                                  columns=['Error Vel East [m]', 'Error Vel North [m]', 'Error Vel Up [m]'])
        error_enu_df = pd.concat([time, pos_enu_df, vel_enu_df], axis=1)

        # Save the error dataframes
        self.log.info(f"Writing the ECEF estimation error to {output_dir / 'error_ecef.txt'}...")
        self.log.info(f"Writing the ENU estimation error to {output_dir / 'error_enu.txt'}...")
        error_ecef_df.to_csv(output_dir / "error_ecef.txt", index=False)
        error_enu_df.to_csv(output_dir / "error_enu.txt", index=False)

        # Overall Estimation Stats
        stats = "Root Mean Square Error:\n" \
                "\tECEF frame\n" \
                f"\t\tx = {rms_pos_ecef['x']} [m]\n" \
                f"\t\ty = {rms_pos_ecef['y']} [m]\n" \
                f"\t\tz = {rms_pos_ecef['z']} [m]\n" \
                f"\t\tvx = {rms_vel_ecef['x']} [m/s]\n" \
                f"\t\tvy = {rms_vel_ecef['y']} [m/s]\n" \
                f"\t\tvz = {rms_vel_ecef['z']} [m/s]\n" \
                "\n\n\tENU frame\n" \
                f"\t\teast = {rms_pos_enu['x']} [m]\n" \
                f"\t\tnorth = {rms_pos_enu['y']} [m]\n" \
                f"\t\tup = {rms_pos_enu['z']} [m]\n" \
                f"\t\tvel east = {rms_vel_enu['x']} [m/s]\n" \
                f"\t\tvel north = {rms_vel_enu['y']} [m/s]\n" \
                f"\t\tvel up = {rms_vel_enu['z']} [m/s]\n" \
                "\n\n\tOverall\n" \
                f"\t\tHorizontal Position Error = {rms_pos_enu['2D']} [m]\n" \
                f"\t\tHorizontal Velocity Error = {rms_vel_enu['2D']} [m/s]\n" \
                f"\t\t3D Position Error = {rms_pos_enu['3D']} [m]\n" \
                f"\t\t3D Velocity Error = {rms_vel_enu['3D']} [m/s]" \
                f"\n"

        self.log.info(f"Writing Root Mean Square Error statistics to {output_dir / 'EstimationStats.txt'}")
        f_rms_stats = open(output_dir / "EstimationStats.txt", "w")
        f_rms_stats.write(stats)
        self.log.info(stats)
        f_rms_stats.close()

    def _residual_analysis(self, output_dir):
        """ Performs residual analysis. Two statistical tests are performed to the GNSS post-fit residuals:
            * Chi Squared Test: Check if the residuals are consistent with white noise (mean of zero and uncorrelated).
            * Shapiro Test: Check if the residuals follow a normal distribution.

        The post-fit residuals of both pseudorange and pseudorange rate observations are evaluated, if available.
        The tests are performed for each pair of available satellites and datatypes
        """
        residuals_file = open(output_dir / "ResidualsStatistics.txt", "w")

        for res_name, data_name, df_name in zip(["pr_postfit_residuals", "pr_rate_postfit_residuals"],
                                                ["Pseudorange", "PseudorangeRate"],
                                                ["pr_estimated_params", "pr_rate_estimated_params"]):
            residuals = None
            empty = False
            try:
                residuals = self.data_manager.get_data(res_name)
            except ValueError:
                empty = True
            else:
                if residuals.is_empty():
                    empty = True

            if empty:
                self.log.warn(f"Skipping residual analysis for {data_name} observables")
                continue
            self.log.info(f"Performing residual analysis for {data_name} observables")

            # Group the data by satellite and data type
            grouped = residuals.data.groupby(['sat', 'data_type'])

            chi_sq_file = open(output_dir / f"ChiSquaredResidualAnalysis_{data_name}.txt", "w")
            shapiro_file = open(output_dir / f"ShapiroTestAnalysis.txt_{data_name}.txt", "w")

            df = self.config.get("performance_evaluation", "residual_configs", "chi_squared_test", df_name)
            chi_sq_alpha = self.config.get("performance_evaluation", "residual_configs", "chi_squared_test",
                                           "significance_level")
            shapiro_alpha = self.config.get("performance_evaluation", "residual_configs", "shapiro_test",
                                            "significance_level")

            self.log.info(f"Chi-Squared test for {data_name} observables: significance_level={chi_sq_alpha}, "
                          f"number of states={df}. Results are stored in {os.path.basename(chi_sq_file.name)}")
            self.log.info(f"Shapiro test for {data_name} observables: significance_level={shapiro_alpha}. "
                          f"Results are stored in {os.path.basename(shapiro_file.name)}")

            # Perform the tests for each satellite and data type combination
            for (sat, data_type), group in grouped:
                residual_array = group.iloc[:, 5].values
                mean = np.mean(residual_array)
                std = np.std(residual_array)

                residuals_file.write(f"{sat}, {data_type}, mean residual = {mean}, std = {std}\n")

                report_chi_sq = error.chi_squared_test(residual_array, df, chi_sq_alpha)
                report_shapiro = error.shapiro_test(residual_array, shapiro_alpha)

                chi_sq_file.write(f"Report for {sat} and {data_type}:\n{report_chi_sq}\n")
                shapiro_file.write(f"Report for {sat} and {data_type}:\n{report_shapiro}\n")

            chi_sq_file.close()
            shapiro_file.close()
        residuals_file.close()

    def _plot(self):
        # fetch estimated data
        # limpar aqui...
        self._plot_obs()

        self._plot_errors(time, self.vel_error["error_ecef"], self.vel_error["cov_ecef"],
                          "Velocity Estimation Error in ECEF", "Velocity Error [m/s]", "X", "Y", "Z")
        self._plot_errors(time, self.vel_error["error_enu"], self.vel_error["cov_enu"],
                          "Velocity Estimation Error in ENU", "Velocity Error [m/s]", "East", "North", "Up")

        self._plot_residuals(time, self.data_manager.get_data("vel_prefit_residuals"), "m/s")
        self._plot_residuals(time, self.data_manager.get_data("vel_postfit_residuals"), "m/s")

        position = self.data_manager.get_data("position")

        plot_satellite_availability(self.data_manager.get_data("prefit_residuals").data, x_label="Time",
                                    y_label="Number of Sats",
                                    title="Satellite Availability")

        self._plot_clock_bias(time)
        self._plot_clock_rate(time)

        plot_3D_trajectory_with_avg_covariance(position.to_data_array(), self.pos_error["cov_ecef"],
                                               true_position=self.true_pos,
                                               x_label="X ECEF [m]", y_label="Y ECEF [m]", z_label="Z ECEF [m]",
                                               title=position.title)

        plot_3D_trajectory_with_avg_covariance(self.pos_error["error_enu"], self.pos_error["cov_enu"],
                                               true_position=[0, 0, 0],
                                               x_label="East [m]", y_label="North [m]", z_label="Up [m]",
                                               title="3D ENU Error")

        self._plot_errors(time, self.pos_error["error_ecef"], self.pos_error["cov_ecef"],
                          "Position Estimation Error in ECEF", "Position Error [m]", "X", "Y", "Z")
        self._plot_errors(time, self.pos_error["error_enu"], self.pos_error["cov_enu"],
                          "Position Estimation Error in ENU", "Position Error [m]", "East", "North", "Up")
        self._plot_dops(time)

        plot_2D_trajectory(self.pos_error["error_enu"], self.pos_error["cov_enu"],
                           x_label="East [m]", y_label="North [m]", title="Horizontal Position Error")

        self._plot_residuals(time, self.data_manager.get_data("prefit_residuals"), "m")
        self._plot_residuals(time, self.data_manager.get_data("postfit_residuals"), "m")

        self._plot_iono(time)

        self._plot_tropo(time)
        self._plot_latlon()

        plot_skyplot(self.data_manager.get_data("satellite_azel").data)

        show_all()

    def _plot_allandev(self, clock_series, tau0):
        # Convert clock estimates to frequency deviations

        f = np.diff(clock_series) / tau0

        # Compute Allan deviation using allantools
        taus, adevs, _, _ = at.oadev(f, rate=1 / tau0, data_type='freq', taus='decade')

        # Plotting
        plt.loglog(taus, adevs)
        plt.xlabel('Tau (s)')
        plt.ylabel('Allan Deviation')
        plt.title('Allan Deviation of GNSS Clock Estimates')
        plt.grid(True)
        # plt.show()

    def _plot_latlon(self):
        latlon = compute_latlon(self.data_manager.get_data("position"))
        true_latlon = cartesian2geodetic(*self.true_pos) * np.array([constants.RAD2DEG, constants.RAD2DEG, 1])
        plot_2D_trajectory(latlon, None, true_pos=true_latlon,
                           x_label="Latitude [deg]", y_label="Longitude [deg]", title="Latitude-Longitude estimation")

    def _plot_residuals(self, time, residuals, units):

        ax = None
        # Group the data by satellite and data type
        grouped = residuals.data.groupby(['sat', 'data_type'])

        # Add a line for each satellite and data type combination
        for (sat, data_type), group in grouped:
            sow_array = group.iloc[:, 1].values
            week_array = group.iloc[:, 0].values
            residual_array = group.iloc[:, 5].values
            this_time = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]
            # TODO: why merge? simply take the `this_time`...
            merged_time = Epoch.merge_time_arrays(time, this_time)
            ax = plot_1D(merged_time, residual_array, ax=ax, x_label="Time", label=f'Residual for {sat} {data_type}',
                         set_legend=True, scatter=True, markersize=5)
        ax.set_ylabel(units)
        ax.set_title(residuals.title)

    def _plot_clock_bias(self, time):
        # TODO: find name of master constellation and add it to label
        # TODO: add clock bias + ISB and plot slave clock bias
        clock_bias = self.data_manager.get_data("clock_bias")
        master_clock = clock_bias.to_data_array()

        isb = self.data_manager.get_data("isb")
        isb_data = isb.to_data_array()
        # Extract master and slave constellations
        master, slave = extract_constellations(isb.data.columns[2])

        master_clock_series = master_clock[:, 0]
        sigma_series = np.sqrt(master_clock[:, 1])

        self._plot_allandev(master_clock_series, (time[1] - time[0]).total_seconds())

        isb_series = isb_data[:, 0]
        isb_sigma_series = np.sqrt(isb_data[:, 1])

        slave_clock_series = master_clock_series + isb_series

        ax = plot_1D(time, master_clock_series, x_label="Time", label=f'clock bias ({master})', color='blue')
        ax = plot_1D(time, slave_clock_series, x_label="Time", ax=ax, label=f'clock bias ({slave})', color='orange')
        ax = plot_1D(time, master_clock_series + sigma_series, ax=ax, x_label="Time", label="+/- sigma [s]",
                     linestyle='--',
                     linewidth=1.0, color='lightblue')
        ax = plot_1D(time, master_clock_series - sigma_series, ax=ax, x_label="Time", linestyle='--',
                     linewidth=1.0, title=f"{clock_bias.title} for {master}", set_legend=True, y_label="clock bias [s]",
                     color='lightblue')
        ax.fill_between(time, master_clock_series + sigma_series, master_clock_series - sigma_series, color='lightblue',
                        alpha=0.3)

        ax2 = plot_1D(time, isb_series, x_label="Time", title=f"{isb.title} (master={master}, slave={slave})",
                      label="ISB", set_legend=True, y_label="ISB [s]", color='blue')
        ax2.fill_between(time, isb_series + isb_sigma_series, isb_series - isb_sigma_series, label="sigma",
                         color='lightblue', alpha=0.3)

    def _plot_clock_rate(self, time):
        try:
            clock_rate = self.data_manager.get_data("clock_bias_rate")

            # Set colors for ±σ signals
            colors = {'GPS': 'blue', 'GAL': 'green'}

            # Get unique constellations
            unique_constellations = clock_rate.data['constellation'].unique()

            # Create sub-dataframes for each constellation
            sub_dataframes = {constellation: clock_rate.data[clock_rate.data['constellation'] == constellation]
                              for constellation in unique_constellations}

            # Plot clock rate and covariance for each constellation
            ax = None
            for constellation, sub_df in sub_dataframes.items():
                rate = sub_df['clock rate']
                sigma = np.sqrt(sub_df['cov'])

                ax = plot_1D(time, rate, ax=ax, x_label="Time", label=f'clock rate for {constellation}',
                             color=colors[constellation])
                ax = plot_1D(time, rate + sigma, ax=ax, x_label="Time", label=f"+/- sigma [s/s] for {constellation}",
                             linestyle='--',
                             linewidth=1.0, color=f'light{colors[constellation]}')
                ax = plot_1D(time, rate - sigma, ax=ax, x_label="Time", linestyle='--',
                             linewidth=1.0, title=clock_rate.title, set_legend=True, y_label="clock rate [s/s]",
                             color=f'light{colors[constellation]}')
                ax.fill_between(time, rate + sigma, rate - sigma, color=f'light{colors[constellation]}',
                                alpha=0.3)

        except Exception as e:
            print(f"Not plotting clock bias rate due to: {e}")

    def _plot_errors(self, time, error, cov, title, y_axis, x_label, y_label, z_label):
        ax = None
        x_error = error[:, 0]
        y_error = error[:, 1]
        z_error = error[:, 2]
        x_sigma = np.array([np.sqrt(x[0, 0]) for x in cov])
        y_sigma = np.array([np.sqrt(x[1, 1]) for x in cov])
        z_sigma = np.array([np.sqrt(x[2, 2]) for x in cov])

        norm = np.linalg.norm(error, axis=1)
        ax = plot_1D(time, norm, ax=ax, label="Norm", linewidth=2.0, linestyle="solid", color="k")
        plot_1D(time, x_error, ax=ax, label=x_label, linewidth=2.0, linestyle="solid", color="r")
        plot_1D(time, y_error, ax=ax, label=y_label, linewidth=2.0, linestyle="solid", color="b")
        plot_1D(time, z_error, ax=ax, label=z_label, linewidth=2.0, linestyle="solid", color="g")

        plot_1D(time, x_sigma, ax=ax, label=f"{x_label} sigma", linewidth=1.0, linestyle="dashed", color="r")
        plot_1D(time, y_sigma, ax=ax, label=f"{y_label} sigma", linewidth=1.0, linestyle="dashed", color="b")
        plot_1D(time, -x_sigma, ax=ax, linewidth=1.0, linestyle="dashed", color="r")
        plot_1D(time, -y_sigma, ax=ax, linewidth=1.0, linestyle="dashed", color="b")
        plot_1D(time, z_sigma, ax=ax, label=f"{z_label} sigma", linewidth=1.0, linestyle="dashed", color="g")
        plot_1D(time, -z_sigma, ax=ax, x_label="Time", linewidth=1.0, linestyle="dashed", color="g",
                y_label=y_axis, title=title, set_legend=True)

    def _plot_dops(self, time):
        dop_ecef = self.data_manager.get_data("dop_ecef").to_data_array()
        dop_enu = self.data_manager.get_data("dop_local").to_data_array()

        x_dop = [dop[0] for dop in dop_ecef]
        y_dop = [dop[1] for dop in dop_ecef]
        z_dop = [dop[2] for dop in dop_ecef]
        t_dop = [dop[3] for dop in dop_ecef]
        geometry_dop = [dop[4] for dop in dop_ecef]
        position_dop = [dop[5] for dop in dop_ecef]

        east_dop = [dop[0] for dop in dop_enu]
        north_dop = [dop[1] for dop in dop_enu]
        up_dop = [dop[2] for dop in dop_enu]
        horizontal_dop = [dop[3] for dop in dop_enu]

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

    def _plot_iono(self, time):
        iono = self.data_manager.get_data("iono")
        # Group the data by satellite and data type
        grouped = iono.data.groupby(['sat'])

        # Add a line for each satellite and data type combination
        for (sat,), group in grouped:
            sow_array = group.iloc[:, 1].values
            week_array = group.iloc[:, 0].values
            iono_array = group.iloc[:, 3].values

            sigma_array = np.sqrt(group.iloc[:, 4].values)
            this_time = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]

            merged_time = Epoch.merge_time_arrays(time, this_time)

            ax = plot_1D(merged_time, iono_array, label="iono", linewidth=2.0, linestyle="solid", color="b")
            plot_1D(merged_time, iono_array + sigma_array, ax=ax, linewidth=1.0, linestyle="dashed", color="lightblue")
            plot_1D(merged_time, iono_array - sigma_array, ax=ax, label=f"sigma", linewidth=1.0, linestyle="dashed",
                    color="lightblue",
                    x_label="Time", y_label="Ionosphere [m]", title=f"{iono.title} for sat {sat}",
                    set_legend=True)
            ax.fill_between(merged_time, iono_array - sigma_array, iono_array + sigma_array,
                            color='lightblue', alpha=0.3)

    def _plot_tropo(self, time):
        tropo = self.data_manager.get_data("tropo_wet")
        tropo_data = tropo.to_data_array()

        tropo_array = tropo_data[:, 0]
        sigma_series = np.sqrt(tropo_data[:, 1])

        ax = plot_1D(time, tropo_array, label="tropo", linewidth=2.0, linestyle="solid", color="b")
        plot_1D(time, tropo_array + sigma_series, ax=ax, linewidth=1.0, linestyle="dashed", color="lightblue")
        plot_1D(time, tropo_array - sigma_series, ax=ax, label=f"sigma", linewidth=1.0, linestyle="dashed",
                color="lightblue", x_label="Time", y_label="Wet Delay [m]", title=f"{tropo.title}",
                set_legend=True)
        ax.fill_between(time, tropo_array - sigma_series, tropo_array + sigma_series,
                        color='lightblue', alpha=0.3)

    def _plot_obs(self):
        observations = self.data_manager.get_data("obs")

        # Group the data by satellite and data type
        grouped = observations.data.groupby(['Satellite', 'Observation'])

        ax_dict = {}

        # Add a line for each satellite and data type combination
        for (sat, data_type), group in grouped:
            key = f"{sat[0]}_{data_type}"
            ax = ax_dict.get(key, None)

            sow_array = group.iloc[:, 1].values
            week_array = group.iloc[:, 0].values
            obs_array = group.iloc[:, 4].values
            this_time = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]
            # TODO: why merge? simply take the `this_time`...
            ax = plot_1D(this_time, obs_array, ax=ax, x_label="Time", label=f'{sat} {data_type}',
                         set_legend=True, scatter=True, markersize=5)
            ax_dict[key] = ax

        for key, ax in ax_dict.items():
            const = RINEX_SATELLITE_SYSTEM[key[0]]
            data_type = re.match(r'^[A-Za-z]+', key[2:]).group()
            units = RINEX_OBS_TYPES_UNITS.get(data_type, 'XX')

            ax.set_ylabel(f"{data_type} [{units}]")
            ax.set_title(f'Observation {key[2:]} for {const}')
