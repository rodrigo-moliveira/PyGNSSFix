import os

import numpy as np
import pandas as pd

from . import error
from .plots import plot_gnss

# TODO List
#   * todos os datamanager.get_data() têm que estar protegidos por um try catch
#   * nos plots quando faço datamanager.get_data() verificar se is_empty()


from .plots.utils import show_all
from ... import constants
from ...data_types.date import Epoch
from ...models.frames import cartesian2geodetic
from ...utils.str_utils import replace_whitespace_with_underscore


# TODO: tirar isto daqui


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
                self._plot(plot_dir)
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

    def _plot(self, plot_dir):
        # TODO: self.log.info("plotting...")
        # self._plot_obs(plot_dir)
        # self._plot_errors(plot_dir)
        # self._plot_residuals(plot_dir)
        # self._plot_clock_bias(plot_dir)
        # self._plot_clock_rate(plot_dir)
        # self._plot_iono(plot_dir)
        # self._plot_tropo(plot_dir)
        # self._plot_sat_availability(plot_dir)
        # self._plot_3D_traj(plot_dir)
        # self._plot_3D_errors(plot_dir)
        # self._plot_2D_errors(plot_dir)
        # self._plot_latlon(plot_dir)
        # self._plot_dops(plot_dir)
        self._skyplot(plot_dir)

        if self.config.get("performance_evaluation", "plot_configs", "show_plots"):
            show_all()
        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            pass

    def _plot_obs(self, plot_dir):
        """ Plot the GNSS observables """
        observations = self.data_manager.get_data("obs")
        ax_dict = plot_gnss.plot_observations(observations)

        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            for ax in ax_dict.values():
                plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
                plot_path = plot_dir / plot_name
                self.log.info(f"Saving figure {plot_path}")
                ax.figure.savefig(plot_path, format='png')

    def _plot_errors(self, plot_dir):
        """ Plot the estimation errors """
        time = self.data_manager.get_data("time")
        sow_array = time.data.iloc[:, 1].values
        week_array = time.data.iloc[:, 0].values
        time_array = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]

        ax1 = plot_gnss.plot_estimation_errors(time_array, self.vel_error["error_ecef"], self.vel_error["cov_ecef"],
                                               "Velocity Estimation Error in ECEF", "Velocity Error [m/s]", "X", "Y",
                                               "Z")
        ax2 = plot_gnss.plot_estimation_errors(time_array, self.vel_error["error_enu"], self.vel_error["cov_enu"],
                                               "Velocity Estimation Error in ENU", "Velocity Error [m/s]", "East",
                                               "North", "Up")
        ax3 = plot_gnss.plot_estimation_errors(time_array, self.pos_error["error_ecef"], self.pos_error["cov_ecef"],
                                               "Position Estimation Error in ECEF", "Position Error [m]", "X", "Y", "Z")
        ax4 = plot_gnss.plot_estimation_errors(time_array, self.pos_error["error_enu"], self.pos_error["cov_enu"],
                                               "Position Estimation Error in ENU", "Position Error [m]", "East",
                                               "North", "Up")
        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            for ax in [ax1, ax2, ax3, ax4]:
                plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
                plot_path = plot_dir / plot_name
                self.log.info(f"Saving figure {plot_path}")
                ax.figure.savefig(plot_path, format='png')

    def _plot_residuals(self, plot_dir):
        """ Plot the GNSS prefit and postfit residuals of the estimation process """
        ax1 = ax2 = ax3 = ax4 = None
        try:
            ax1 = plot_gnss.plot_estimation_residuals(self.data_manager.get_data("pr_rate_prefit_residuals"), "m/s")
        except ValueError:
            self.log.warn("Pseudorange Rate Prefit Residuals not found. Skipping this plot")
        try:
            ax2 = plot_gnss.plot_estimation_residuals(self.data_manager.get_data("pr_rate_postfit_residuals"), "m/s")
        except ValueError:
            self.log.warn("Pseudorange Rate Postfit Residuals not found. Skipping this plot")
        try:
            ax3 = plot_gnss.plot_estimation_residuals(self.data_manager.get_data("pr_prefit_residuals"), "m")
        except ValueError:
            self.log.warn("Pseudorange Prefit Residuals not found. Skipping this plot")
        try:
            ax4 = plot_gnss.plot_estimation_residuals(self.data_manager.get_data("pr_postfit_residuals"), "m")
        except ValueError:
            self.log.warn("Pseudorange Postfit Residuals not found. Skipping this plot")

        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            for ax in [ax1, ax2, ax3, ax4]:
                if ax is not None:
                    plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
                    plot_path = plot_dir / plot_name
                    self.log.info(f"Saving figure {plot_path}")
                    ax.figure.savefig(plot_path, format='png')

    def _plot_clock_bias(self, plot_dir):
        """
        Performs the following plots, if available:
            * Plot clock bias (master and slave clocks)
            * Plot the estimated ISB
            * Plot the Allan Deviation of the clock bias
        """
        clock_bias = isb = None
        try:
            clock_bias = self.data_manager.get_data("clock_bias")
            isb = self.data_manager.get_data("isb")
        except ValueError:
            pass

        if clock_bias is None or clock_bias.is_empty():
            self.log.warn("Clock Bias state not found. Skipping this plot")
            return
        if isb is None or isb.is_empty():
            self.log.warn("ISB state not found. This is expected if the scenario is single constellation. "
                          "Skipping this plot")

        ax1, ax2 = plot_gnss.plot_clock_bias(clock_bias, isb)
        ax3 = plot_gnss.plot_allan_deviation(clock_bias)
        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            for ax in [ax1, ax2, ax3]:
                if ax is not None:
                    plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
                    plot_path = plot_dir / plot_name
                    self.log.info(f"Saving figure {plot_path}")
                    ax.figure.savefig(plot_path, format='png')

    def _plot_clock_rate(self, plot_dir):
        """ Plot estimated clock bias rate """
        try:
            clock_rate = self.data_manager.get_data("clock_bias_rate")
            if clock_rate.is_empty():
                raise ValueError
        except ValueError:
            self.log.warn("Clock Bias Rate state not found. Skipping this plot")
            return

        ax = plot_gnss.plot_clock_bias_rate(clock_rate)

        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            if ax is not None:
                plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
                plot_path = plot_dir / plot_name
                self.log.info(f"Saving figure {plot_path}")
                ax.figure.savefig(plot_path, format='png')

    def _plot_iono(self, plot_dir):
        """ Plot Iono estimated states for all available satellites """
        try:
            iono = self.data_manager.get_data("iono")
            if iono.is_empty():
                raise ValueError
        except ValueError:
            self.log.warn("Satellite Iono state not found. Skipping this plot")
            return

        ax_list = plot_gnss.plot_iono_states(iono)
        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            for ax in ax_list:
                plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
                plot_path = plot_dir / plot_name
                self.log.info(f"Saving figure {plot_path}")
                ax.figure.savefig(plot_path, format='png')

    def _plot_tropo(self, plot_dir):
        """ Plot the estimated troposphere wet delay """
        try:
            tropo = self.data_manager.get_data("tropo_wet")
            if tropo.is_empty():
                raise ValueError
        except ValueError:
            self.log.warn("Satellite Tropo state not found. Skipping this plot")
            return

        ax = plot_gnss.plot_tropo_wet_delay(tropo)
        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
            plot_path = plot_dir / plot_name
            self.log.info(f"Saving figure {plot_path}")
            ax.figure.savefig(plot_path, format='png')

    def _plot_sat_availability(self, plot_dir):
        """ Plot satellite availability """
        try:
            prefit_residuals = self.data_manager.get_data("pr_prefit_residuals")
            position = self.data_manager.get_data("position")
            if prefit_residuals.is_empty() or position.is_empty():
                raise ValueError
        except ValueError:
            self.log.warn("Satellite information (prefit-residuals) not found. Skipping this plot")
            return

        obs_time = float(position.data.iloc[1, 1]) - float(position.data.iloc[0, 1])
        ax = plot_gnss.plot_satellite_availability(prefit_residuals, obs_time)
        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
            plot_path = plot_dir / plot_name
            self.log.info(f"Saving figure {plot_path}")
            ax.figure.savefig(plot_path, format='png')

    def _plot_3D_traj(self, plot_dir):
        """ plot the 3D estimated trajectory and covariance """
        position_df = self.data_manager.get_data("position")
        positions = position_df.data.loc[:, ['x', 'y', 'z']].values
        ax = plot_gnss.plot_3D_trajectory_with_avg_covariance(positions, self.pos_error["cov_ecef"],
                                                              true_position=self.true_pos,
                                                              x_label="X ECEF [m]", y_label="Y ECEF [m]",
                                                              z_label="Z ECEF [m]", title=position_df.title)

        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
            plot_path = plot_dir / plot_name
            self.log.info(f"Saving figure {plot_path}")
            ax.figure.savefig(plot_path, format='png')

    def _plot_3D_errors(self, plot_dir):
        """ Plot the absolute 3D trajectory and covariance, as well as the 3D errors in ENU and ECEF frames """
        ax1 = plot_gnss.plot_3D_trajectory_with_avg_covariance(self.pos_error["error_enu"], self.pos_error["cov_enu"],
                                                               true_position=[0, 0, 0],
                                                               x_label="East [m]", y_label="North [m]",
                                                               z_label="Up [m]",
                                                               title="3D ENU Error")
        ax2 = plot_gnss.plot_3D_trajectory_with_avg_covariance(self.pos_error["error_ecef"], self.pos_error["cov_ecef"],
                                                               true_position=[0, 0, 0],
                                                               x_label="X [m]", y_label="Y [m]",
                                                               z_label="Z [m]",
                                                               title="3D ECEF Error")
        for ax in [ax1, ax2]:
            if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
                plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
                plot_path = plot_dir / plot_name
                self.log.info(f"Saving figure {plot_path}")
                ax.figure.savefig(plot_path, format='png')

    def _plot_2D_errors(self, plot_dir):
        """ Plot the 2D errors in the ENU frame """
        ax = plot_gnss.plot_2D_trajectory(self.pos_error["error_enu"], self.pos_error["cov_enu"], true_pos=[0, 0],
                                          x_label="East [m]", y_label="North [m]", title="Horizontal Position Error",
                                          label="Error Estimates")

        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
            plot_path = plot_dir / plot_name
            self.log.info(f"Saving figure {plot_path}")
            ax.figure.savefig(plot_path, format='png')

    def _plot_latlon(self, plot_dir):
        """ Plot the Latitude-Longitude in a 2D figure"""
        position_df = self.data_manager.get_data("position")
        latlon = []

        def to_latlon(row):
            _lla = cartesian2geodetic(row.x, row.y, row.z)
            latlon.append([_lla[0] * constants.RAD2DEG, _lla[1] * constants.RAD2DEG])

        position_df.data.apply(to_latlon, axis=1)
        true_latlon = cartesian2geodetic(*self.true_pos) * np.array([constants.RAD2DEG, constants.RAD2DEG, 1])
        ax = plot_gnss.plot_2D_trajectory(latlon, None, true_pos=true_latlon, label="Position Estimates",
                                          x_label="Latitude [deg]", y_label="Longitude [deg]",
                                          title="Latitude-Longitude estimation")

        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
            plot_path = plot_dir / plot_name
            self.log.info(f"Saving figure {plot_path}")
            ax.figure.savefig(plot_path, format='png')

    def _plot_dops(self, plot_dir):
        """ Plot the DOP figures """
        try:
            dop_ecef = self.data_manager.get_data("dop_ecef")
            dop_enu = self.data_manager.get_data("dop_local")
            if dop_ecef.is_empty() or dop_enu.is_empty():
                raise ValueError
        except ValueError:
            self.log.warn("DOPs not found. Skipping this plot")
            return

        ax1, ax2, ax3 = plot_gnss.plot_dops(dop_ecef, dop_enu)

        for ax in [ax1, ax2, ax3]:
            if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
                plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
                plot_path = plot_dir / plot_name
                self.log.info(f"Saving figure {plot_path}")
                ax.figure.savefig(plot_path, format='png')

    def _skyplot(self, plot_dir):
        """ Plot the skyplot of the satellites in view """
        try:
            satellite_azel = self.data_manager.get_data("satellite_azel")
            if satellite_azel.is_empty():
                raise ValueError
        except ValueError:
            self.log.warn("Satellite Azimuth-Elevation information not found. Skipping this plot")
            return

        ax = plot_gnss.plot_skyplot(satellite_azel)

        if self.config.get("performance_evaluation", "plot_configs", "save_plots"):
            plot_name = f"{replace_whitespace_with_underscore(ax.get_title())}.png"
            plot_path = plot_dir / plot_name
            self.log.info(f"Saving figure {plot_path}")
            ax.figure.savefig(plot_path, format='png')
