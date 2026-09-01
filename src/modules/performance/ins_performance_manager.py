import os
import numpy as np
from src.common_log import get_logger, POST_PROC_LOG
from src.io.config import config_dict
from src.models.frames import lld2ecef
from src.models.plots.plot_gnss import plot_3D_trajectory
from src.models.plots.utils import show_all, plot_1D
from src.utils.str_utils import replace_whitespace_with_underscore


class INSPerformanceManager:
    """
    TODO: update this docstrings
    TODO: this needs some work
    INS Performance manager class.

    This class runs the post-processing / performance manager stage of the previously executed GNSS run.

    Three main tasks are performed:
        * Processing estimation errors
        * Processing residual analysis
        * Plotting data

    See Also:
        `algs/post_processing_gnss.py`: This script executes the PerformanceManager algorithm.
    """
    def __init__(self, data_manager):
        """
        Constructor of the PerformanceManager class.

        Args:
            data_manager(src.data_mng.ins.ins_data_mng.InsDataManager): Data Manager (InsDataManager)
                with all the information about the INS algorithm executed

        """
        self.data_manager = data_manager
        self.log = get_logger(POST_PROC_LOG)


    def process(self, output_dir):
        """ Main function of the performance manager. """

        # create post-processing folder, if not yet created
        post_proc_dir = f"{output_dir}/post_proc"
        if not os.path.exists(post_proc_dir):
            self.log.info(f"Creating directory to store computed errors {post_proc_dir}...")
            os.makedirs(post_proc_dir)

        # create plot folder, if not yet created
        plot_dir = f"{output_dir}/plot"
        if not os.path.exists(plot_dir):
            self.log.info(f"Creating directory to store plots {plot_dir}...")
            os.makedirs(plot_dir)

        #if config_dict.get("performance_evaluation", "run_error_computation"):
        #    self.log.info("Running Task: Processing estimation errors...")
        #    try:
        #        self._process_errors(post_proc_dir)
        #    except Exception as e:
        #        self.log.error(f"Error when processing estimation errors due to: {e}")

        #if config_dict.get("performance_evaluation", "run_residual_analysis"):
        #    self.log.info("Running Task: Processing residual analysis...")
        #    try:
        #        self._residual_analysis(post_proc_dir)
        #    except Exception as e:
        #        self.log.error(f"Error when running residual analysis due to: {e}")

        # plots
        if config_dict.get("performance_evaluation", "run_plot_manager"):
            self.log.info("Running Task: Plotting data...")
            try:
                self._plot(plot_dir)
            except Exception as e:
                self.log.error(f"Error when plotting data due to: {e}")

        self.log.info("Performance Manager finished execution.")

    def show_plots(self):
        """ Method to show all the plots. """
        if config_dict.get("performance_evaluation", "plot_configs", "show_plots"):
            show_all()

    def _plot(self, plot_dir):
        self.log.info("Plotting position data...")
        self._plot_pos("ref_pos", plot_dir, "ref")
        #self._plot_pos("pos", plot_dir, "est")

        self.log.info("Plotting velocity data...")
        self._plot_vel("ref_vel", plot_dir, "ref")
        #self._plot_vel("vel", plot_dir, "est")

        self.log.info("Plotting attitude data...")
        self._plot_att("ref_att", plot_dir, "ref")

        self.log.info("Plotting IMU data...")
        self._plot_imu("ref_gyro", plot_dir, "true")
        self._plot_imu("ref_accel", plot_dir, "true")
        self._plot_imu("gyro", plot_dir, "meas")
        self._plot_imu("accel", plot_dir, "meas")
        self._plot_imu_error(plot_dir)

        # Add Allan Variance plot


    def _plot_pos(self, data_str, plot_dir, label_prefix):
        try:
            pos_data = self.data_manager.get_data(data_str)
            if pos_data.is_empty():
                raise ValueError
        except:
            self.log.warning(f"Position data {data_str} dataframe not found or is empty. Skipping plot_pos")
            return
        try:
            ax_latlon, ax_alt, ax_ecef, ax_3d = self._plot_pos_inner(pos_data, label_prefix)
        except Exception as e:
            self.log.error(f"Unexpected error when performing plot_pos function for {data_str}: {e}")
            return

        # Save Figure
        for ax in [ax_latlon, ax_alt, ax_ecef, ax_3d]:
            self._save_figure(plot_dir, ax)

    # TODO: the functions here can be refactored

    def _plot_pos_inner(self, pos_df, label_prefix, ax_latlon=None, ax_alt=None, ax_ecef=None, ax_3d=None):
        # LLD plot
        lld = np.array(pos_df.to_data_array())
        time = np.array(pos_df.to_time_array())
        legend = pos_df.legend if pos_df.legend is not None else ["", "", ""]
        units = pos_df.units if pos_df.units is not None else ["", "", ""]
        ax_latlon = plot_1D(time, lld[:, 0], ax=ax_latlon, label=f"{label_prefix}_{legend[0]}")
        ax_latlon = plot_1D(time, lld[:, 1], ax=ax_latlon, title="Latitude, Longitude", label=f"{label_prefix}_{legend[1]}",
                            set_legend=True, y_label=units[0], x_label="Simulation Time (sec)")
        ax_alt = plot_1D(time, lld[:, 2], title="Altitude (Down)", label=f"{label_prefix}_{legend[2]}", set_legend=True,
                         y_label=units[2], x_label="Simulation Time (sec)", ax=ax_alt)

        # ECEF plot
        pos_ecef = lld2ecef(lld)
        ax_ecef = plot_1D(time, pos_ecef[:, 0], ax=ax_ecef, label=f"{label_prefix}_ref_x")
        ax_ecef = plot_1D(time, pos_ecef[:, 1], ax=ax_ecef, label=f"{label_prefix}_ref_y")
        ax_ecef = plot_1D(time, pos_ecef[:, 2], ax=ax_ecef, title="ECEF Cartesian Coordinates", label=f"{label_prefix}_ref_z",
                          set_legend=True, y_label="m", x_label="Simulation Time (sec)")

        # 3D ECEF plot
        ax_3d = plot_3D_trajectory(pos_ecef, label=label_prefix, x_label="x [m]", y_label="y [m]", z_label="z [m]",
                                   title="3D trajectory", ax=ax_3d)
        return ax_latlon, ax_alt, ax_ecef, ax_3d

    def _plot_vel(self, data_str, plot_dir, label_prefix):
        try:
            vel_data = self.data_manager.get_data(data_str)
            if vel_data.is_empty():
                raise ValueError
        except:
            self.log.warning(f"Velocity data {data_str} dataframe not found or is empty. Skipping plot_vel")
            return
        try:
            ax = self.plot_generic(vel_data, label_prefix)
        except Exception as e:
            self.log.error(f"Unexpected error when performing plot_vel for {data_str} function: {e}")
            return

        # Save Figure
        self._save_figure(plot_dir, ax)

    def _plot_att(self, data_str, plot_dir, label_prefix):
        try:
            att_data = self.data_manager.get_data(data_str)
            if att_data.is_empty():
                raise ValueError
        except:
            self.log.warning(f"Attitude (Euler) data {data_str} dataframe not found or is empty. Skipping plot_att")
            return
        try:
            ax = self.plot_generic(att_data, label_prefix)
        except Exception as e:
            self.log.error(f"Unexpected error when performing plot_att for {data_str} function: {e}")
            return

        # Save Figure
        self._save_figure(plot_dir, ax)

    def plot_generic(self, est_data, label_prefix):
        time = np.array(est_data.to_time_array())
        data_matrix = np.array(est_data.to_data_array())
        ax = None

        if not est_data.is_empty():
            ax = plot_1D(time, data_matrix[:, 0], ax=ax, label=f"{label_prefix}_{est_data.legend[0]}",
                            title=est_data.title,
                            x_label="Simulation Time (sec)", y_label=est_data.units[0],
                            set_legend=True)
            ax = plot_1D(time, data_matrix[:, 1], ax=ax, label=f"{label_prefix}_{est_data.legend[1]}",
                            title=est_data.title,
                            x_label="Simulation Time (sec)", y_label=est_data.units[1],
                            set_legend=True)
            ax = plot_1D(time, data_matrix[:, 2], ax=ax, label=f"{label_prefix}_{est_data.legend[2]}",
                            title=est_data.title,
                            x_label="Simulation Time (sec)", y_label=est_data.units[2],
                            set_legend=True)

        return ax

    def _plot_imu(self, data_str, plot_dir, label_prefix):
        try:
            imu_data = self.data_manager.get_data(data_str)
            if imu_data.is_empty():
                raise ValueError
        except:
            self.log.warning(f"IMU data {data_str} dataframe not found or is empty. Skipping plot_imu")
            return
        try:
            ax = self.plot_generic(imu_data, label_prefix)
        except Exception as e:
            self.log.error(f"Unexpected error when performing plot_imu for {data_str} function: {e}")
            return

        # Save Figure
        self._save_figure(plot_dir, ax)


    def _plot_imu_error(self, plot_dir):
        try:
            true_accel = self.data_manager.get_data("accel")
            meas_accel = self.data_manager.get_data("ref_accel")
            true_gyro = self.data_manager.get_data("gyro")
            meas_gyro = self.data_manager.get_data("ref_gyro")
            if true_accel.is_empty() or meas_accel.is_empty() or true_gyro.is_empty() or meas_gyro.is_empty():
                raise ValueError
        except:
            self.log.warning(f"IMU data dataframe not found or is empty. Skipping plot_imu_error")
            return
        try:
            ax_accel = self._plot_error(true_accel, meas_accel)
            ax_gyro = self._plot_error(true_gyro, meas_gyro)
        except Exception as e:
            self.log.error(f"Unexpected error when performing plot_imu_error for function: {e}")
            return

        # Save Figure
        self._save_figure(plot_dir, ax_accel)
        self._save_figure(plot_dir, ax_gyro)

    def _plot_error(self, true_df, meas_df):
        time = np.array(true_df.to_time_array())
        true_data = np.array(true_df.to_data_array())
        meas_data = np.array(meas_df.to_data_array())
        error_data = true_data - meas_data
        ax = None

        ax = plot_1D(time, error_data[:, 0], ax=ax, label=f"{true_df.legend[0]}",
                     x_label="Simulation Time (sec)", y_label=true_df.units[0],
                     set_legend=True)
        ax = plot_1D(time, error_data[:, 1], ax=ax, label=f"{true_df.legend[1]}",
                     x_label="Simulation Time (sec)", y_label=true_df.units[1],
                     set_legend=True)
        ax = plot_1D(time, error_data[:, 2], ax=ax, label=f"{true_df.legend[2]}",
                     title=f"Measurement Error - {meas_df.title}",
                     x_label="Simulation Time (sec)", y_label=true_df.units[2],
                     set_legend=True)
        return ax


    def _save_figure(self, plot_dir, ax):
        """ Saves the figure to a file. """
        if config_dict.get("performance_evaluation", "plot_configs", "save_plots") and ax is not None:
            plot_path = f"{plot_dir}/{replace_whitespace_with_underscore(ax.get_title())}.png"
            if not os.path.isfile(plot_path):
                self.log.info(f"Saving figure {plot_path}")
                ax.figure.savefig(plot_path, format='png')
            else:
                self.log.warning(f"File {plot_path} already exists.")