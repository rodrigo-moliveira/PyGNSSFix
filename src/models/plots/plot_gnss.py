""" Module with functions for GNSS plotting """
import re

import allantools as at
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import Axes3D

from src.data_mng.csv.csv_data import CSVData
from src.io.rinex_parser.utils import RINEX_SATELLITE_SYSTEM, RINEX_OBS_TYPES_UNITS
from src.utils.str_utils import extract_constellations
from src.models.plots import skyplot
from .utils import *
from ...data_types.date import Epoch


def plot_observations(observations: CSVData):
    """
    Plots the GNSS Observables.

    Args:
        observations(CSVData): input CSVData object with the GNSS Observables dataset
    Returns:
        dict: dict with the different matplotlib.Axes for each constellation and observable datatype
    """
    if not observations.is_empty():

        # Group the data by satellite and data type
        grouped = observations.data.groupby(['sat', 'data_type'])

        ax_dict = {}

        # Add a line for each satellite and data type combination
        for (sat, data_type), group in grouped:
            key = f"{sat[0]}_{data_type}"
            ax = ax_dict.get(key, None)

            sow_array = group.iloc[:, 1].values
            week_array = group.iloc[:, 0].values
            obs_array = group['value'].values
            time_array = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]

            ax = plot_1D(time_array, obs_array, ax=ax, x_label="Time", label=f'{sat} {data_type}',
                         set_legend=True, scatter=True, markersize=5)
            ax_dict[key] = ax

        for key, ax in ax_dict.items():
            const = RINEX_SATELLITE_SYSTEM[key[0]]
            data_type = re.match(r'^[A-Za-z]+', key[2:]).group()
            units = RINEX_OBS_TYPES_UNITS.get(data_type, 'XX')

            ax.set_ylabel(f"{data_type} [{units}]")
            ax.set_title(f'Observation {key[2:]} for {const}')
        return ax_dict


def plot_estimation_errors(time_array, error, cov, title, y_axis, dim1_label, dim2_label, dim3_label):
    """
    Plots the position/velocity estimation errors

    Args:
        time_array(list): list of Epoch objects of size N
        error(numpy.ndarray): matrix of shape Nx3 with the estimation error for each epoch
        cov( list[numpy.ndarray]): list with length N with the estimation covariance matrices (size 3x3)
        title(str): title for the figure
        y_axis(str): string with the label for the y-axis of the figure
        dim1_label(str): label with the description of the 1st dimension of the vector (x-axis)
        dim2_label(str): label with the description of the 2nd dimension of the vector (y-axis)
        dim3_label(str): label with the description of the 3rd dimension of the vector (z-axis)
    Returns:
        matplotlib.pyplot.Axes: returns the Axes object for this figure
    """
    ax = None
    x_error = error[:, 0]
    y_error = error[:, 1]
    z_error = error[:, 2]
    x_sigma = np.array([np.sqrt(x[0, 0]) for x in cov])
    y_sigma = np.array([np.sqrt(x[1, 1]) for x in cov])
    z_sigma = np.array([np.sqrt(x[2, 2]) for x in cov])

    norm = np.linalg.norm(error, axis=1)
    ax = plot_1D(time_array, norm, ax=ax, label="Norm", linewidth=2.0, linestyle="solid", color="k")
    plot_1D(time_array, x_error, ax=ax, label=dim1_label, linewidth=2.0, linestyle="solid", color="r")
    plot_1D(time_array, y_error, ax=ax, label=dim2_label, linewidth=2.0, linestyle="solid", color="b")
    plot_1D(time_array, z_error, ax=ax, label=dim3_label, linewidth=2.0, linestyle="solid", color="g")

    plot_1D(time_array, x_sigma, ax=ax, label=f"{dim1_label} sigma", linewidth=1.0, linestyle="dashed", color="r")
    plot_1D(time_array, y_sigma, ax=ax, label=f"{dim2_label} sigma", linewidth=1.0, linestyle="dashed", color="b")
    plot_1D(time_array, -x_sigma, ax=ax, linewidth=1.0, linestyle="dashed", color="r")
    plot_1D(time_array, -y_sigma, ax=ax, linewidth=1.0, linestyle="dashed", color="b")
    plot_1D(time_array, z_sigma, ax=ax, label=f"{dim3_label} sigma", linewidth=1.0, linestyle="dashed", color="g")
    plot_1D(time_array, -z_sigma, ax=ax, x_label="Time", linewidth=1.0, linestyle="dashed", color="g",
            y_label=y_axis, title=title, set_legend=True)
    return ax


def plot_estimation_residuals(residuals: CSVData, units):
    """
    Plot the GNSS prefit and postfit residuals of the estimation process.

    Args:
        residuals(CSVData): input CSVData object with the GNSS Observable Residuals dataset
        units(str): string with the units for the residuals (ex: m or m/s)
    Returns:
        matplotlib.pyplot.Axes: returns the Axes object for this figure
    """
    ax = None

    if not residuals.is_empty():

        # Group the data by satellite and data type
        grouped = residuals.data.groupby(['sat', 'data_type'])

        # Add a line for each satellite and data type combination
        for (sat, data_type), group in grouped:
            sow_array = group.iloc[:, 1].values
            week_array = group.iloc[:, 0].values
            residual_array = group.iloc[:, 5].values
            time_array = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]
            ax = plot_1D(time_array, residual_array, ax=ax, x_label="Time", label=f'Residual for {sat} {data_type}',
                         set_legend=True, scatter=True, markersize=5)
        ax.set_ylabel(units)
        ax.set_title(residuals.title)
    return ax


def plot_clock_bias(clock_bias: CSVData, isb: CSVData):
    """
    Plot the estimated receiver clock bias state. If the ISB is available:
        * the clock for the slave constellation is also plotted together with the master clock
        * an additional plot with the ISB estimated state is depicted as well.

    Args:
        clock_bias(CSVData): input CSVData object with the estimated receiver clock bias (master clock)
        isb(str): input CSVData object with the ISB of the master clock wrt the slave clock
    Returns:
        matplotlib.pyplot.Axes: returns the Axes objects for the two figures
    """
    ax1, ax2 = None, None
    has_slave = False
    if clock_bias is None or clock_bias.is_empty():
        return
    if isb is not None and not isb.is_empty():
        has_slave = True

    # Extract master and slave constellations
    if has_slave:
        master, slave = extract_constellations(clock_bias.data.columns[2], isb.data.columns[2])
    else:
        master, slave = extract_constellations(clock_bias.data.columns[2], "")

    sow_array = clock_bias.data.iloc[:, 1].values
    week_array = clock_bias.data.iloc[:, 0].values
    time_array = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]

    # fetch clock data
    master_clock_array = clock_bias.data.iloc[:, 2].values
    clock_cov_array = clock_bias.data.iloc[:, 3].values
    clock_sigma_array = np.sqrt(clock_cov_array)

    if has_slave:
        # fetch ISB data
        isb_array = isb.data.iloc[:, 2].values
        isb_cov_array = isb.data.iloc[:, 3].values
        isb_sigma_array = np.sqrt(isb_cov_array)

        # compute slave clock array
        slave_clock_array = master_clock_array + isb_array
    else:
        isb_array = None
        isb_sigma_array = None
        slave_clock_array = None

    # plot master + slave clocks (if available)
    ax1 = plot_1D(time_array, master_clock_array, x_label="Time", label=f'clock bias (master={master})', color='blue')
    if has_slave:
        ax1 = plot_1D(time_array, slave_clock_array, x_label="Time", ax=ax1, label=f'clock bias (slave={slave})',
                      color='orange')
    ax1 = plot_1D(time_array, master_clock_array + clock_sigma_array, ax=ax1, x_label="Time", label="+/- sigma [s]",
                  linestyle='--', linewidth=1.0, color='lightblue')
    ax1 = plot_1D(time_array, master_clock_array - clock_sigma_array, ax=ax1, x_label="Time", linestyle='--',
                  linewidth=1.0, title=f"{clock_bias.title} for master={master}", set_legend=True,
                  y_label="clock bias [s]", color='lightblue')
    ax1.fill_between(time_array, master_clock_array + clock_sigma_array, master_clock_array - clock_sigma_array,
                     color='lightblue', alpha=0.3)

    if has_slave:
        # plot isb in a separate plot
        ax2 = plot_1D(time_array, isb_array, x_label="Time", title=f"{isb.title} (master={master}, slave={slave})",
                      label="ISB", set_legend=True, y_label="ISB [s]", color='blue')
        ax2.fill_between(time_array, isb_array + isb_sigma_array, isb_array - isb_sigma_array, label="sigma",
                         color='lightblue', alpha=0.3)
    return ax1, ax2


def plot_allan_deviation(clock_bias):
    """
    Plot the Allan Deviation of the estimated receiver clock bias state.
    The Allan Deviation is computed using the function :py:func:`allantools.oadev`.
    
    Args:
        clock_bias(CSVData): input CSVData object with the estimated receiver clock bias (master clock)
    Returns:
        matplotlib.pyplot.Axes: returns the Axes object for the figure
    """
    if clock_bias is None or clock_bias.is_empty():
        return
    sow_array = clock_bias.data.iloc[:, 1].values
    week_array = clock_bias.data.iloc[:, 0].values
    time_array = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]
    master_clock_array = clock_bias.data.iloc[:, 2].values

    # Convert clock estimates to frequency deviations
    tau0 = (time_array[1] - time_array[0]).total_seconds()
    f = np.diff(master_clock_array) / tau0

    # Compute Allan deviation using allantools
    taus, adevs, _, _ = at.oadev(f, rate=1 / tau0, data_type='freq', taus='decade')

    # Plotting
    ax = loglog(taus, adevs, x_label='Tau (s)', y_label='Allan Deviation',
                title='Allan Deviation of GNSS Clock Estimates', grid="True")
    return ax


def plot_clock_bias_rate(clock_rate: CSVData):
    """
    Plot the estimated receiver clock bias rate state.

    Args:
        clock_rate(CSVData): input CSVData object with the estimated receiver clock bias rate (master clock)
    Returns:
        matplotlib.pyplot.Axes: returns the Axes object for the figure
    """
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
        sow_array = sub_df.iloc[:, 1].values
        week_array = sub_df.iloc[:, 0].values
        time_array = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]
        rate = sub_df['clock rate']
        sigma = np.sqrt(sub_df['cov'])

        ax = plot_1D(time_array, rate, ax=ax, x_label="Time", label=f'clock rate for {constellation}',
                     color=colors[constellation])
        ax = plot_1D(time_array, rate + sigma, ax=ax, x_label="Time", label=f"+/- sigma [s/s] for {constellation}",
                     linestyle='--',
                     linewidth=1.0, color=f'light{colors[constellation]}')
        ax = plot_1D(time_array, rate - sigma, ax=ax, x_label="Time", linestyle='--',
                     linewidth=1.0, title=clock_rate.title, set_legend=True, y_label="clock rate [s/s]",
                     color=f'light{colors[constellation]}')
        ax.fill_between(time_array, rate + sigma, rate - sigma, color=f'light{colors[constellation]}',
                        alpha=0.3)
    return ax


def plot_tropo_wet_delay(tropo: CSVData):
    """
    Plot the estimated receiver clock bias rate state.

    Args:
        tropo(CSVData): input CSVData object with the estimated tropo state
    Returns:
        matplotlib.pyplot.Axes: returns the Axes object for the figure
    """
    sow_array = tropo.data.iloc[:, 1].values
    week_array = tropo.data.iloc[:, 0].values
    time_array = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]
    tropo_array = tropo.data['tropo_wet']
    sigma = np.sqrt(tropo.data['cov'])

    ax = plot_1D(time_array, tropo_array, label="tropo", linewidth=2.0, linestyle="solid", color="b")
    plot_1D(time_array, tropo_array + sigma, ax=ax, linewidth=1.0, linestyle="dashed", color="lightblue")
    plot_1D(time_array, tropo_array - sigma, ax=ax, label=f"sigma", linewidth=1.0, linestyle="dashed",
            color="lightblue", x_label="Time", y_label="Wet Delay [m]", title=f"{tropo.title}",
            set_legend=True)
    ax.fill_between(time_array, tropo_array - sigma, tropo_array + sigma,
                    color='lightblue', alpha=0.3)
    return ax


def plot_iono_states(iono: CSVData):
    """
    Plot the estimated iono states for each satellite.

    Args:
        iono(CSVData): input CSVData object with the estimated states for all available satellites
    Returns:
        list[matplotlib.pyplot.Axes]: returns a list of Axes objects for each available satellite iono state
    """
    grouped = iono.data.groupby(['sat'])
    ax_list = []

    # Add a line for each satellite and data type combination
    for (sat,), group in grouped:
        sow_array = group.iloc[:, 1].values
        week_array = group.iloc[:, 0].values
        iono_array = group.iloc[:, 3].values

        sigma_array = np.sqrt(group.iloc[:, 4].values)
        time_array = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]

        ax = plot_1D(time_array, iono_array, label="iono", linewidth=2.0, linestyle="solid", color="b")
        plot_1D(time_array, iono_array + sigma_array, ax=ax, linewidth=1.0, linestyle="dashed", color="lightblue")
        plot_1D(time_array, iono_array - sigma_array, ax=ax, label=f"sigma", linewidth=1.0, linestyle="dashed",
                color="lightblue",
                x_label="Time", y_label="Ionosphere [m]", title=f"{iono.title} for sat {sat}",
                set_legend=True)
        ax.fill_between(time_array, iono_array - sigma_array, iono_array + sigma_array,
                        color='lightblue', alpha=0.3)
        ax_list.append(ax)
    return ax_list


def plot_satellite_availability(residuals: CSVData, obs_time):
    """
    Plot the satellite availability.

    Args:
        residuals(CSVData): input CSVData object with the pseudorange residuals for each datatype, satellite and epoch
        obs_time(float): default time between observations, in seconds
    Returns:
        matplotlib.pyplot.Axes: returns the Axes object
    """
    df = residuals.data.copy()

    # Combine Satellite and Signal for the y-axis labels
    df['Satellite_Signal'] = df['sat'] + '_' + df['data_type']

    def to_epoch_column(_row):
        _week = _row.iloc[0]
        _sow = _row.iloc[1]
        return Epoch.from_gnss_time(_week, _sow, scale="GPST")

    # Convert time into a format suitable for plotting
    df['Time'] = df.apply(to_epoch_column, axis=1)

    # Identify start and end points of intervals
    intervals = []

    for satellite_signal in df['Satellite_Signal'].unique():
        subset = df[df['Satellite_Signal'] == satellite_signal]
        subset = subset.sort_values(by='Time')

        start_time = subset.iloc[0]['Time']
        end_time = subset.iloc[0]['Time']

        for i in range(1, len(subset)):
            current_time = subset.iloc[i]['Time']
            if (current_time - end_time).total_seconds() > obs_time:
                intervals.append((start_time, end_time, satellite_signal))
                start_time = current_time
            end_time = current_time

        intervals.append((start_time, end_time, satellite_signal))

    # Create a reduced DataFrame with intervals
    interval_df = pd.DataFrame(intervals, columns=['Start_Time', 'End_Time', 'Satellite_Signal'])

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(15, 10))

    # Get unique Satellite_Signal combinations for the y-axis
    satellite_signal_combinations = interval_df['Satellite_Signal'].unique()
    y_labels = {label: i for i, label in enumerate(satellite_signal_combinations)}

    # Define the two alternating colors
    colors = ['blue', 'red']

    # Map each Satellite_Signal to a specific color
    color_map = {satellite_signal: colors[i % len(colors)] for i, satellite_signal in
                 enumerate(satellite_signal_combinations)}

    # Plot each interval as a line with the color based on Satellite_Signal
    for _, row in interval_df.iterrows():
        y = y_labels[row['Satellite_Signal']]
        color = color_map[row['Satellite_Signal']]  # Consistent color for each Satellite_Signal
        ax.plot([row['Start_Time'], row['End_Time']], [y, y], color=color, linewidth=2)

    # Customize the y-axis with colored labels
    ax.set_yticks(range(len(y_labels)))
    ax.set_yticklabels(y_labels.keys())

    # Set color for each y-axis label
    for label in ax.get_yticklabels():
        satellite_signal = label.get_text()
        label.set_color(color_map[satellite_signal])

    # Set axis labels and title
    ax.set_xlabel('Time')
    ax.set_ylabel('Number of Sats')
    ax.set_title('Satellite Signal Availability Over Time')
    ax.grid(True)

    return ax


def plot_3D_trajectory_with_avg_covariance(data_points, cov,
                                           true_position=None, x_label: str = "X-axis",
                                           y_label: str = "Y-axis",
                                           z_label: str = "Z-axis",
                                           title: str = "3D Trajectory with Covariance") -> Axes3D:
    """
    Plots the 3D trajectory of the estimated positions. If available, an average covariance ellipsoid and the true
    static position are also depicted.

    Args:
        data_points (list[tuple[float, float, float]]): List of estimated positions
        cov (list[np.ndarray]): List of estimated covariance matrices
        true_position (tuple[float, float, float]): True static position
        x_label(str): Label for the X-axis
        y_label(str): Label for the Y-axis
        z_label(str): Label for the Z-axis
        title(str): Title of the plot
    Returns:
        Axes3D: returns the Axes object with the 3D plot
    Raises:
        ValueError: an exception is raised if the input arguments do not have the required shape
    """
    # Validate data_points
    if not all(len(point) == 3 for point in data_points):
        raise ValueError("Each data point must have 3 elements (pos_x, pos_y, pos_z).")

    # Validate true_position
    if true_position and len(true_position) != 3:
        raise ValueError("True position must be a 3-dimensional point.")

    # If no axis is provided, create a new one
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # Plot estimated positions
    positions = np.array([[x[0], x[1], x[2]] for x in data_points])
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
               label="Estimated Position", marker='o', color='k', s=20)

    # Calculate the average covariance matrix
    mean_cov_matrix = np.mean(np.array([np.array(c) for c in cov]), axis=0)

    # Plot the average covariance ellipsoid at the mean position
    mean_position = np.mean(positions, axis=0)
    plot_covariance_ellipsoid(ax, mean_position, mean_cov_matrix)

    # Plot mean estimated position
    ax.scatter(mean_position[0], mean_position[1], mean_position[2],
               label="Mean Estimated Position", s=80, marker='D', color='b')

    # Plot true position if provided
    if true_position:
        ax.scatter(true_position[0], true_position[1], true_position[2],
                   label="True Static Position", marker='*', color='r', s=80)

    # Set labels and title
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.set_title(title)

    # Format the ticks
    ax.ticklabel_format(useOffset=False)

    # Show legend
    ax.legend()

    return ax


def plot_covariance_ellipsoid(ax, pos, cov, n_std=1.0):
    """
    Plot a 3D covariance ellipsoid in the provided Axes3D object.

    Args:
        ax (Axes3D): Axes on which to plot.
        pos (tuple[float, float, float]): Center position of the ellipsoid.
        cov (numpy.ndarray): Covariance matrix (3x3).
        n_std (float): Number of standard deviations to scale the ellipsoid. For example if n_std is 1, the 1-sigma
                        surface of the covariance is plotted
    """
    # Eigenvalues and eigenvectors of the covariance matrix
    eigvals, eigvecs = np.linalg.eigh(cov)
    # Sort the eigenvalues and eigenvectors
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # Compute radii of the ellipsoid
    radii = n_std * np.sqrt(eigvals)

    # Generate points on a unit sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 50)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    # Transform the unit sphere into the ellipsoid
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            [x[i, j], y[i, j], z[i, j]] = np.dot(eigvecs,
                                                 [x[i, j] * radii[0], y[i, j] * radii[1], z[i, j] * radii[2]]) + pos

    # Plot the surface
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='c', alpha=0.1, edgecolor='k', linewidth=0.5)


def plot_2D_trajectory(data, cov=None, true_pos=None, label="", x_label="", y_label="", title="", ax=None):
    """
    Plots the 2D trajectory of the estimated positions with an average covariance ellipse.

    Args:
        data (list[tuple[float, float]]): List of estimated positions in the ENU frame.
        cov (list[numpy.ndarray]): List of the estimated covariance matrices.
        true_pos (tuple[float, float]): True position
        label (str): Label for the estimated positions
        x_label (str): Label for the X-axis
        y_label (str): Label for the Y-axis
        title (str): Title of the plot
        ax (matplotlib.pyplot.Axes): pre-existing Axes object to plot write the plots

    Returns:
        matplotlib.pyplot.Axes: The axes with the plot
    """
    if ax is None:
        fig, ax = plt.subplots()

    x = [d[0] for d in data]
    y = [d[1] for d in data]

    mean_error_matrix = np.mean(np.array([np.array(c) for c in data]), axis=0)

    # Plot the estimated positions
    ax.scatter(x, y, linewidth=0.2, marker='o', label=label)
    ax.scatter(mean_error_matrix[0], mean_error_matrix[1], linewidth=2, marker='D', label="Mean Estimate")
    if true_pos is not None:
        ax.scatter(true_pos[0], true_pos[1], linewidth=2, marker='*', label="True Position")

    if cov is not None:
        # Calculate the mean covariance matrix
        mean_cov_matrix = np.mean(np.array([np.array(c) for c in cov]), axis=0)

        # Plot the covariance ellipse
        plot_covariance_ellipse(ax, mean_cov_matrix[0:2, 0:2], center=mean_error_matrix[0:2], color='r')

    # Set labels and title
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    # Show legend
    ax.legend(loc='upper right')

    return ax


def plot_covariance_ellipse(ax, cov_matrix, center=(0, 0), color='r'):
    """
    Plots a 2D covariance ellipse based on the covariance matrix.

    Args:
        ax (matplotlib.axes.Axes): The axes to plot on.
        cov_matrix (numpy.ndarray): The 2x2 covariance matrix.
        center (tuple[float, float]): Center of the ellipse
        color (str): Color of the ellipse.
    """
    # Calculate the eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(cov_matrix)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # Calculate the angle and width/height of the ellipse
    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
    width, height = 2 * np.sqrt(eigvals)

    # Create and add the ellipse patch
    ellipse = Ellipse(center, width, height, angle, edgecolor=color, facecolor='none', linewidth=2)
    ax.add_patch(ellipse)
    ellipse.set_label('Covariance Ellipse')


def plot_dops(dop_ecef: CSVData, dop_enu: CSVData):
    """
    Plots the DOP parameters in three plots:
        * geometry, position and horizontal DOPs
        * ECEF x, y, z and time DOPs
        * East, North and UP DOPs.

    Args:
        dop_ecef (CSVData): input CSVData object with the DOP ECEF dataframe
        dop_enu (CSVData): input CSVData object with the DOP ENU dataframe
    Returns:
        list[matplotlib.pyplot.Axes]: The list with the plotted axes
    """
    sow_array = dop_ecef.data.iloc[:, 1].values
    week_array = dop_ecef.data.iloc[:, 0].values
    time_array = [Epoch.from_gnss_time(week, sow, scale="GPST") for (sow, week) in zip(sow_array, week_array)]

    ax1 = plot_1D(time_array, dop_ecef.data['dop_geometry'].values, label="geometry")
    plot_1D(time_array, dop_ecef.data['dop_position'].values, ax=ax1, label="position")
    plot_1D(time_array, dop_enu.data['dop_horizontal'].values, ax=ax1, x_label="Time", y_label="DOPs [m]",
            label="horizontal", title="Dilution of Precision", set_legend=True)

    ax2 = plot_1D(time_array, dop_ecef.data['dop_x'].values, label="x-axis")
    plot_1D(time_array, dop_ecef.data['dop_y'].values, ax=ax2, label="y-axis")
    plot_1D(time_array, dop_ecef.data['dop_time'].values, ax=ax2, label="time")
    plot_1D(time_array, dop_ecef.data['dop_z'].values, ax=ax2, x_label="Time", y_label="DOPs [m]", label="z-axis",
            title="Dilution of Precision in ECEF", set_legend=True)

    ax3 = plot_1D(time_array, dop_enu.data['dop_east'].values, label="east")
    plot_1D(time_array, dop_enu.data['dop_north'].values, ax=ax3, label="north")
    plot_1D(time_array, dop_enu.data['dop_up'].values, ax=ax3, x_label="Time", y_label="DOPs [m]", label="up",
            title="Dilution of Precision in ENU", set_legend=True)

    return ax1, ax2, ax3


def plot_skyplot(azel: CSVData):
    """
    Plot the skyplot of the satellites in view.

    Args:
        azel (CSVData): input CSVData object with the Satellite Azimuth/Elevation dataframe
    Returns:
        matplotlib.pyplot.Axes: The Axes object with the skyplot
    """
    ax = None
    grouped = azel.data.groupby(['sat'])

    # Add a line for each satellite and data type combination
    for (sat,), group in grouped:
        az_array = group['azimuth'].values
        el_array = group['elevation'].values
        traj = [[az, el] for az, el in zip(az_array, el_array)]
        ax = skyplot.plot(traj, sat, north_to_east_ccw=False, ax=ax)

    ax.set_title("Sky Plot")
    ax.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))

    return ax
