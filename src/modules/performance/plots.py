import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import Axes3D
from typing import List, Tuple, Optional

from src.modules.performance.skyplot import plot_sky
from src.data_types.date import Epoch


def format_time_for_plot(time_in):
    time_out = []
    for t in time_in["Epoch_timetag"]:
        args = t.split()
        epoch = Epoch.strptime(args[0], format=Epoch.ISO_FORMAT, scale=args[1])
        time_out.append(epoch)
    return time_out


# TODO: fix kwargs as below
def plot_observables(observation_data, satellite, datatype, **kwargs):
    epochs = observation_data.get_epochs()

    times = []
    observables = []

    for epoch in epochs:
        try:
            value = observation_data.get_observable_at_epoch(satellite, epoch, datatype)
            times.append(epoch.to_datetime())
            observables.append(float(value))

        except:
            continue

    # plot
    if "ax" not in kwargs:
        fig, ax = plt.subplots()
    else:
        ax = kwargs["ax"]

    ax.scatter(times, observables, s=1.0)
    ax.set_xlabel(kwargs.get("x_label", "Time"))
    ax.set_ylabel(kwargs.get("y_label", str(datatype)))
    ax.set_title(kwargs.get("title", ""))

    return ax


def plot_skyplot(azel):
    ax = None

    _data = {}
    for epoch, sat_info in azel.items():
        for sat, info in sat_info.items():

            if str(sat) not in _data:
                _data[str(sat)] = []
            _data[str(sat)].append(info)

    for sat, data in _data.items():
        ax = plot_sky(data, sat, north_to_east_ccw=False, style_kwargs={'s': 10}, ax=ax)

    ax.set_title("Sky Plot")
    ax.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))


def plot_covariance_ellipsoid(ax: Axes3D, pos: Tuple[float, float, float], cov: np.ndarray, n_std: float = 1.0,
                              **kwargs):
    """
    Plot a 3D covariance ellipsoid.

    Parameters:
    - ax (Axes3D): Axes on which to plot.
    - pos (Tuple[float, float, float]): Center position of the ellipsoid.
    - cov (np.ndarray): Covariance matrix (3x3).
    - n_std (float): Number of standard deviations to scale the ellipsoid.
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
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='c', alpha=0.1, edgecolor='k', linewidth=0.5, **kwargs)


def plot_covariance_ellipsoid_diagonal(ax: Axes3D, pos: Tuple[float, float, float],
                                       variances: Tuple[float, float, float], n_std: float = 1.0, **kwargs):
    """
    Plot a 3D covariance ellipsoid using only variances (diagonal covariance matrix).

    Parameters:
    - ax (Axes3D): Axes on which to plot.
    - pos (Tuple[float, float, float]): Center position of the ellipsoid.
    - variances (Tuple[float, float, float]): Variances along x, y, and z axes.
    - n_std (float): Number of standard deviations to scale the ellipsoid.
    """
    radii = n_std * np.sqrt(variances)

    # Generate points on a unit sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 50)
    x = np.outer(np.cos(u), np.sin(v)) * radii[0]
    y = np.outer(np.sin(u), np.sin(v)) * radii[1]
    z = np.outer(np.ones_like(u), np.cos(v)) * radii[2]

    # Shift the ellipsoid to the correct position
    x += pos[0]
    y += pos[1]
    z += pos[2]

    # Plot the surface with reduced visibility
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='c', alpha=0.01,
                    edgecolor='lightblue', linewidth=0.5, **kwargs)


def plot_3D_trajectory_with_avg_covariance(
        data_points: List[Tuple[float, float, float]],
        cov: List, # TODO: add here the proper list float arguments
        true_position: Optional[Tuple[float, float, float]] = None,
        x_label: str = "X-axis", y_label: str = "Y-axis", z_label: str = "Z-axis",
        title: str = "3D Trajectory with Covariance", est_marker: str = 'o', true_marker: str = '*',
        mean_marker='D',  est_color: str = 'k', true_color: str = 'r', mean_color='b',
        est_size: int = 20, true_size: int = 80, mean_size=80, ax: Optional[Axes3D] = None) -> None:
    """
    Plots the 3D trajectory of estimated positions with an average covariance ellipsoid against the true static position.

    Parameters:
    - data_points (List[Tuple[float, float, float, float, float, float, float, float, float]]): List of estimated positions and their covariances.
    - true_position (Optional[Tuple[float, float, float]]): True static position.
    - x_label (str): Label for the X-axis.
    - y_label (str): Label for the Y-axis.
    - z_label (str): Label for the Z-axis.
    - title (str): Title of the plot.
    - est_marker (str): Marker style for estimated positions.
    - true_marker (str): Marker style for the true position.
    - est_color (str): Color for estimated positions.
    - true_color (str): Color for the true position.
    - est_size (int): Marker size for estimated positions.
    - true_size (int): Marker size for the true position.
    - ax (Optional[Axes3D]): Pre-existing axes for the plot.
    """
    # Validate data_points
    #if not all(len(point) < 3 for point in data_points):
    #    raise ValueError(
    #        "Each data point must have 3 elements "
    #        "(pos_x, pos_y, pos_z).")

    # Validate true_position
    if true_position and len(true_position) != 3:
        raise ValueError("True position must be a 3-dimensional point.")

    # If no axis is provided, create a new one
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

    # Plot estimated positions
    positions = np.array([[x[0], x[1], x[2]] for x in data_points])
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
               label="Estimated Position", marker=est_marker, color=est_color, s=est_size)

    # Calculate the average covariance matrix
    mean_cov_matrix = np.mean(np.array([np.array(c) for c in cov]), axis=0)

    # Plot the average covariance ellipsoid at the mean position
    mean_position = np.mean(positions, axis=0)
    plot_covariance_ellipsoid(ax, mean_position, mean_cov_matrix)

    # Plot mean estimated position
    ax.scatter(mean_position[0], mean_position[1], mean_position[2],
               label="Mean Estimated Position", s=mean_size, marker=mean_marker, color=mean_color)

    # Plot true position if provided
    if true_position:
        ax.scatter(true_position[0], true_position[1], true_position[2],
                   label="True Static Position", marker=true_marker, color=true_color, s=true_size)

    # Set labels and title
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    ax.set_title(title)

    # Format the ticks
    ax.ticklabel_format(useOffset=False)

    # Show legend
    ax.legend()


def grid():
    plt.grid(True)


# TODO: fix kwargs as in the other
def plot_1D(x, y, **kwargs):
    # try to fetch axis to insert the plot. If no ax is provided, create a new one
    ax = kwargs.get("ax", None)
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(x, y, linewidth=kwargs.get("linewidth", 2.0), label=kwargs.get("label", None),
            marker=kwargs.get("marker", ''), markersize=kwargs.get("markersize", 5),
            linestyle=kwargs.get("linestyle", "solid"), color=kwargs.get('color', None))
    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_title(kwargs.get("title", ""))

    if "set_legend" in kwargs:
        #ax.legend()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if kwargs.get("tight_layout", False) is True:
        plt.tight_layout()

    if kwargs.get("equal", False) is True:
        ax.axis("equal")

    if "y_scale" in kwargs:
        plt.yscale(kwargs.get("y_scale", "linear"))

    return ax


def loglog(x, y, **kwargs):
    # try to fetch axis to insert the plot. If no ax is provided, create a new one
    ax = kwargs.get("ax", None)
    if ax is None:
        fig, ax = plt.subplots()

    ax.loglog(x, y, linewidth=kwargs.get("linewidth", 2.0), label=kwargs.get("label", None),
              marker=kwargs.get("marker", ''), markersize=kwargs.get("markersize", 5))
    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_title(kwargs.get("title", ""))

    if "set_legend" in kwargs:
        ax.legend()

    if kwargs.get("tight_layout", False) is True:
        plt.tight_layout()

    return ax


def plot_satellite_availability(sat_info, **kwargs):
    from matplotlib.ticker import MaxNLocator
    time = []
    sats = []

    for epoch, data in sat_info.items():
        time.append(float(epoch[1]))
        sats.append(list(data.keys()))

    availability = [len(_y) for _y in sats]

    # plot
    fig, ax = plt.subplots()

    ax.plot(time, availability, linewidth=2.0)
    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_title(kwargs.get("title", ""))

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))


def plot_2D_trajectory(data, cov, **kwargs):
    """
    Plots the 2D trajectory of estimated positions with an average covariance ellipse.

    Parameters:
    - data (list of tuples): List of estimated positions in the ENU frame.
    - cov (pandas.Series): Series of covariance matrices.
    - kwargs: Additional keyword arguments for customization.
      - ax (matplotlib.axes.Axes): Pre-existing axes for the plot.
      - marker (str): Marker style for estimated positions.
      - label (str): Label for the estimated positions.
      - x_label (str): Label for the X-axis.
      - y_label (str): Label for the Y-axis.
      - title (str): Title of the plot.

    Returns:
    - ax (matplotlib.axes.Axes): The axes with the plot.
    """
    # Try to fetch axis to insert the plot. If no ax is provided, create a new one
    ax = kwargs.get("ax", None)
    if ax is None:
        fig, ax = plt.subplots()

    x = [d[0] for d in data]
    y = [d[1] for d in data]

    # Plot the estimated positions
    ax.scatter(x, y, linewidth=0.2, marker=kwargs.get("marker", 'o'), label=kwargs.get("label", "Estimated Positions"))
    ax.scatter(0, 0, linewidth=0.2, marker='*', label="True Position")

    # Calculate the mean covariance matrix
    mean_cov_matrix = np.mean(np.array([np.array(c) for c in cov]), axis=0)

    # Plot the covariance ellipse
    plot_covariance_ellipse(ax, mean_cov_matrix, kwargs.get("ellipse_color", 'r'))

    # Set labels and title
    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_title(kwargs.get("title", ""))

    # Show legend
    ax.legend(loc='upper right')

    return ax


def plot_covariance_ellipse(ax, cov_matrix, color='r'):
    """
    Plots a 2D covariance ellipse based on the covariance matrix.

    Parameters:
    - ax (matplotlib.axes.Axes): The axes to plot on.
    - cov_matrix (numpy.ndarray): The 3x3 covariance matrix.
    - color (str): Color of the ellipse.
    """
    # Extract the 2x2 covariance matrix for the EN plane (East-North)
    cov_2d = cov_matrix[:2, :2]

    # Calculate the eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eigh(cov_2d)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # Calculate the angle and width/height of the ellipse
    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
    width, height = 2 * np.sqrt(eigvals)

    # Create and add the ellipse patch
    ellipse = Ellipse((0, 0), width, height, angle, edgecolor=color, facecolor='none', linewidth=2)
    ax.add_patch(ellipse)
    ellipse.set_label('Covariance Ellipse')


def show_all():
    plt.show()
