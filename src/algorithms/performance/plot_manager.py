import matplotlib.pyplot as plt

from PositioningSolver.src.gnss.data_types.ObservationData import ObservationData
from ..data_types.basics.Epoch import Epoch
from ..data_types.containers.TimeSeries import TimeSeries
from ..plots.skyplot import plot_sky

from ..utils.errors import NonExistentObservable


def plot_observables(observation_data: ObservationData, satellite, datatype, **kwargs):
    epochs = observation_data.get_epochs()

    times = []
    observables = []

    for epoch in epochs:
        try:
            value = observation_data.get_observable_at_epoch(satellite, epoch, datatype)
            times.append(epoch.to_datetime())
            observables.append(float(value))

        except NonExistentObservable:
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


def plot_skyplot(sat_info_series):
    ax = None

    _data = {}
    for epoch, sat_info in sat_info_series.items():
        for sat, info in sat_info.items():

            if str(sat) not in _data:
                _data[str(sat)] = []

            _data[str(sat)].append([info.az, info.el])

    for sat, data in _data.items():
        ax = plot_sky(data, sat, north_to_east_ccw=False, style_kwargs={'s': 10}, ax=ax)

    ax.set_title("Sky Plot")
    ax.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))


def plot_3D_trajectory(data_points, **kwargs):
    # try to fetch axis to insert the plot. If no ax is provided, create a new one
    ax = kwargs.get("ax", None)
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

    ax.plot([x[0] for x in data_points], [x[1] for x in data_points], [x[2] for x in data_points],
            label=kwargs.get("label", ""))
    if "true_position" in kwargs:
        true = kwargs.get("true_position")
        ax.scatter(true[0], true[1], true[2], marker='*', label="True State")

    ax.ticklabel_format(useOffset=False)

    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_zlabel(kwargs.get("z_label", ""))
    ax.set_title(kwargs.get("title", ""))

    if "label" in kwargs:
        ax.legend()
    return ax


def plot_1D_TimeSeries(series: TimeSeries, **kwargs):
    time, data = series.export2time_data()

    # convert x from Epoch to datetime Objects
    if isinstance(time[0], Epoch):
        time = [i.to_datetime() for i in time]

    # try to fetch axis to insert the plot. If no ax is provided, create a new one
    ax = kwargs.get("ax", None)
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(time, data, linewidth=2.0)
    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_title(kwargs.get("title", ""))

    return ax


def grid():
    plt.grid(True)


def plot_1D(x, y, **kwargs):
    # convert x from Epoch to datetime Objects
    if isinstance(x[0], Epoch):
        x = [i.to_datetime() for i in x]

    # try to fetch axis to insert the plot. If no ax is provided, create a new one
    ax = kwargs.get("ax", None)
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(x, y, linewidth=kwargs.get("linewidth", 2.0), label=kwargs.get("label", None),
            marker=kwargs.get("marker", ''), markersize=kwargs.get("markersize", 5))
    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_title(kwargs.get("title", ""))

    if "set_legend" in kwargs:
        ax.legend()

    if kwargs.get("tight_layout", False) is True:
        plt.tight_layout()

    if kwargs.get("equal", False) is True:
        ax.axis("equal")

    if "y_scale" in kwargs:
        plt.yscale(kwargs.get("y_scale", "linear"))

    return ax


def loglog(x, y, **kwargs):
    # convert x from Epoch to datetime Objects
    if isinstance(x[0], Epoch):
        x = [i.to_datetime() for i in x]

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


def plot_satellite_availability(sat_info: TimeSeries, **kwargs):
    from matplotlib.ticker import MaxNLocator

    x = list(sat_info.keys())
    y = list(sat_info.values())
    availability = [len(_y) for _y in y]

    # convert x from Epoch to datetime Objects
    if isinstance(x[0], Epoch):
        x = [i.to_datetime() for i in x]

    # plot
    fig, ax = plt.subplots()

    ax.plot(x, availability, linewidth=2.0)
    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_title(kwargs.get("title", ""))

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))


def plot_2D_trajectory(x, y, **kwargs):
    # try to fetch axis to insert the plot. If no ax is provided, create a new one
    ax = kwargs.get("ax", None)
    if ax is None:
        fig, ax = plt.subplots()

    # plot
    ax.scatter(x, y, linewidth=0.2, marker=kwargs.get("marker", 'o'), label=kwargs.get("label", None))

    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_title(kwargs.get("title", ""))

    if "set_legend" in kwargs:
        ax.legend()

    return ax


def show_all():
    plt.show()
