import matplotlib.pyplot as plt
from src.data_types.date import Epoch

def plot_1D(x, y, **kwargs):
    # try to fetch axis to insert the plot. If no ax is provided, create a new one
    ax = kwargs.get("ax", None)
    if ax is None:
        fig, ax = plt.subplots()

    if kwargs.get("scatter", False):
        ax.scatter(x, y, label=kwargs.get("label", None), s=kwargs.get("markersize", 5),
                color=kwargs.get('color', None))
    else:
        ax.plot(x, y, linewidth=kwargs.get("linewidth", 2.0), label=kwargs.get("label", None),
                marker=kwargs.get("marker", ''), markersize=kwargs.get("markersize", 5),
                linestyle=kwargs.get("linestyle", "solid"), color=kwargs.get('color', None))
    ax.set_xlabel(kwargs.get("x_label", ""))
    ax.set_ylabel(kwargs.get("y_label", ""))
    ax.set_title(kwargs.get("title", ""))

    if "set_legend" in kwargs:
        #ax.legend()
        ax.legend(loc='upper right')#, bbox_to_anchor=(1, 0.5))

    if kwargs.get("tight_layout", False) is True:
        plt.tight_layout()

    if kwargs.get("equal", False) is True:
        ax.axis("equal")

    if "y_scale" in kwargs:
        plt.yscale(kwargs.get("y_scale", "linear"))

    return ax


# USED FUNCTIONS ABOVE


def format_time_for_plot(time_in):
    time_out = []
    for t in time_in["Epoch_timetag"]:
        args = t.split()
        epoch = Epoch.strptime(args[0], format=Epoch.ISO_FORMAT, scale=args[1])
        time_out.append(epoch)
    return time_out


# TODO: fix kwargs as below
def plot_observables(observation_data, satellite, datatype, **kwargs):
    # TODO: to be deleted
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
def grid():
    plt.grid(True)


# TODO: fix kwargs as in the other



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
    if kwargs.get("grid", False):
        ax.grid(True)

    if "set_legend" in kwargs:
        ax.legend()

    if kwargs.get("tight_layout", False) is True:
        plt.tight_layout()

    return ax




def show_all():
    plt.show()




