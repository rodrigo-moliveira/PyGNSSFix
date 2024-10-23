""" Utility module for interfacing with the matplotlib plot functions. """
import matplotlib.pyplot as plt


def plot_1D(x, y, **kwargs):
    """
    Plots a 1D Time series, where `x` is the time axis and `y` is the data axis.

    Args:
        x(numpy.ndarray or list): array-like vector with the time points
        y(numpy.ndarray or list): array-like vector with the 1D data points
        kwargs: configuration dict with the following available properties:
            ax(matplotlib.pyplot.Axes or None): Axes instance to draw the plot. If None, creates a new one.
            scatter(bool): if True, uses the `scatter` function, otherwise uses the `plot` function
            markersize(float): `markersize` property of `matplotlib`
            marker(str): `marker` property of `matplotlib`
            color(str): `color` property of `matplotlib`
            linewidth(float): `linewidth` property of `matplotlib`
            linestyle(str) `linestyle` property of `matplotlib`
            set_legend(bool): if True, the legend is drawn
            tight_layout(bool): if True, the `tight_layout` property is enabled
            equal(bool): if True, the axes proportions are set to equal
            label(str): label (description) of this data series
            x_label(str): label of the x-axis
            y_label(str): label of the y-axis
            y_scale(str): scale of the y axis
            title(str): title name

    For more information about the configuration properties in `kwargs`, consult the `matplotlib` documentation.
    """
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
        # ax.legend()
        ax.legend(loc='upper right')  # , bbox_to_anchor=(1, 0.5))

    if kwargs.get("tight_layout", False) is True:
        plt.tight_layout()

    if kwargs.get("equal", False) is True:
        ax.axis("equal")

    if "y_scale" in kwargs:
        plt.yscale(kwargs.get("y_scale", "linear"))

    return ax


def loglog(x, y, **kwargs):
    """
    Plots a logarithmic graph, where `x` is the time axis and `y` is the data axis.

    Args:
        x(numpy.ndarray or list): array-like vector with the time points
        y(numpy.ndarray or list): array-like vector with the 1D data points
        kwargs: configuration dict with the following available properties:
            ax(matplotlib.pyplot.Axes or None): Axes instance to draw the plot. If None, creates a new one.
            markersize(float): `markersize` property of `matplotlib`
            marker(str): `marker` property of `matplotlib`
            linewidth(float): `linewidth` property of `matplotlib`
            set_legend(bool): if True, the legend is drawn
            tight_layout(bool): if True, the `tight_layout` property is enabled
            label(str): label (description) of this data series
            x_label(str): label of the x-axis
            y_label(str): label of the y-axis
            title(str): title name

    For more information about the configuration properties in `kwargs`, consult the `matplotlib` documentation.
    """
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
    """ Show all plots in cache. Calls function `matplotlib.pyplot.show` """
    plt.show()
