import re

import allantools as at

from src.data_mng.csv.csv_data import CSVData
from src.data_types.date import Epoch
from src.io.rinex_parser.utils import RINEX_SATELLITE_SYSTEM, RINEX_OBS_TYPES_UNITS
from .utils import *


def plot_observations(observations: CSVData):
    """
    Plots the GNSS Observables

    Arguments:
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

    Arguments:
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
    Plot the GNSS prefit and postfit residuals of the estimation process

    Arguments:
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
