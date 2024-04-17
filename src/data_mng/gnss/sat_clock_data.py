from collections import OrderedDict

from src.data_types.date import Epoch
from src.data_types.gnss import Satellite
from src.data_mng import TimeSeries
from src.models.gnss_models.clock_obs import *
from src.io.rinex_parser import RinexClockReader
from src.common_log import IO_LOG, get_logger
from src.utils.interpolation import linear_interpolation_scipy

__all__ = ["SatelliteClocks"]


class SatelliteClocks:
    """
    Satellite Clocks DataFrame class
    this class stores and manages satellite clocks, either provided from:
        * navigation message (RINEX NAV)
        * precise products (RINEX CLOCK)

    NavigationData objects contain _data and _header attributes.
    """

    def __init__(self):
        self._data = OrderedDict()
        self.nav_data = None
        self.use_precise_products = False

    def init(self, nav_data, clock_files, use_precise_products, first_obs_epoch, last_obs_epoch):
        self.nav_data = nav_data
        self.use_precise_products = use_precise_products

        if use_precise_products:
            log = get_logger(IO_LOG)
            log.info("Launching RinexClockReader")
            for file in clock_files:
                RinexClockReader(file, self, first_obs_epoch, last_obs_epoch)

    def __str__(self):
        """Print the navigation data to a string, for debug purposes"""
        myStr = f"SatelliteClocks:\n"
        for sat, data in self._data.items():
            myStr += str(sat) + f"\n{str(data)}\n"

        return myStr

    def set_data(self, epoch: Epoch, satellite: Satellite, clock: float):
        """
        TODO: update docstring
        method to set a navigation data point for a given epoch and satellite
        Args:
            epoch (Epoch)
            satellite (Satellite)
            clock (float)
        """

        if not isinstance(epoch, Epoch):
            raise AttributeError(f'First argument should be a valid Epoch object. Type {type(epoch)} was provided '
                                 f'instead')
        if not isinstance(satellite, Satellite):
            raise AttributeError(f'Second argument should be a valid Satellite object. Type {type(satellite)} '
                                 f'was provided instead')
        if not isinstance(clock, float):
            raise AttributeError(f'Third argument should be a float object. Type {type(clock)} '
                                 f'was provided instead')

        if satellite in self._data:
            self._data[satellite].set_data(epoch, clock)
        else:
            timeseries = TimeSeries()
            timeseries.set_data(epoch, clock)
            self._data[satellite] = timeseries

    # Getters
    def get_data(self):
        return self._data

    def get_sat_data(self, sat):
        try:
            return self._data[sat]
        except KeyError:
            raise KeyError(f"Satellite {str(sat)} has no available navigation data")

    def get_clock(self, sat, epoch):
        if self.use_precise_products:
            return self.get_clock_precise(sat, epoch)
        else:
            pass  # return broadcast clock
    def get_clock_precise(self, sat, epoch):
        # Apply linear interpolation. Get node before and after the provided epoch
        # see https://en.wikipedia.org/wiki/Linear_interpolation
        # apply interpolation if precise clock or navigation broadcast function otherwise.
        # try/catch here... return None or exception if anything fails...
        # log to the logger.
        # add to docstrings: raises Exception if problem is found. Try catch must be on the **caller** side

        sat_data = self.get_sat_data(sat)

        # check if the epoch is already contained in the TimeSeries or if an interpolation is needed
        if sat_data.has_epoch(epoch):
            return sat_data.get_data_for_epoch(epoch)
        else:
            # an interpolation is needed
            knots = sat_data.get_epoch_knots(epoch)
            return linear_interpolation_scipy(epoch.mjd, [x.mjd for x in knots],
                                              [sat_data.get_data_for_epoch(x) for x in knots])
