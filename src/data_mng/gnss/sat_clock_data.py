from collections import OrderedDict

from src.data_types.date import Epoch
from src.data_types.gnss import Satellite
from src.data_mng import TimeSeries
from src.models.gnss_models import broadcast_clock
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
    """

    def __init__(self):
        """Constructor of the `SatelliteClocks` object"""
        self._data = OrderedDict()  # contains the RINEX Clock data
        self.nav_data = None  # contains the broadcast navigation data
        self.use_precise_products = False

    def init(self, nav_data, clock_files, use_precise_products, first_epoch=None, last_epoch=None):
        """
        Initialize this SatelliteClocks instance with navigation and precise data.
        The argument `use_precise_products` allows to select which satellite clocks to output.
        If it is set to True, then the RINEX Clocks are used, if set to False the navigation clocks
        are used instead

        Args:
            nav_data(src.data_mng.gnss.navigation_data.NavigationData): the `NavigationData` with navigation data
                (must be already initialized)
            clock_files(list): list of RINEX clock files from the user configuration
            use_precise_products(bool): True to use precise clocks and False to use broadcast navigation clocks
            first_epoch(str): the first epoch to store RINEX clocks
            last_epoch(str): the last epoch to store RINEX clocks
        """
        self.nav_data = nav_data
        self.use_precise_products = use_precise_products

        if use_precise_products:
            log = get_logger(IO_LOG)
            log.info("Launching RinexClockReader")
            for file in clock_files:
                RinexClockReader(file, self, first_epoch, last_epoch)

    def __str__(self):
        """Print the RINEX clock data to a string, for debug purposes"""
        myStr = f"SatelliteClocks:\n"
        for sat, data in self._data.items():
            myStr += str(sat) + f"\n{str(data)}\n"

        return myStr

    def set_data(self, epoch: Epoch, satellite: Satellite, clock: float):
        """
        method to set a satellite clock data point for a given epoch and satellite
        Args:
            epoch (Epoch)
            satellite (Satellite)
            clock (float) : satellite clock bias from the RINEX Clock file
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
        """
        Compute the satellite clock.
        If `use_precise_products` is True then the RINEX Clock is returned. Otherwise,
        the clock bias is computed from the provided navigation data.

        Args:
            sat(src.data_types.gnss.satellite.Satellite)
            epoch(src.data_types.date.date.Epoch)
        Returns:
            float: the clock bias for the provided satellite and epoch
        """
        if self.use_precise_products:
            return self.get_clock_precise(sat, epoch)
        else:
            return self.get_clock_broadcast(sat, epoch)

    def get_satellites(self):
        """Returns a list with all available satellites (RINEX clock)"""
        return list(self._data.keys())

    def get_epochs(self):
        """Returns a list with all available epochs (RINEX clock)"""
        for key, val in self._data.items():
            return val.get_all_epochs()

    def get_clock_precise(self, sat, epoch):
        """
        Compute the precise satellite clocks from the provided RINEX Clock data.
        Applies linear interpolation when needed (when the provided epoch is in-between data points).

        Args:
            sat(src.data_types.gnss.satellite.Satellite)
            epoch(src.data_types.date.date.Epoch)
        Returns:
            float: the clock bias for the provided satellite and epoch

        Raises an exception if the provided satellite does not have valid data
        """
        sat_data = self.get_sat_data(sat)

        # check if the epoch is already contained in the TimeSeries or if an interpolation is needed
        if sat_data.has_epoch(epoch):
            return sat_data.get_data_for_epoch(epoch)
        else:
            # an interpolation is needed
            knots = sat_data.get_epoch_knots(epoch)
            return linear_interpolation_scipy(epoch.mjd, [x.mjd for x in knots],
                                              [sat_data.get_data_for_epoch(x) for x in knots])

    def get_clock_broadcast(self, sat, epoch):
        """
        Compute the broadcast satellite clocks from the provided RINEX Navigation data.

        Args:
            sat(src.data_types.gnss.satellite.Satellite)
            epoch(src.data_types.date.date.Epoch)
        Returns:
            float: the clock bias for the provided satellite and epoch

        Raises an exception if the provided satellite does not have valid data
        """
        nav_message = self.nav_data.get_closest_message(sat, epoch)
        dt_sat, _ = broadcast_clock(nav_message.af0, nav_message.af1, nav_message.af2,
                                    nav_message.toc.gnss_time[1], epoch.gnss_time[1])
        return dt_sat
