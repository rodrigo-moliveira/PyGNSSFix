""" Satellite Clocks Module
This module implements the Satellite Clocks manager, that manages the computation of satellite clock corrections in the
reconstruction of GNSS observations.
"""
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
    Satellite Clocks Manager class.
    This class stores and manages satellite clocks, either provided from:
        * navigation message (RINEX Navigation files)
        * precise products (RINEX Clock files)

    Attributes:
        precise_data(OrderedDict): attribute that stores the precise clock products. It is an :py:class:`OrderedDict`
            defined as:
                * keys -> :py:class:`Satellite` instances.
                * values -> :py:class:`TimeSeries` instances with the clock products.
        nav_data(src.data_mng.gnss.navigation_data.NavigationData): attribute that stores the broadcast navigation data
        use_precise_products(bool): True to use precise clocks and False to use broadcast clocks
        interp_order(int): interpolation order for satellite clocks
    """

    def __init__(self):
        """ Constructor of the `SatelliteClocks` instance """
        self.precise_data = OrderedDict()  # contains the RINEX Clock data
        self.nav_data = None  # contains the broadcast navigation data
        self.use_precise_products = False
        self.interp_order = 0

    def init(self, nav_data, clock_files, use_precise_products, interp_order=1):
        """
        Initialize this SatelliteClocks instance with navigation and precise data.
        The argument `use_precise_products` allows to select which satellite clocks to output.
        If it is set to True, then the RINEX Clocks are used, if set to False the navigation clocks
        are used instead.

        Args:
            nav_data(src.data_mng.gnss.navigation_data.NavigationData): the `NavigationData` with navigation data
                (must be already initialized)
            clock_files(list): list of RINEX clock files from the user configuration
            use_precise_products(bool): True to use precise clocks and False to use broadcast navigation clocks
            interp_order(int): interpolation order for satellite clocks (default order is 1)
        """
        self.nav_data = nav_data
        self.use_precise_products = use_precise_products
        self.interp_order = interp_order  # default interpolation order for clocks is 1

        if use_precise_products:
            log = get_logger(IO_LOG)
            log.info("Launching RinexClockReader")
            for file in clock_files:
                RinexClockReader(file, self)

    def __str__(self):
        """ Print the RINEX clock data to a string, for debug purposes """
        myStr = f"SatelliteClocks:\n"
        for sat, data in self.precise_data.items():
            myStr += str(sat) + f"\n{str(data)}\n"

        return myStr

    def set_data(self, epoch: Epoch, satellite: Satellite, clock: float):
        """
        Method to set a satellite precise clock correction for a given epoch and satellite

        Args:
            epoch (Epoch): epoch to set the correction
            satellite (Satellite): satellite to set the correction
            clock (float) : value of the precise clock bias correction [s].

        Raises:
            AttributeError: this exception is raised if one of the input arguments is not of the valid type.
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

        if satellite in self.precise_data:
            self.precise_data[satellite].set_data(epoch, clock)
        else:
            timeseries = TimeSeries()
            timeseries.set_data(epoch, clock)
            self.precise_data[satellite] = timeseries

    # Getters
    def get_data(self) -> OrderedDict:
        """ Returns the precise clock products. """
        return self.precise_data

    def get_sat_data(self, sat) -> TimeSeries:
        """
        Returns:
            TimeSeries: returns the precise clock products (saved in a TimeSeries) for the queried satellite

        Raises:
              KeyError: an exception is raised if the satellite is not available.
        """
        try:
            return self.precise_data[sat]
        except KeyError:
            raise KeyError(f"The loaded precise clock products do not have data for satellite {str(sat)}")

    def get_clock(self, sat, epoch):
        """
        Compute the satellite clock and drift.
        If the attribute `use_precise_products` is True then the clocks are computed from the precise products.
        Otherwise, they are computed from the broadcast navigation data.

        Args:
            sat(src.data_types.gnss.Satellite): satellite to compute the clock bias and drift
            epoch(src.data_types.date.Epoch): epoch for the computation
        Returns:
            tuple[float, float]: the clock bias (sec) and drift (sec/sec) for the provided satellite and epoch
        """
        if self.use_precise_products:
            # NOTE: for now, no computation of drift for precise data is performed
            return self.get_clock_precise(sat, epoch), 0.0
        else:
            return self.get_clock_broadcast(sat, epoch)

    def get_satellites(self) -> list[Satellite]:
        """ Returns a list with all available satellites (defined in the precise products) """
        return list(self.precise_data.keys())

    def get_epochs(self):
        """ Returns a list with all available epochs (defined in the precise products) """
        for key, val in self.precise_data.items():
            return val.get_all_epochs()

    def get_clock_precise(self, sat, epoch):
        """
        Compute the precise satellite clocks from the provided RINEX Clock data.
        Applies linear interpolation when needed (when the queried epoch is in-between data points).

        Args:
            sat(src.data_types.gnss.Satellite): satellite for the computation
            epoch(src.data_types.date.Epoch): Epoch for the computation

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
            n_points = (self.interp_order + 1) // 2
            knots = sat_data.get_n_items(epoch, n_points)

            return linear_interpolation_scipy((epoch - knots[0]).total_seconds(),
                                              [(x - knots[0]).total_seconds() for x in knots],
                                              [sat_data.get_data_for_epoch(x) for x in knots])

    def get_clock_broadcast(self, sat, epoch):
        """
        Compute the broadcast satellite clocks from the provided RINEX Navigation data.

        Args:
            sat(src.data_types.gnss.Satellite): satellite for the computation
            epoch(src.data_types.date.Epoch): Epoch for the computation
        Returns:
            tuple[float, float]: the clock bias and clock drift for the provided satellite and epoch

        Raises an exception if the provided satellite does not have valid data
        """
        nav_message = self.nav_data.get_closest_message(sat, epoch)
        dt, drift = broadcast_clock(nav_message.af0, nav_message.af1, nav_message.af2,
                                    nav_message.toc.gnss_time[1], epoch.gnss_time[1])
        return dt, drift

    def get_nav_message(self, sat, epoch):
        """
        Fetches the closest navigation message valid for the satellite and epoch.

        Args:
            sat(src.data_types.gnss.Satellite): satellite to fetch the data
            epoch(src.data_types.date.Epoch): epoch to fetch the data
        Returns:
            src.data_mng.gnss.navigation_data.NavigationPoint: navigation message

        Raises an exception if the provided satellite does not have valid data
        """
        return self.nav_data.get_closest_message(sat, epoch)
