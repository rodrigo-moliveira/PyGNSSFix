from collections import OrderedDict
import numpy

from src.data_types.date import Epoch
from src.data_types.gnss import Satellite
from src.data_mng import TimeSeries
from src.models.gnss_models import compute_nav_sat_eph
from src.io.rinex_parser import SP3OrbitReader
from src.common_log import IO_LOG, get_logger
from src.utils.interpolation import lagrange_interpolation

__all__ = ["SatelliteOrbits"]


class SatelliteOrbits:
    """
    Satellite Orbits DataFrame class
    this class stores and manages satellite orbits, either provided from:
        * navigation message (RINEX NAV)
        * precise products (SP3)
    """

    def __init__(self):
        """Constructor of the `SatelliteOrbits` object"""
        self._data = OrderedDict()  # contains the SP3 data
        self.nav_data = None  # contains the broadcast navigation data
        self.use_precise_products = False
        self.interp_order = 0

    def init(self, nav_data, sp3_files, use_precise_products, interp_order=9, first_epoch=None, last_epoch=None):
        """
        Initialize this SatelliteOrbits instance with navigation and precise data.
        The argument `use_precise_products` allows to select which satellite orbits to output.
        If it is set to True, then the precise SP3 data is used, if set to False the
        navigation orbits are used instead

        Args:
            nav_data(src.data_mng.gnss.navigation_data.NavigationData): the `NavigationData` with navigation data
                (must be already initialized)
            sp3_files(list): list of SP3 files from the user configuration
            use_precise_products(bool): True to use precise orbits and False to use broadcast navigation orbits
            interp_order(int): interpolation order for satellite clocks (default order is 9)
            first_epoch(str): the first epoch to save data from the SP3 files
            last_epoch(str): the last epoch to save data from the SP3 files
        """
        self.nav_data = nav_data
        self.use_precise_products = use_precise_products
        self.interp_order = interp_order

        if use_precise_products:
            log = get_logger(IO_LOG)
            log.info("Launching SP3OrbitReader")
            for file in sp3_files:
                SP3OrbitReader(file, self, first_epoch, last_epoch)

    def __str__(self):
        """Print the SP3 data to a string, for debug purposes"""
        myStr = f"SatelliteOrbits:\n"
        for sat, data in self._data.items():
            myStr += str(sat) + f"\n{str(data)}\n"

        return myStr

    def set_data(self, epoch: Epoch, satellite: Satellite, position: numpy.ndarray):
        """
        method to set a satellite orbit data point for a given epoch and satellite
        Args:
            epoch (Epoch)
            satellite (Satellite)
            position (numpy.ndarray) : satellite orbit from the SP3 Clock file (ECEF cartesian components)
        """

        if not isinstance(epoch, Epoch):
            raise AttributeError(f'First argument should be a valid Epoch object. Type {type(epoch)} was provided '
                                 f'instead')
        if not isinstance(satellite, Satellite):
            raise AttributeError(f'Second argument should be a valid Satellite object. Type {type(satellite)} '
                                 f'was provided instead')
        if not isinstance(position, numpy.ndarray):
            raise AttributeError(f'Third argument should be a numpy.ndarray object. Type {type(position)} '
                                 f'was provided instead')

        if satellite in self._data:
            self._data[satellite].set_data(epoch, position)
        else:
            timeseries = TimeSeries()
            timeseries.set_data(epoch, position)
            self._data[satellite] = timeseries

    # Getters
    def get_data(self):
        return self._data

    def get_sat_data(self, sat):
        try:
            return self._data[sat]
        except KeyError:
            raise KeyError(f"Satellite {str(sat)} has no available navigation data")

    def get_satellites(self):
        """Returns a list with all available satellites (from SP3 file)"""
        return list(self._data.keys())

    def get_epochs(self):
        """Returns a list with all available epochs (from SP3 file)"""
        for key, val in self._data.items():
            return val.get_all_epochs()

    def get_orbit(self, sat, epoch):
        """
        Compute the satellite orbits (position in ECEF frame).
        If `use_precise_products` is True then the SP3 Orbit is returned. Otherwise,
        the clock bias is computed from the provided navigation data.

        Args:
            sat(src.data_types.gnss.satellite.Satellite)
            epoch(src.data_types.date.date.Epoch)
        Returns:
            float: the clock bias for the provided satellite and epoch
        """
        if self.use_precise_products:
            return self.get_orbit_precise(sat, epoch)
        else:
            return self.get_orbit_broadcast(sat, epoch)

    def get_orbit_precise(self, sat, epoch):
        """
        Compute the precise satellite orbits from the provided SP3 Orbit data.
        Applies lagrange interpolation when needed (when the provided epoch is in-between data points).

        Args:
            sat(src.data_types.gnss.satellite.Satellite)
            epoch(src.data_types.date.date.Epoch)
        Returns:
            numpy.ndarray: the position for the provided satellite and epoch

        Raises an exception if the provided satellite does not have valid data
        """
        sat_data = self.get_sat_data(sat)

        # check if the epoch is already contained in the TimeSeries or if an interpolation is needed
        if sat_data.has_epoch(epoch):
            return sat_data.get_data_for_epoch(epoch)
        else:
            # an interpolation is needed
            n_points = (self.interp_order + 1)//2
            knots = sat_data.get_n_items(epoch, n_points)
            poly = lagrange_interpolation([x.mjd for x in knots],
                                          [sat_data.get_data_for_epoch(x) for x in knots],
                                          self.interp_order)

            return poly(epoch.mjd)

    def get_orbit_broadcast(self, sat, epoch):
        """
        Compute the broadcast satellite orbits from the provided RINEX Navigation data.

        Args:
            sat(src.data_types.gnss.satellite.Satellite)
            epoch(src.data_types.date.date.Epoch)
        Returns:
            float: the clock bias for the provided satellite and epoch

        Raises an exception if the provided satellite does not have valid data
        """
        nav_message = self.nav_data.get_closest_message(sat, epoch)
        position, _, _ = compute_nav_sat_eph(nav_message, epoch)
        return position
