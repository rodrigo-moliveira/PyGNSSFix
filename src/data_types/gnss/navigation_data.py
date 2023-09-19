from collections import OrderedDict
from src.errors import TimeSeriesError
from src.data_types.date.date import Epoch
from src.data_types.gnss.satellite import Satellite
from src.data_mng.timeseries import TimeSeries
from src.data_mng.container import Container


class NavigationHeader(Container):
    """
    NavigationHeader class, inherits from Container
    stores relevant data from the header section of a navigation file
    """
    __slots__ = ["rinex_version", "satellite_system",
                 "iono_corrections", "leap_seconds", "first_epoch"]

    def __init__(self):
        super().__init__()
        for attr in self.__slots__:
            setattr(self, attr, None)
        self.iono_corrections = {}


class NavigationPoint(Container):
    """
    NavigationPointGPS class, inherits from Container
    stores the data contained in a single navigation message for GPS satellites
    """
    __slots__ = ["satellite", "toc", "af0", "af1", "af2",
                 "IODE", "crs", "deltaN", "M0",
                 "cuc", "eccentricity", "cus", "sqrtA",
                 "cic", "RAAN0", "cis",
                 "i0", "crc", "omega", "RAANDot",
                 "iDot", "codesL2", "toe", "flagL2",
                 "SV_URA", "SV_health", "TGD", "IODC",
                 "TransmissionTime"]

    def __init__(self):
        super().__init__()
        for attr in self.__slots__:
            setattr(self, attr, None)

    def __str__(self):
        _allAttrs = ""
        for atr in self.__slots__:
            if atr == "satellite" or atr == "toc":
                continue
            _allAttrs += atr + "=" + str(getattr(self, atr)) + ", "
        _allAttrs = _allAttrs[0:-2]
        return f'{type(self).__name__}({_allAttrs})'


class NavigationData:
    """
    NavigationDataMap
    this class stores data from rinex navigation files
    """

    def __init__(self):
        self._data = OrderedDict()
        self._header = NavigationHeader()

    def __str__(self):

        myStr = f"{repr(self._header)}\n"
        for sat, data in self._data.items():
            myStr += str(sat) + f"\n{str(data)}\n"

        return myStr

    def set_data(self, epoch: Epoch, satellite: Satellite, nav_message: NavigationPoint):
        """
        method to set a navigation data point for a given epoch and satellite
        Args:
            epoch (Epoch)
            satellite (Satellite)
            nav_message (NavigationPoint)
        """

        if not isinstance(epoch, Epoch):
            raise AttributeError(f'First argument should be a valid Epoch object. Type {type(epoch)} was provided '
                                 f'instead')
        if not isinstance(satellite, Satellite):
            raise AttributeError(f'Second argument should be a valid Satellite object. Type {type(satellite)} '
                                 f'was provided instead')
        if not issubclass(type(nav_message), NavigationPoint):
            raise AttributeError(f'Third argument should be a valid NavigationData object. Type {type(nav_message)} '
                                 f'was provided instead')

        if satellite in self._data:
            self._data[satellite].set_data(epoch, nav_message)
        else:
            timeseries = TimeSeries()
            timeseries.set_data(epoch, nav_message)
            self._data[satellite] = timeseries

    def set_header(self, header: NavigationHeader):
        """
        method to set the navigation header.
        Note: the header is valid from the first epoch to the last. Multiple files can be used to increase the cover
        Args:
            header (NavigationHeader)
        """
        if not issubclass(type(header), NavigationHeader):
            raise AttributeError(f'Third argument should be a valid NavigationData object. Type {type(header)} '
                                 f'was provided instead')

        self._header = header

    # Getters
    def get_data(self):
        return self._data

    def get_sat_data(self, sat):
        try:
            return self._data[sat]
        except KeyError:
            raise KeyError(f"Satellite {str(sat)} has no available navigation data")

    def get_sat_data_for_epoch(self, sat, epoch):
        """
        gets navigation message closest to the provided epoch

        Args:
            sat (Satellite)
            epoch (Epoch)
        Return:
            NavigationPointGPS: navigation data for the given satellite closest to the provided epoch
        Raises:
            TimeSeriesError
        """
        try:
            _epoch = self._data[sat].get_closest_epoch(epoch)
            return self._data[sat].get_data_for_epoch(_epoch)
        except TimeSeriesError as e:
            raise TimeSeriesError(f"satellite {str(sat)} has no available navigation data for epoch {str(epoch)}, {e}")

    @property
    def header(self):
        return self._header
