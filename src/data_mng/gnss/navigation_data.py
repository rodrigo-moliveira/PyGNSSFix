from collections import OrderedDict
from src.errors import TimeSeriesError, EphemerideError
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
                 "iono_corrections", "leap_seconds", "time_correction"]

    def __init__(self):
        super().__init__()
        for attr in self.__slots__:
            setattr(self, attr, None)
        self.iono_corrections = {}
        self.time_correction = {}


class NavigationPoint(Container):

    def __int__(self):
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

    @property
    def constellation(self):
        return None


class NavigationPointGPS(NavigationPoint):
    """
    NavigationPointGPS class, inherits from Container
    stores the data contained in a single navigation message for GPS satellites
    """
    __slots__ = ["satellite", "toc", "af0", "af1", "af2",       # sv / Toc / sv clk
                 "IODE", "crs", "deltaN", "M0",                 # BROADCAST ORBIT 1
                 "cuc", "eccentricity", "cus", "sqrtA",         # BROADCAST ORBIT 2
                 "cic", "RAAN0", "cis",                         # BROADCAST ORBIT 3
                 "i0", "crc", "omega", "RAANDot",               # BROADCAST ORBIT 4
                 "iDot", "codesL2", "toe", "flagL2",            # BROADCAST ORBIT 5
                 "SV_URA", "SV_health", "TGD", "IODC",          # BROADCAST ORBIT 6
                 "TransmissionTime"]                            # BROADCAST ORBIT 7

    @property
    def constellation(self):
        return "GPS"


class NavigationPointGAL(NavigationPoint):
    """
    NavigationPointGPS class, inherits from Container
    stores the data contained in a single navigation message for GAL satellites
    """
    __slots__ = ["satellite", "toc", "af0", "af1", "af2",       # sv / Toc / sv clk
                 "IODnav", "crs", "deltaN", "M0",               # BROADCAST ORBIT 1
                 "cuc", "eccentricity", "cus", "sqrtA",         # BROADCAST ORBIT 2
                 "cic", "RAAN0", "cis",                         # BROADCAST ORBIT 3
                 "i0", "crc", "omega", "RAANDot",               # BROADCAST ORBIT 4
                 "iDot", "dataSrc", "toe",                      # BROADCAST ORBIT 5
                 "SISA", "SV_health", "BGDE1E5a", "BGDE1E5b",       # BROADCAST ORBIT 6
                 "TransmissionTime",                            # BROADCAST ORBIT 7
                 "nav_type"]

    def __init__(self):
        super().__init__()
        self.nav_type = None

    def find_message_type(self):
        data_source = int(self.dataSrc)

        bit0 = data_source & 0b1  # flag for INAV
        bit1 = (data_source >> 1) & 0b1  # flag for FNAV
        bit2 = (data_source >> 2) & 0b1  # flag for INAV
        bit8 = (data_source >> 8) & 0b1  # flag for E1-E5a clock
        bit9 = (data_source >> 9) & 0b1  # flag for E1-E5b clock

        hasINAV = bool((bit0 | bit2) & bit9)
        hasFNAV = bool(bit1 & bit8)

        if hasINAV and hasFNAV:
            raise EphemerideError(f'Invalid galileo navigation message type. It cannot be INAV and FNAV '
                                  f'simultaneously. Please check field Data Source ({bin(data_source)})')
        if hasINAV:
            message_type = 'INAV'
        elif hasFNAV:
            message_type = 'FNAV'
        else:
            raise EphemerideError(f'Invalid galileo navigation message type ({bin(data_source)}). '
                                  f'No INAV or FNAV message available. Please check field Data Source')
        self.nav_type = message_type

        return self.nav_type

    @property
    def constellation(self):
        return "GAL"


class NavigationData:
    """
    NavigationDataMap
    this class stores data from rinex_parser navigation files
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
            try:
                self._data[satellite].set_data(epoch, nav_message)
            except TimeSeriesError as e:
                pass  # TODO: add warning to logger
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

    def get_closest_message(self, sat, epoch):
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
        except (TimeSeriesError, KeyError) as e:
            raise TimeSeriesError(f"satellite {str(sat)} has no available navigation data for epoch {str(epoch)}, {e}")

    @property
    def header(self):
        return self._header
