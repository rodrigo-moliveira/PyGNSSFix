""" Navigation Data Module.
This module implements classes to store the Navigation Messages from RINEX Navigation files
"""
from collections import OrderedDict

from src.errors import TimeSeriesError, EphemerideError
from src.data_types.date import Epoch
from src.data_types.gnss import Satellite
from src.data_mng import TimeSeries, Container

__all__ = ["NavigationData"]


class NavigationHeader(Container):
    """
    NavigationHeader class, inherits from the :py:class:`Container` class.
    Stores relevant data from the header section of a RINEX Navigation file.

    Attributes:
        iono_corrections(dict)
        leap_seconds(int)
        time_correction(dict)
    """
    __slots__ = ["iono_corrections", "leap_seconds", "time_correction"]

    def __init__(self):
        super().__init__()
        for attr in self.__slots__:
            setattr(self, attr, None)
        self.iono_corrections = {}
        self.time_correction = {}
        self.leap_seconds = 0


class NavigationPoint(Container):
    """
    NavigationPoint class, inherits from the :py:class:`Container` class.

    This is the base class to store GPS and GAL specific navigation messages.
    The navigation messages refer to the data contained in the RINEX Navigation files.
    """

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
    NavigationPointGPS class, inherits from the :py:class:`NavigationPoint` class.
    Stores the data contained in a single navigation message for GPS satellites.
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
    NavigationPointGPS class, inherits from the :py:class:`NavigationPoint` class.
    stores the data contained in a single navigation message for GAL satellites.
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
        """
        Checks if this ephemeride point is I/NAV or F/NAV and fills the attribute `nav_type` accordingly.

        See RINEX Navigation documentation (https://igs.org/wg/rinex/#documents-formats) for more information about the
        `Data Source` field.

        Raises:
            EphemerideError: an exception is raised if the `Data Source` (attribute `dataSrc`) field in the navigation
                message is not a valid one.
        """
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
    Navigation Data Frame.
    This class stores the navigation data for all satellites, that is read from RINEX Navigation files
    """

    def __init__(self):
        self._data = OrderedDict()
        self._header = NavigationHeader()

    def __str__(self):
        """ Print the navigation data to a string, for debug purposes """
        myStr = f"{repr(self._header)}\n"
        for sat, data in self._data.items():
            myStr += str(sat) + f"\n{str(data)}\n"

        return myStr

    def set_data(self, epoch: Epoch, satellite: Satellite, nav_message: NavigationPoint):
        """
        Method to set a navigation data point for a given epoch and satellite.

        Args:
            epoch (Epoch): the epoch correspondent to the Time of Clock (TOC) of the message
            satellite (Satellite): the satellite associated with the message
            nav_message (NavigationPoint): the message object (instance of `NavigationPoint`)

        Raises:
            AttributeError: an exception is raised if one of the input arguments is not of the specified type.
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
        Method to set the navigation header.

        Args:
            header (NavigationHeader): instance of `NavigationHeader` to be set.

        Raises:
            AttributeError: an exception is raised if the input header is not an instance of the required
                type (`NavigationHeader`)
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
        Gets navigation message of the given satellite closest to the provided epoch.

        Args:
            sat (Satellite)
            epoch (Epoch)

        Returns:
            NavigationPoint: navigation data for the given satellite closest to the provided epoch

        Raises:
            TimeSeriesError: an exception is raised if no valid navigation message was fetched for the given satellite
        """
        try:
            _epoch = self._data[sat].get_closest_epoch(epoch)
            return self._data[sat].get_data_for_epoch(_epoch)
        except (TimeSeriesError, KeyError) as e:
            raise TimeSeriesError(f"satellite {str(sat)} has no available navigation data for epoch {str(epoch)}, {e}")

    @property
    def header(self):
        return self._header

    def is_empty(self) -> bool:
        """ Returns True if the navigation data is empty. """
        return len(self._data) == 0
