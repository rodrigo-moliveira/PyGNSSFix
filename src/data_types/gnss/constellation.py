from .service_utils import AvailableConstellations, ConstellationToCodeMap
from ...errors import UnknownConstellationError

__all__ = ["Constellation", "get_constellation"]


class Constellation(str):
    """
    Class Constellation, inherits from string
    Represents a constellation system (GPS or GAL)
    This class is a subclass of str with special utility methods
    """

    def __new__(cls, content):

        if content.upper() in AvailableConstellations:
            return super(Constellation, cls).__new__(cls, content.upper())
        else:
            raise ValueError("Unknown Satellite System {}. Possible systems are GPS or GAL".format(content.upper()))

    def __repr__(self):
        """A repr is useful for debugging"""
        return f'{type(self).__name__}({super().__repr__()})'

    def __getattribute__(self, name):
        if name in dir(str):  # only handle str methods here

            def method(self, *args, **kwargs):
                value = getattr(super(), name)(*args, **kwargs)
                # not every string method returns a str:
                if isinstance(value, str):
                    return type(self)(value)
                elif isinstance(value, list):
                    return [type(self)(i) for i in value]
                elif isinstance(value, tuple):
                    return tuple(type(self)(i) for i in value)
                else:  # dict, bool, or int
                    return value

            return method.__get__(self)  # bound method
        else:  # delegate to parent
            return super().__getattribute__(name)

    def is_gps(self):
        return self == "GPS"

    def is_gal(self):
        return self == "GAL"

    def get_system(self):
        return type(self)(self)

    def get_system_short(self):
        """
        Return:
            str : "G" for sensors, "E" for galileo "U" for unkown
        """
        return ConstellationToCodeMap.get(self, "UNKNOWN")


_GPS = Constellation("GPS")
_GAL = Constellation("GAL")


def get_constellation(system: str):
    """
    Method to return the Constellation object

    Args:
        system (str) : short name for the constellation
    Return:
        SatelliteSystem : the corresponding system object
    """
    if system.upper() == "GPS" or system.upper() == "G":
        return _GPS
    elif system.upper() == "GAL" or system.upper() == "E":
        return _GAL
    else:
        raise UnknownConstellationError(
            "No Satellite System matched the descriptor {}.".format(system)
        )
