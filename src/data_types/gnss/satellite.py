""" Satellite Module """
from .constellation import Constellation, get_constellation

__all__ = ["get_satellite", "Satellite"]

# list of all satellites for the SatelliteFactory
__cache__ = []


class Satellite:
    """
    Class Satellite
    Representation of a satellite with the important variables associated with it

    Attributes:
        nID (int): The integer identifier of the satellite
        sat_system (Constellation): The instance of the corresponding SatelliteSystem (constellation)
    """

    def __init__(self, identifier, sat_system=None):
        """
        Constructor of Satellite instances

        Args:
            identifier (str or int): the satellite identifier, either an int nID (ex: 03), or a string composed of
            constellation + nID (ex: "G03")
            sat_system (Constellation): instance of Constellation. Must be provided if `identifier` is of type int,
                and may be skipped if it is of type str (since constellation is inferred from the string)
        Raises:
            AttributeError: an exception is raised if one of the input arguments if not of the specified type
        """

        if isinstance(identifier, int) and isinstance(sat_system, Constellation):
            self.nID = identifier
            self.sat_system = sat_system

        elif isinstance(identifier, str) and (len(identifier) == 3):
            self.nID = int(identifier[1:])
            self.sat_system = get_constellation(identifier[0])

        else:
            raise AttributeError("Unable to initialize satellite with arguments {} and {}. See documentation."
                                 "".format(identifier, sat_system))

        new_sat(self)

    def __repr__(self):
        """ A repr is useful for debugging """
        # sat_str = str(self.nID) if self.nID > 9 else '0' + str(self.nID)
        # return f'{type(self).__name__}({sat_str}, {self.sat_system})'
        return str(self)

    def __str__(self):
        sat_str = str(self.nID) if self.nID > 9 else '0' + str(self.nID)
        return "{}{}".format(self.sat_system.get_system_short(), sat_str)

    def __eq__(self, other):
        return self.nID == other.nID and self.sat_system == other.sat_system

    def __hash__(self):
        # used for hashing Satellites in lists or dict keys
        return hash(str(self))

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not (self == other)


def new_sat(sat: Satellite):
    """ Saves the newly created Satellite instance to the internal cache """
    if sat not in __cache__:
        __cache__.append(sat)


def get_satellite(sat_string: str) -> Satellite:
    """ Satellite Factory. Creates and returns the Satellite instance from the provided string name """
    for sat in __cache__:
        if str(sat) == sat_string:
            return sat
    else:
        return Satellite(sat_string)
