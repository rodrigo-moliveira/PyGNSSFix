from PositioningSolver.src.gnss.data_types.Constellation import SatelliteSystem, SatelliteSystemFactory

# list of all satellites for the SatelliteFactory
__all_sats__ = []


def new_sat(sat):
    if sat not in __all_sats__:
        __all_sats__.append(sat)


def SatelliteFactory(sat_string):
    for sat in __all_sats__:
        if str(sat) == sat_string:
            return sat
    else:
        return Satellite(sat_string)


class Satellite:
    """
    Class Satellite
    Representation of a satellite with the important variables associated with it

    Attributes
        ----------
        nID : int
            The integer identifier of the satellite
        sat_system : SatelliteSystem
            The instance of the corresponding SatelliteSystem (constellation)

    """

    def __init__(self, identifier, sat_system: SatelliteSystem = None):
        """
        Constructor of Satellite instances

        Args:
            identifier (str, int):
                the satellite identifier, either an int (ex: 03) - just
                the ID, or a string composed of constellation + ID (ex: "G03")
            sat_system (SatelliteSystem):
                instance of SatelliteSystem. Must be provided if ´´identifier´´ is of type int,
                and may be skipped if it is of type str (since constellation is inferred from the string)
        """

        if isinstance(identifier, int) and isinstance(sat_system, SatelliteSystem):
            self._nID = identifier
            self._sat_system = sat_system

        elif isinstance(identifier, str) and (len(identifier) == 3):
            self._nID = int(identifier[1:])
            self._sat_system = SatelliteSystemFactory(identifier[0])

        else:
            raise TypeError("Unable to initialize satellite with arguments {} and {}. See documentation."
                            "".format(identifier, sat_system))

        new_sat(self)

    @property
    def nID(self):
        return self._nID

    @property
    def sat_system(self):
        return self._sat_system

    def __repr__(self):
        """A repr is useful for debugging"""
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
