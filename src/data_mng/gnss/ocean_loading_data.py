""" Ocean Loading Data Manager Module
"""

from src.io.config import config_dict
from src.io.rinex_parser.blq_reader import BLQReader


class OceanLoadingData:
    """
    This module stores the ocean loading data for the computation of Earth Deformation Effects of receiver positions.
    This data can be downloaded in  http://holt.oso.chalmers.se/loading/

    Ocean tide loading is the deformation of the Earth due to the weight of the ocean tides. The water in the ocean
    tides moves back and forth and these mass redistributions cause periodic loading of the ocean bottom. Since the
    Earth is not completely rigid, it deforms under this load and this is called ocean tide loading. One can observe it
    as variations at your station in vertical and horizontal displacement, in gravity, tilt and in strain. The ocean
    tides are produced by the gravitational pull of the Moon and Sun and since their orbits have more than one
    periodicities due to the eccentricity, evection and the lot, the ocean tides can be described as a sum of several
    ocean tides with each having their own period. The 11 periods, also called harmonics, with the largest amplitude
    are mostly used to compute the ocean tide loading.
    """
    def __init__(self):
        self._enabled = config_dict.get("model", "earth_deformation_effects", "enable_ocean_loading", fallback=False)
        self._station = config_dict.get("inputs", "ocean_loading", "station", fallback="")
        self._radial_amplitude = None
        self._west_amplitude = None
        self._south_amplitude = None
        self._radial_phase = None
        self._west_phase = None
        self._south_phase = None

    def init(self, file):
        """ Initialize the Ocean Loading Manager.

        Args:
            file (str): file path of the ocean loading file
        """
        if self.enabled:
            BLQReader(file, self, self._station)

    @property
    def radial_amplitude(self):
        """ Fetch Radial Amplitude parameters (meters). """
        return self._radial_amplitude

    @radial_amplitude.setter
    def radial_amplitude(self, data):
        if len(data) != 11:
            raise AttributeError(f"Provided data array (radial amplitude) should have 11 elements, but has {len(data)}")
        self._radial_amplitude = data

    @property
    def west_amplitude(self):
        """ Fetch West Amplitude parameters (meters). """
        return self._west_amplitude

    @west_amplitude.setter
    def west_amplitude(self, data):
        if len(data) != 11:
            raise AttributeError(f"Provided data array (west amplitude) should have 11 elements, but has {len(data)}")
        self._west_amplitude = data

    @property
    def south_amplitude(self):
        """ Fetch South Amplitude parameters (meters). """
        return self._south_amplitude

    @south_amplitude.setter
    def south_amplitude(self, data):
        if len(data) != 11:
            raise AttributeError(f"Provided data array (south amplitude) should have 11 elements, but has {len(data)}")
        self._south_amplitude = data

    @property
    def radial_phase(self):
        """ Fetch Radial Phase parameters (degrees). """
        return self._radial_phase

    @radial_phase.setter
    def radial_phase(self, data):
        if len(data) != 11:
            raise AttributeError(f"Provided data array (radial phase) should have 11 elements, but has {len(data)}")
        self._radial_phase = data

    @property
    def west_phase(self):
        """ Fetch West Phase parameters (degrees). """
        return self._west_phase

    @west_phase.setter
    def west_phase(self, data):
        if len(data) != 11:
            raise AttributeError(f"Provided data array (west phase) should have 11 elements, but has {len(data)}")
        self._west_phase = data

    @property
    def south_phase(self):
        """ Fetch South Phase parameters (degrees). """
        return self._south_phase

    @south_phase.setter
    def south_phase(self, data):
        if len(data) != 11:
            raise AttributeError(f"Provided data array (south phase) should have 11 elements, but has {len(data)}")
        self._south_phase = data

    @property
    def enabled(self):
        """
        Returns:
            bool: True if the Ocean Loading Parameter Correction is enabled, False otherwise.
        """
        return self._enabled

    def __str__(self):
        """ String representation of the Ocean Loading Parameters. """
        myStr = f"Ocean Loading Parameters:\n"
        myStr += f"Station: {self._station}\n"
        myStr += f"Radial Amplitude Array: {self.radial_amplitude}\n"
        myStr += f"West Amplitude Array: {self.west_amplitude}\n"
        myStr += f"South Amplitude Array: {self.south_amplitude}\n"
        myStr += f"Radial Phase Array: {self.radial_phase}\n"
        myStr += f"West Phase Array: {self.west_phase}\n"
        myStr += f"South Phase Array: {self.south_phase}"
        return myStr
