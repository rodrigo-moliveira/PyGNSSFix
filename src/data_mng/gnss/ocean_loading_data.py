""" Ocean Loading Data Manager Module
This module stores the ocean loading data for the computation of Earth Deformation Effects of receiver positions.
This data can be downloaded in  http://holt.oso.chalmers.se/loading/

Ocean tide loading is the deformation of the Earth due to the weight of the ocean tides. The water in the ocean tides
moves back and forth and these mass redistributions cause periodic loading of the ocean bottom. Since the Earth is not
completely rigid, it deforms under this load and this is called ocean tide loading. One can observe it as variations at
your station in vertical and horizontal displacement, in gravity, tilt and in strain. The ocean tides are produced by
the gravitational pull of the Moon and Sun and since their orbits have more than one periodicities due to the
eccentricity, evection and the lot, the ocean tides can be described as a sum of several ocean tides with each having
their own period. The 11 periods, also called harmonics, with the largest amplitude are mostly used to compute the ocean
tide loading.
"""

from src.io.config import config_dict
from src.io.rinex_parser.blq_reader import BLQReader


class OceanLoadingData:
    def __init__(self):
        self._enabled = config_dict.get("model", "earth_deformation_effects", "enable_ocean_loading")
        self._station = config_dict.get("inputs", "ocean_loading", "station")
        self._data = {}

    def init(self, file):
        """ Initialize the Ocean Loading Manager.

        Args:
            file (str): file path of the ocean loading file
        """
        if self.enabled:
            BLQReader(file, self, self._station)

    def set_data(self):
        pass

    @property
    def enabled(self):
        """
        Returns:
            bool: True if the Ocean Loading Parameter Correction is enabled, False otherwise.
        """
        return self._enabled
