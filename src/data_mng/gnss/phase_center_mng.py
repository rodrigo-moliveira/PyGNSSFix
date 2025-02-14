from src.data_types.gnss import Satellite
from src.data_types.gnss.antena import Antenna


class PhaseCenterManager:
    # TODO: write documentation
    def __init__(self):
        self._receiver_antenna = Antenna()
        self._satellite_antennas = {}

    def add_satellite_antenna(self, satellite: Satellite, antenna: Antenna):
        self._satellite_antennas[satellite] = antenna

    def get_satellite_antenna(self, satellite: Satellite):
        if satellite not in self._satellite_antennas:
            raise AttributeError(f"Antenna for satellite {satellite} not loaded from ANTEX file.")
        return self._satellite_antennas.get(satellite)

    def set_receiver_antenna(self, antenna: Antenna):
        self._receiver_antenna = antenna

    def get_receiver_antenna(self):
        return self._receiver_antenna
