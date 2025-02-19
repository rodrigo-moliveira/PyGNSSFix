""" Satellite Phase Center Manager Module
This module implements the Satellite Phase Center manager, that manages the computation of phase center
offsets and variations for the receiver and satellite antennas in the reconstruction of GNSS observations. """

from src.data_types.gnss import Satellite
from src.data_types.gnss.antena import ReceiverAntenna, SatelliteAntenna


class PhaseCenterManager:
    """
    Phase Center Manager Class.
    This class stores and manages the phase center offsets and variations for the receiver and satellite antennas.
    provided from ANTEX files.

    Attributes:
        _receiver_antenna(ReceiverAntenna): attribute that stores the receiver antenna phase center information.
        _satellite_antennas(dict): attribute that stores the satellite antenna phase center information.
            It is a dictionary defined as:
                * keys -> :py:class:`Satellite` instances.
                * values -> :py:class:`SatelliteAntenna` instances with the antenna phase center
    """
    def __init__(self):
        self._receiver_antenna = ReceiverAntenna()
        self._satellite_antennas = {}

    def add_satellite_antenna(self, satellite: Satellite, antenna: SatelliteAntenna):
        """ Add a satellite antenna to the manager. """
        self._satellite_antennas[satellite] = antenna

    def get_satellite_antenna(self, satellite: Satellite):
        """ Get the antenna for a given satellite. """
        if satellite not in self._satellite_antennas:
            raise AttributeError(f"Antenna for satellite {satellite} not loaded from ANTEX file.")
        return self._satellite_antennas.get(satellite)

    def set_receiver_antenna(self, antenna: ReceiverAntenna):
        """ Set the receiver antenna for the manager. """
        self._receiver_antenna = antenna

    def get_receiver_antenna(self):
        """ Get the receiver antenna for the manager. """
        return self._receiver_antenna
