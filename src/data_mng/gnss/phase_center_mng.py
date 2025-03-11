""" Satellite Phase Center Manager Module
This module implements the Satellite Phase Center manager, that manages the computation of phase center
offsets and variations for the receiver and satellite antennas in the reconstruction of GNSS observations. """

from src.data_types.gnss import Satellite
from src.data_types.gnss.antena import ReceiverAntenna, SatelliteAntenna
from src.io.config import config_dict, EnumPCVModel, EnumAlgorithmPNT


class PhaseCenterManager:
    """
    Phase Center Manager Class.
    This class stores and manages the phase center offsets and variations for the receiver and satellite antennas.
    provided from ANTEX files.

    Furthermore, it provides the receiver reference point corrections as well, either provided from RINEX OBS
    (ANTENNA: DELTA H/E/N) or directly by the user in the configuration file.

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

        if config_dict.get("gnss_alg") == EnumAlgorithmPNT.SPS:
            self._enabled = False  # return with self._enabled set to False
            return

        # fetch from the user configurations
        self._enabled = config_dict.get("model", "phase_center_corrections", "enabled")
        rec_arp_enabled = config_dict.get("model", "phase_center_corrections", "receiver", "ARP_enabled")
        rec_pco_enabled = config_dict.get("model", "phase_center_corrections", "receiver", "PCO_enabled")
        rec_pcv_enabled = config_dict.get("model", "phase_center_corrections", "receiver", "PCV_enabled")
        rec_arp_offset = config_dict.get("model", "phase_center_corrections", "receiver", "ARP_offset")
        rec_pcv_model = EnumPCVModel(config_dict.get("model", "phase_center_corrections", "receiver", "PCV_model"))

        if rec_arp_offset is not None:
            if len(rec_arp_offset) != 3:
                raise AttributeError("Receiver ARP offset must be a 3-element list.")
            self._receiver_antenna.arp_offset = rec_arp_offset
        self._receiver_antenna.arp_enabled = rec_arp_enabled
        self._receiver_antenna.pco_enabled = rec_pco_enabled
        self._receiver_antenna.pcv_enabled = rec_pcv_enabled
        self._receiver_antenna.pcv_model = rec_pcv_model

    def add_satellite_antenna(self, satellite: Satellite, antenna: SatelliteAntenna):
        """ Add a satellite antenna to the manager. """
        self._satellite_antennas[satellite] = antenna
        # TODO: add here extra configurations for SatelliteAntenna PCO and PCV enabled.

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

    def __str__(self):
        """ String representation of the Phase Center Manager. """
        myStr = f"Phase Center Offsets and Variations:\n"
        myStr += f"\n-------------Receiver Antenna-------------\n{str(self._receiver_antenna)}"
        for sat, ant in self._satellite_antennas.items():
            myStr += f"-------------Antenna for {sat}-------------\n{str(ant)}"
        return myStr

    @property
    def enabled(self):
        """
        Returns:
            bool: True if the Phase Center Corrections are enabled, False otherwise.
        """
        return self._enabled
