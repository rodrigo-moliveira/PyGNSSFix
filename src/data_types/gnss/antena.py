""" Antenna Module """
import numpy

from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.satellite import Satellite
from src.constants import PI
from src.utils.interpolation import linear_interpolation_scipy


class PhaseCenter:
    def __init__(self):
        self.pcv_noazi = None  # PCV values for non-azimuth dependent PCV
        self.pcv_azi = None  # PCV values for azimuth dependent PCV
        self.pco = None  # PCO values for the phase center offset

    def __str__(self):
        """ String representation of the PhaseCenter. """
        return f"PCO: {self.pco}\nPCV No Azimuth: {self.pcv_noazi}\nPCV Azimuth:\n{self.pcv_azi}\n"


class Antenna:
    """ Antenna base class.

    Attributes:
        ant_type (str): Antenna type
        zenith_vec (numpy.ndarray): Zenith angles for the antenna pattern (in degrees)
        azimuth_vec (numpy.ndarray): Azimuth angles for the antenna pattern (in degrees)
        freq_data (dict): dict of PhaseCenter data for the antenna for each frequency
    """
    def __init__(self):
        self.ant_type = ""
        self.zenith_vec = None
        self.azimuth_vec = None
        self.freq_data = {}
        self.pco_enabled = False
        self.pcv_enabled = False
        self.pcv_model = None

    def __str__(self):
        """ String representation of the Antenna. """
        with numpy.printoptions(threshold=numpy.inf, linewidth=numpy.inf):
            myStr = f"Antenna Type: {self.ant_type}\nZenith Vector: {self.zenith_vec}\n" \
                    f"Azimuth Vector: {self.azimuth_vec}\nPCO Enabled: {self.pco_enabled}\n" \
                    f"PCV Enabled: {self.pcv_enabled}\nPCV Model: {self.pcv_model}\n\n"
            for freq, data in self.freq_data.items():
                myStr += f"Frequency: {freq}\n{str(data)}\n"
            return myStr

    @property
    def ant_type(self):
        """
        Returns:
            str: Antenna type
        """
        return self._ant_type

    @ant_type.setter
    def ant_type(self, value: str):
        """ Set the antenna type. """
        if not isinstance(value, str):
            raise AttributeError("Antenna type must be a string.")
        self._ant_type = value

    def set_freq_data(self, freq_type: DataType, freq_data: PhaseCenter):
        """ Set the frequency data for the antenna.

        Args:
            freq_type (DataType): Frequency type
            freq_data (PhaseCenter): Phase center data for the frequency
         """
        if not isinstance(freq_type, DataType):
            raise AttributeError("freq_type must be a DataType")
        if not isinstance(freq_data, PhaseCenter):
            raise AttributeError("freq_data must be a PhaseCenter")
        self.freq_data[freq_type] = freq_data

    def get_pcv(self, freq_type: DataType, azimuth, elevation):
        # perform interpolation
        # Only zenith, first
        zenith_angle = 90 - elevation*180/PI
        # zen_angles = self.freq_data[freq_type]
        pco_values = self.freq_data[freq_type].pcv_noazi
        #knot_zenith, idx_zen = binary_search(list(self.zenith_vec), zenith_angle, 1, ret_index=True, extrapolation=True)

        # TODO: simply check this
        #pco_values2 = [pco_values[i] for i in idx_zen]
        #intep1 = linear_interpolation_scipy(zenith_angle, knot_zenith, pco_values2)
        intep2 = linear_interpolation_scipy(zenith_angle, self.zenith_vec, pco_values)
        return float(intep2)

class SatelliteAntenna(Antenna):
    """ Satellite Antenna class, inherits from Antenna base class.

    Attributes:
        satellite (Satellite): Satellite object
    """
    def __init__(self, satellite: Satellite):
        super().__init__()
        self.satellite = satellite

    def __str__(self):
        """ String representation of the Satellite Antenna. """
        return f"Satellite: {self.satellite}\n{super().__str__()}"

    @property
    def satellite(self):
        """
        Returns:
            Satellite: The instance of the satellite
        """
        return self._satellite

    @satellite.setter
    def satellite(self, value: Satellite):
        """ Set the satellite object. """
        if not isinstance(value, Satellite):
            raise AttributeError("satellite must be a Satellite.")
        self._satellite = value


class ReceiverAntenna(Antenna):
    """ Receiver Antenna class, inherits from Antenna base class.

    Attributes:
        serial_no (int): Antenna serial number
        arp_offset (numpy.ndarray): Antenna Reference Point offset from ANTENNA: DELTA H/E/N (meters
    """
    def __init__(self):
        super().__init__()
        self.serial_no = 0  # Antenna serial number
        self.arp_offset = [0, 0, 0]
        self.arp_enabled = False

    @property
    def serial_no(self):
        """
        Returns:
             int: Antenna Serial Number
        """
        return self._serial_no

    @serial_no.setter
    def serial_no(self, value: int):
        """ Set the antenna serial number. """
        if not isinstance(value, int):
            raise AttributeError("Antenna Serial Number must be an integer.")
        self._serial_no = value

    @property
    def arp_offset(self):
        """
        Returns:
             numpy.ndarray: Antenna Reference Point offset from ANTENNA: DELTA H/E/N (meters)
        """
        return self._arp_offset

    @arp_offset.setter
    def arp_offset(self, value: numpy.ndarray):
        """ Set the antenna reference point offset. """
        # convert list or tuple to numpy array
        if isinstance(value, list) or isinstance(value, tuple):
            value = numpy.array(value)

        if not isinstance(value, numpy.ndarray):
            raise AttributeError("arp_offset must be a numpy.ndarray")
        if value.size != 3:
            raise AttributeError("arp_offset must have a size of 3")
        self._arp_offset = value

    def __str__(self):
        """ String representation of the Receiver Antenna. """
        return f"Antenna Type: {self.ant_type}\nSerial No.: {self.serial_no}\nARP Enabled: {self.arp_enabled}\n" \
               f"ARP Offset: {self.arp_offset}\n" \
               f"{super().__str__()}"
