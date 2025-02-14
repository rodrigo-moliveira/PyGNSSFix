import numpy

from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.satellite import Satellite

# TODO: add documentation


class PhaseCenter:
    def __init__(self):
        self.pcv_noazi = None  # PCV values for non-azimuth dependent PCV
        self.pcv_azi = None  # PCV values for azimuth dependent PCV
        self.pco = None  # PCO values for the phase center offset


class Antenna:
    def __init__(self):
        self._ant_type = None  # Antenna type
        self.zenith_vec = None  # Zenith angles for the antenna pattern (in degrees)
        self.azimuth_vec = None  # Azimuth angles for the antenna pattern (in degrees)
        self.freq_data = {}  # PhaseCenter data for the antenna for each frequency

    @property
    def ant_type(self):
        return self._ant_type

    @ant_type.setter
    def ant_type(self, value: str):
        if not isinstance(value, str):
            raise AttributeError("Antenna type must be a string.")
        self._ant_type = value

    def set_freq_data(self, freq_type: DataType, freq_data: PhaseCenter):
        if not isinstance(freq_type, DataType):
            raise AttributeError("freq_type must be a DataType")
        if not isinstance(freq_data, PhaseCenter):
            raise AttributeError("freq_data must be a PhaseCenter")
        self.freq_data[freq_type] = freq_data


class SatelliteAntenna(Antenna):
    def __init__(self, satellite: Satellite):
        super().__init__()
        self.satellite = satellite

    @property
    def satellite(self):
        return self._satellite

    @satellite.setter
    def satellite(self, value: Satellite):
        if not isinstance(value, Satellite):
            raise AttributeError("satellite must be a Satellite.")
        self._satellite = value


class ReceiverAntenna(Antenna):
    def __init__(self):
        super().__init__()
        self._serial_no = None  # Antenna serial number
        self._arp_offset = None  # Antenna Reference Point offset from ANTENNA: DELTA H/E/N (meters)

    @property
    def serial_no(self):
        return self._serial_no

    @serial_no.setter
    def serial_no(self, value: int):
        if not isinstance(value, int):
            raise AttributeError("Antenna Serial Number must be an integer.")
        self._serial_no = value

    @property
    def arp_offset(self):
        return self._arp_offset

    @arp_offset.setter
    def arp_offset(self, value: numpy.ndarray):
        # convert list or tuple to numpy array
        if isinstance(value, list) or isinstance(value, tuple):
            value = numpy.array(value)

        if not isinstance(value, numpy.ndarray):
            raise AttributeError("arp_offset must be a numpy.ndarray")
        if value.size != 3:
            raise AttributeError("arp_offset must have a size of 3")
        self._arp_offset = value
