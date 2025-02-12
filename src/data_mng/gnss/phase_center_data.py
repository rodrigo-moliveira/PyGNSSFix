import numpy

from src.data_mng.gnss.antena_pcv import AntennaPCV


class PhaseCenterData:
    # TODO: write documentation
    def __init__(self):
        self._receiver_antenna_type = None  # ANT # / TYPE
        self._receiver_antenna_serial_no = None
        self._receiver_arp_offset = None # Antenna Reference Point offset from ANTENNA: DELTA H/E/N
        self.receiver_pco = {}  # Receiver Antenna Phase Center Offset. Dict indexed by frequency
        self.receiver_pcv = AntennaPCV()  # TODO: here create a new class for PCV (with azimuth and non-azimuth tables).

    @property
    def receiver_antenna_type(self):
        return self._receiver_antenna_type

    @receiver_antenna_type.setter
    def receiver_antenna_type(self, value: str):
        if not isinstance(value, str):
            raise AttributeError("receiver_antenna_type must be a string")
        self._receiver_antenna_type = value

    @property
    def receiver_arp_offset(self):
        return self._receiver_arp_offset

    @receiver_arp_offset.setter
    def receiver_arp_offset(self, value: numpy.ndarray):
        # convert list or tuple to numpy array
        if isinstance(value, list) or isinstance(value, tuple):
            value = numpy.array(value)

        if not isinstance(value, numpy.ndarray):
            raise AttributeError("receiver_arp_offset must be a numpy.ndarray")
        if value.size != 3:
            raise AttributeError("receiver_arp_offset must have a size of 3")
        self._receiver_arp_offset = value

    @property
    def receiver_antenna_serial_no(self):
        return self._receiver_antenna_serial_no

    @receiver_antenna_serial_no.setter
    def receiver_antenna_serial_no(self, value: int):
        if not isinstance(value, int):
            raise AttributeError("receiver_antenna_serial_no must be an integer")
        self._receiver_antenna_serial_no = value

