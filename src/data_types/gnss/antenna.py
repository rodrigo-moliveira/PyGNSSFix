""" Antenna and Phase Center Module. """
import numpy as np

from src.data_types.gnss.data_type import DataType
from src.data_types.gnss.satellite import Satellite
from src.constants import RAD2DEG
from src.io.config.enums import EnumPCVModel
from src.utils.interpolation import linear_interpolation_scipy, binary_search


class PhaseCenter:
    """ Phase Center class. This class holds the PCO and PCV values for a single antenna.
    PCV data can be azimuth dependent or non-azimuth dependent.

    Attributes:
        pcv_noazi (list): List of PCV values in millimeters for non-azimuth dependent PCV
        pcv_azi (numpy.ndarray): Numpy matrix with PCV values in millimeters for azimuth dependent PCV.
            Row axis is azimuth and column axis is zenith.
        pco (list): List with PCO values for the phase center offset, in the format of the ANTEX file, that is:
            * NORTH / EAST / UP in millimeters for receiver antennas
            * X / Y / Z in millimeters for satellite antennas (body-fixed frame coordinates)
    """
    def __init__(self):
        self.pcv_noazi = None  # PCV values for non-azimuth dependent PCV
        self.pcv_azi = None  # PCV values for azimuth dependent PCV
        self.pco = None  # PCO values for the phase center offset

    def __str__(self):
        """ String representation of the PhaseCenter. """
        return f"PCO: {self.pco} [mm]\nPCV No Azimuth: {self.pcv_noazi} [mm]\nPCV Azimuth:\n{self.pcv_azi} [mm]\n"


class Antenna:
    """ Antenna base class.

    Attributes:
        ant_type (str): Antenna type
        zenith_vec (numpy.ndarray): Zenith angles for the antenna pattern (in degrees)
        azimuth_vec (numpy.ndarray): Azimuth angles for the antenna pattern (in degrees)
        freq_data (dict): dict of PhaseCenter data for the antenna for each frequency
        pco_enabled (bool): PCO enabled flag
        pcv_enabled (bool): PCV enabled flag
        pcv_model (EnumPCVModel): PCV model
    """
    def __init__(self):
        """
        Constructor for the Antenna.
        """
        self.ant_type = ""
        self.zenith_vec = None
        self.azimuth_vec = None
        self.freq_data = {}
        self.pco_enabled = False
        self.pcv_enabled = False
        self._pcv_model = None
        self._trace_file = None
        self._warning_thrown = False
        self.write_trace = False

    def __del__(self):
        """ Destructor for the Antenna: close the trace file. """
        if self._trace_file is not None:
            self._trace_file.close()

    def __str__(self):
        """ String representation of the Antenna. """
        with np.printoptions(threshold=np.inf, linewidth=np.inf):
            myStr = f"Antenna Type: {self.ant_type}\nZenith Vector: {self.zenith_vec}\n" \
                    f"Azimuth Vector: {self.azimuth_vec}\nPCO Enabled: {self.pco_enabled}\n" \
                    f"PCV Enabled: {self.pcv_enabled}\nPCV Model: {self._pcv_model}\n\n"
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

    @property
    def pcv_model(self):
        """
        Returns:
             EnumPCVModel: the PCV model for the antenna.
         """
        return self._pcv_model

    @pcv_model.setter
    def pcv_model(self, value: EnumPCVModel):
        """ Set the PCV model for the antenna. """
        if not isinstance(value, EnumPCVModel):
            raise AttributeError("PCV model must be an EnumPCVModel.")
        self._pcv_model = value

    @property
    def trace_file(self):
        """ Returns the trace file for the antenna. """
        return self._trace_file

    @trace_file.setter
    def trace_file(self, value):
        """ Set the trace file for the antenna. """
        if self._pcv_model == EnumPCVModel.AZI_DEPENDENT:
            # only create trace file if the PCV model is azimuth dependent (interpolations required)
            if value is not None:
                self.write_trace = True
                self._trace_file = open(value, "w")
            else:
                self.write_trace = False
                self._trace_file = None

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
        """ Get the PCV value for the antenna.
        This function performs interpolation to get the PCV value for the given azimuth and elevation angles.

        For receiver antennas, the elevation and azimuth angles refer to the local (ENU) frame.
            - Azimuth counts clockwise from the North toward the East.

        For satellite antennas, the elevation and azimuth angles refer to the body-fixed frame.
            - Azimuth counts clockwise from the y-axis toward the x-axis when looking in the direction of the
                negative z-axis or toward deep space.

        Args:
            freq_type (DataType): Frequency type to get the PCV for
            azimuth (float): Azimuth angle in radians
            elevation (float): Elevation angle in radians

        Returns:
            float: PCV value in meters

        Raises:
            KeyError: If the frequency type is not found in the antenna data
        """
        # convert elevation to zenith angle
        zenith_deg = 90 - elevation*RAD2DEG  # in degrees
        azimuth_deg = azimuth*RAD2DEG  # in degrees

        if freq_type is None or freq_type not in self.freq_data:
            raise KeyError(f"Frequency {freq_type} not found in the antenna data.")

        pcv_out = 0.0
        if self.pcv_model == EnumPCVModel.AZI_DEPENDENT and self.freq_data[freq_type].pcv_azi is not None:
            # only apply azimuth dependent PCV if data is available
            pcv_matrix = self.freq_data[freq_type].pcv_azi
            pcv_out = self._interpolate_pcv_matrix(zenith_deg, azimuth_deg, pcv_matrix)

        elif self.freq_data[freq_type].pcv_noazi is not None:
            # else, apply non-azimuth dependent PCV
            if self.pcv_model == EnumPCVModel.AZI_DEPENDENT:
                if not self._warning_thrown:
                    self._warning_thrown = True
                    from src.common_log import get_logger, MODEL_LOG
                    log = get_logger(MODEL_LOG)
                    log.warning(f"Resorting to non-azimuth dependent PCV model for the antenna "
                                f"{self.ant_type} and frequency {freq_type}.")
            pcv_values = self.freq_data[freq_type].pcv_noazi
            pcv_out = linear_interpolation_scipy(zenith_deg, self.zenith_vec, pcv_values)
        else:
            if not self._warning_thrown:
                self._warning_thrown = True
                from src.common_log import get_logger, MODEL_LOG
                log = get_logger(MODEL_LOG)
                log.warning(f"No PCV data available for the antenna {self.ant_type} and frequency {freq_type}. "
                            f"Setting PCV to 0.")
        return float(pcv_out)
    
    def _interpolate_pcv_matrix(self, zenith: float, azimuth: float, pcv_matrix: np.ndarray) -> float:
        """ Interpolate the PCV matrix to get the PCV value for the given zenith and azimuth angles.

        Args:
            zenith (float): Zenith angle in degrees
            azimuth (float): Azimuth angle in degrees
            pcv_matrix (numpy.ndarray): PCV matrix with azimuth and zenith angles

        Returns:
            float: Interpolated PCV value in meters

        Raises:
            ValueError: If the interpolation coefficients are not within the expected range
        """
        # find the closest latitude and longitude grid points and indexes in array
        grid_zen, idx_zen = binary_search(list(self.zenith_vec), zenith, 1, ret_index=True, extrapolation=True)
        grid_az, idx_az = binary_search(list(self.azimuth_vec), azimuth, 1, ret_index=True, extrapolation=True)

        if self.write_trace:
            self.trace_file.write(f"[GRID POINT] PCV for zenith {grid_zen[0]} and azimuth {grid_az[0]}: "
                                  f"{pcv_matrix[idx_az[0], idx_zen[0]]}\n")
            self.trace_file.write(f"[GRID POINT] PCV for zenith {grid_zen[1]} and azimuth {grid_az[1]}: "
                                  f"{pcv_matrix[idx_az[1], idx_zen[1]]}\n")
            self.trace_file.write(f"[GRID POINT] PCV for zenith {grid_zen[0]} and azimuth {grid_az[1]}: "
                                  f"{pcv_matrix[idx_az[1], idx_zen[0]]}\n")
            self.trace_file.write(f"[GRID POINT] PCV for zenith {grid_zen[1]} and azimuth {grid_az[0]}: "
                                  f"{pcv_matrix[idx_az[0], idx_zen[1]]}\n")

        zen_coeff = (zenith - grid_zen[0]) / (grid_zen[1] - grid_zen[0])
        az_coeff = (azimuth - grid_az[0]) / (grid_az[1] - grid_az[0])

        if 0 <= zen_coeff <= 1 and 0 <= az_coeff <= 1:
            pcv = (1 - zen_coeff) * (1 - az_coeff) * pcv_matrix[idx_az[0], idx_zen[0]] + \
                   zen_coeff * (1 - az_coeff) * pcv_matrix[idx_az[0], idx_zen[1]] + \
                   (1 - zen_coeff) * az_coeff * pcv_matrix[idx_az[1], idx_zen[0]] + \
                   zen_coeff * az_coeff * pcv_matrix[idx_az[1], idx_zen[1]]
            if self.write_trace:
                self.trace_file.write(f"\t[INTERPOLATED] PCV for zenith {zenith} and azimuth {azimuth}: {pcv}\n\n")
        else:
            # nearest neighbor extrapolation
            if zen_coeff <= 0.5 and az_coeff <= 0.5:
                pcv = pcv_matrix[idx_az[0], idx_zen[0]]
            elif zen_coeff <= 0.5 and az_coeff > 0.5:
                pcv = pcv_matrix[idx_az[1], idx_zen[0]]
            elif zen_coeff > 0.5 and az_coeff <= 0.5:
                pcv = pcv_matrix[idx_az[0], idx_zen[1]]
            elif zen_coeff > 0.5 and az_coeff > 0.5:
                pcv = pcv_matrix[idx_az[1], idx_zen[1]]
            else:
                raise ValueError(f"Error in interpolation coefficients: zen_coeff={zen_coeff}, az_coeff={az_coeff}")
            if self.write_trace:
                self.trace_file.write(f"\t[EXTRAPOLATED] PCV zenith {zenith} and azimuth {azimuth}: {pcv} "
                                      f"(pulled to the nearest neighbor zen_coeff={zen_coeff}, "
                                      f"az_coeff={az_coeff})\n")
        return pcv


class SatelliteAntenna(Antenna):
    """ Satellite Antenna class, inherits from Antenna base class.

    Attributes:
        satellite (Satellite): Satellite object
    """
    def __init__(self, satellite: Satellite):
        """ Constructor for the Satellite Antenna.

        Args:
            satellite (Satellite): Satellite object
        """
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
        arp_offset (numpy.ndarray): Antenna Reference Point offset from ANTENNA: DELTA H/E/N (meters)
        arp_enabled (bool): ARP enabled flag
    """
    def __init__(self):
        """ Constructor for the Receiver Antenna. """
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
    def arp_offset(self, value: np.ndarray):
        """ Set the antenna reference point offset. """
        # convert list or tuple to numpy array
        if isinstance(value, list) or isinstance(value, tuple):
            value = np.array(value)

        if not isinstance(value, np.ndarray):
            raise AttributeError("arp_offset must be a numpy.ndarray")
        if value.size != 3:
            raise AttributeError("arp_offset must have a size of 3")
        self._arp_offset = value

    def __str__(self):
        """ String representation of the Receiver Antenna. """
        return f"Antenna Type: {self.ant_type}\nSerial No.: {self.serial_no}\nARP Enabled: {self.arp_enabled}\n" \
               f"ARP Offset: {self.arp_offset}\n" \
               f"{super().__str__()}"
