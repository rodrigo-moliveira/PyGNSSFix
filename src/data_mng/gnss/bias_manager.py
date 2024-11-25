"""
Satellite Code and Phase Biases Manager Module
This module implements the Satellite Biases manager, that manages the computation of satellite code and phase bias
corrections in the reconstruction of GNSS observations.
"""
from collections import OrderedDict

from src.data_types.date import Epoch
from src.data_types.gnss import Satellite, DataType, data_type_from_rinex, get_data_type
from src.data_mng import Container
from src.io.rinex_parser import SinexBiasReader
from src.common_log import IO_LOG, get_logger
from src.io.config import EnumSatelliteBias, config_dict
from src import constants

__all__ = ["BiasManager", "BiasEntry"]

from src.models.gnss_models import get_bgd_correction

E1 = get_data_type("E1", "GAL")

_warning_cache = set()

class BiasEntry(Container):
    """
    `BiasEntry` class, inherits from the `Container` class.
    This class is a container that stores a valid satellite code or phase bias entry from a Bias SINEX file.

    Attributes:
        epoch_start(Epoch): Start time for the bias estimate.
        epoch_end(Epoch): End time for the bias estimate.
        bias_type(str): Bias type identifier. Available types are DSB or OSB.
        sat(Satellite): Satellite for the bias estimate
        service_list(tuple[str,str]): Observables used for estimating the biases.
            The OBS2 field remains blank in case of absolute (OSB) estimates.
        value(float): Estimated (offset) value of the bias parameter.
        units(str): Bias estimates are given in the specified unit. Unit has to be `ns` (nanoseconds) for code
            biases; phase biases may be given in `cyc` (cycles).
    """
    __slots__ = ["epoch_start", "epoch_end", "bias_type", "sat", "service_list", "datatype_list", "value", "units"]

    def __init__(self, epoch_start: Epoch, epoch_end: Epoch, bias_type: str, sat: Satellite,
                 service_list: tuple, value: float, units: str):
        """ Base Constructor with no arguments. The attributes are filled in the `compute` method.  """
        super().__init__()
        # Parameter type check
        if not isinstance(epoch_start, Epoch):
            raise AttributeError(f"Parameter 'epoch_start' must be of type Epoch. Type {type(epoch_start)} was "
                                 f"provided instead")
        if not isinstance(epoch_end, Epoch):
            raise AttributeError(f"Parameter 'epoch_end' must be of type Epoch. Type {type(epoch_end)} was "
                                 f"provided instead")
        if not isinstance(bias_type, str):
            raise AttributeError(f"Parameter 'bias_type' must be of type string. Type {type(bias_type)} was "
                                 f"provided instead")
        if not isinstance(sat, Satellite):
            raise AttributeError(f"Parameter 'sat' must be of type Satellite. Type {type(sat)} was "
                                 f"provided instead")
        if not isinstance(service_list, tuple):
            raise AttributeError(f"Parameter 'service_list' must be of type tuple. Type {type(service_list)} was "
                                 f"provided instead")
        if not isinstance(value, float):
            raise AttributeError(f"Parameter 'value' must be of type float. Type {type(value)} was "
                                 f"provided instead")
        if not isinstance(units, str):
            raise AttributeError(f"Parameter 'units' must be of type Satellite. Type {type(units)} was "
                                 f"provided instead")

        if units not in ("ns", "cyc"):
            raise AttributeError(f"Unknown units {units}. Legal values are 'ns' or 'cyc'.")
        if bias_type not in ("OSB", "DSB"):
            raise AttributeError(f"Unknown bias type {bias_type}. Legal values are 'OSB' or 'DSB'.")

        self.epoch_start = epoch_start
        self.epoch_end = epoch_end
        self.service_list = service_list
        self.sat = sat
        self.units = units
        self.value = value
        self.bias_type = bias_type
        self.datatype_list = list()
        for service in self.service_list:
            self.datatype_list.append(data_type_from_rinex(service, sat.sat_system))



class BiasManager:
    """
    Satellite Code and Phase Biases Manager class.
    This class stores and manages satellite biases. There are three types of biases (from three products) considered:
        * Broadcast biases: TGDs for GPS and BGDs for GAL
        * Differential Code Biases (DCBs): Bias-SINEX files with DSB (Relative) biases
        * Observable-Specific Biases (OSBs): Bias-SINEX files with OSB (Absolute) biases

    For DCBs and OSBs, in the current implementation, it is assumed that we have 1 Bias value available per satellite,
    and per observation (or observation differential), and it is valid for the full simulation.
    That is, we do not have to deal with a TimeSeries of biases with different validity intervals.

    Attributes:
        dcb_data(OrderedDict): attribute that stores the DCB bias products. It is an :py:class:`OrderedDict`
            defined as:
                * keys -> :py:class:`Satellite` instances.
                * values -> :py:class:`list` list of :py:class:`BiasEntry` objects for this satellite.
        nav_data(src.data_mng.gnss.navigation_data.NavigationData): attribute that stores the broadcast navigation data
        bias_enum(EnumSatelliteBias): Enumeration that selects which type of bias is configured by the user
    """

    def __init__(self):
        """ Constructor of the `BiasManager` instance """
        self.dcb_data = OrderedDict()  # internal dict to store DCB data
        self.osb_data = OrderedDict()  # internal dict to store OSB data
        self.nav_data = None  # contains the broadcast navigation data
        self.bias_enum = None  # enumeration with the bias type configured by the user


    def init(self, nav_data, bias_files, bias_enum: EnumSatelliteBias):
        """
        Initialize this BiasManager instance.
        The enumeration in `bias_enum` selects the active bias type to be considered.

        Args:
            nav_data(src.data_mng.gnss.navigation_data.NavigationData or None): the `NavigationData`
                with navigation data (must be already initialized)
            bias_files(list): list of Bias SINEX files from the user configuration
            bias_enum(EnumSatelliteBias): Enumeration that selects which type of bias is configured by the user
        """
        self.nav_data = nav_data
        self.bias_enum = bias_enum

        log = get_logger(IO_LOG)
        if self.bias_enum in (EnumSatelliteBias.DCB, EnumSatelliteBias.OSB):
            log.info("Launching Satellite Code and Phase Bias Reader")
            for file in bias_files:
                SinexBiasReader(file, self, self.bias_enum)

    def __str__(self):
        """ Print the bias data to a string, for debug purposes """
        if self.bias_enum == EnumSatelliteBias.BROADCAST:
            return f"Code bias from broadcast ephemerides (TGD and BGD)"
        elif self.bias_enum == EnumSatelliteBias.DCB:
            myStr = "Differential Code Bias (DCB) Data:"
            data = self.dcb_data
        else:
            myStr = "Observable Specific Code and Phase Bias (OSB) Data:"
            data = self.osb_data
        for sat, bias_list in data.items():
            myStr += f"Satellite {sat}:\n"
            for bias in bias_list:
                myStr += f"\t{repr(bias)}\n"
        return myStr

    def add_bias(self, epoch_start: Epoch, epoch_end: Epoch, bias_type: str, sat: Satellite,
                 service_list: tuple, value: float, units: str):
        """
        Method to set a satellite bias from the SINEX file (either DCB or OSB) for the given satellite

        Args:
            epoch_start(Epoch): Start time for the bias estimate.
            epoch_end(Epoch): End time for the bias estimate.
            bias_type(str): Bias type identifier. Available types are DSB or OSB.
            sat(Satellite): Satellite for the bias estimate
            service_list(tuple[str,str]): Observables used for estimating the biases.
                The OBS2 field remains blank in case of absolute (OSB) estimates.
            value(float): Estimated (offset) value of the bias parameter.
            units(str): Bias estimates are given in the specified unit. Unit has to be `ns` (nanoseconds) for code
                biases; phase biases may be given in `cyc` (cycles).

        Raises:
            AttributeError: this exception is raised if one of the input arguments is not of the valid type.
        """

        bias_entry = BiasEntry(epoch_start, epoch_end, bias_type, sat, service_list, value, units)

        if bias_type == "DSB":
            if sat not in self.dcb_data:
                self.dcb_data[sat] = list()
            self.dcb_data[sat].append(bias_entry)
        elif bias_type == "OSB":
            if sat not in self.osb_data:
                self.osb_data[sat] = list()
            self.osb_data[sat].append(bias_entry)

    # Getters
    def get_dcb_data(self) -> OrderedDict:
        """ Returns the precise DCB products. """
        return self.dcb_data

    def get_osb_data(self) -> OrderedDict:
        """ Returns the precise OSB products. """
        return self.dcb_data

    def get_nav_message(self, sat, epoch):
        """
        Fetches the closest navigation message valid for the satellite and epoch.

        Args:
            sat(src.data_types.gnss.Satellite): satellite to fetch the data
            epoch(src.data_types.date.Epoch): epoch to fetch the data
        Returns:
            src.data_mng.gnss.navigation_data.NavigationPoint: navigation message

        Raises an exception if the provided satellite does not have valid data
        """
        return self.nav_data.get_closest_message(sat, epoch)

    def bias_correction(self, epoch, sat, datatype):
        bias = 0
        if self.bias_enum == EnumSatelliteBias.BROADCAST:
            bias = self._bias_correction_broadcast(epoch, sat, datatype)
        elif self.bias_enum == EnumSatelliteBias.DCB:
            bias = self._bias_correction_dcb(epoch, sat, datatype)
        elif self.bias_enum == EnumSatelliteBias.OSB:
            pass
        return bias

    def _bias_correction_broadcast(self, epoch, sat, datatype):
        nav_message = self.get_nav_message(sat, epoch)
        return get_bgd_correction(datatype, nav_message)

    def _bias_correction_dcb(self, epoch, sat:Satellite, datatype:DataType):
        print("getting DCB for", epoch, sat, datatype)
        if sat not in self.dcb_data:
            if sat not in _warning_cache:
                print("warning: no DCB data available for satellite", sat)  # Add as log message.
                _warning_cache.add(sat)

        dcb_list = self.dcb_data[sat]
        #for dcb in dcb_list:
        #    print(repr(dcb))
        constellation = sat.sat_system
        user_service = config_dict.get_services()[constellation]['user_service']
        ionofree_service = config_dict.get_services()[constellation]['clock_product_service']
        # print("services", config_dict.get_services()[constellation])

        # NOTA: isto pode ser feito previamente, na construção do BiasManager.
        for dcb in dcb_list:
            if BiasManager.check_list(dcb.service_list, ionofree_service):
                if_dcb = dcb
                break
        else:
            print(f"ERROR: could not find DCB for the iono-free combination {ionofree_service}")
            exit()
        #print("dcb entry for iono free combination", ionofree_service, " is ", repr(if_dcb))

        # find the user service associated with the provided datatype
        for service in user_service:
            if str(datatype.freq_number) == service[0]:
                print("-> datatype", datatype, "is associated to matched service", service)
                break
        else:
            print("ERROR: the provided datatype does not match any user service (should never happen...)")
            exit(-1)
        #print("ionofree_service", ionofree_service)
        if service in ionofree_service:
            print("user service is part of iono-free combination! Only 1 DCB is required")
            # fetch required DCB
            # print(f"Fix {datatype} of service {service} with DCB {repr(if_dcb)}")
            return BiasManager.get_dcb_for_service(if_dcb, service)
        else:
            print("user service is NOT part of iono-free combination. DCB convertion is required")
            for dcb in dcb_list:
                if BiasManager.check_list(dcb.service_list, user_service):
                    user_dcb = dcb
                    break
            else:
                print(f"ERROR: could not find DCB for the iono-free combination {ionofree_service}")
                exit()

            bias = BiasManager.correct_service(if_dcb, user_dcb)
            return bias + BiasManager.get_dcb_for_service(user_dcb, service)
        #base_freq =
        return 0

    @staticmethod
    def check_list(service_list, ionofree_service):
        # returns True if the two observation lists match (also in code attribute)
        if len(service_list) != 2 or len(ionofree_service) != 2:
            raise "error"
        for service in ionofree_service:
            if f"C{service}" not in service_list:
                return False
        return True

    @staticmethod
    def correct_service(if_dcb, user_dcb):
        print(f"correcting bias. IF dcb = {repr(if_dcb)}. User dcb = {repr(user_dcb)}")
        if_datatypes = if_dcb.datatype_list
        user_datatypes = user_dcb.datatype_list
        mu_base = (E1.freq_value / if_datatypes[0].freq.freq_value) ** 2  # this is always 1.
        mu_second = (E1.freq_value / if_datatypes[1].freq.freq_value) ** 2

        mu_user = (E1.freq_value / user_datatypes[1].freq.freq_value) ** 2

        return mu_base / (mu_second - mu_base) * if_dcb.value * 1E-9 - mu_base / (mu_user - mu_base) * user_dcb.value * 1E-9

    @staticmethod
    def get_dcb_for_service(dcb, service):
        # TODO: check signal (+/-)
        # If we reach here, the provided service **must** be one of the services of the dcb
        print(f"Fix datatype of service {service} with DCB {repr(dcb)}")
        datatypes = dcb.datatype_list
        mu1 = (E1.freq_value / datatypes[0].freq.freq_value) ** 2  # this is always 1.
        mu2 = (E1.freq_value / datatypes[1].freq.freq_value) ** 2
        id = dcb.service_list.index(f"C{service}")
        if id == 0:
            scale = mu1 / (mu2 - mu1)
        else:
            scale = mu2 / (mu2 - mu1)
        return scale * dcb.value * 1E-9
