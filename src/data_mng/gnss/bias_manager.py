"""
Satellite Code and Phase Biases Manager Module
This module implements the Satellite Biases manager, that manages the computation of satellite code and phase bias
corrections in the reconstruction of GNSS observations.
"""
from collections import OrderedDict

from src.data_types.date import Epoch
from src.data_types.gnss import Satellite, DataType, data_type_from_rinex, get_data_type
from src.data_mng import Container
from src.errors import ReconstructionError
from src.io.rinex_parser import SinexBiasReader
from src.common_log import IO_LOG, get_logger
from src.io.config import EnumSatelliteBias, config_dict

__all__ = ["BiasManager", "BiasEntry"]

from src.models.gnss_models import get_bgd_correction

E1 = get_data_type("E1", "GAL")


class BiasEntry(Container):
    """
    `BiasEntry` class, inherits from the `Container` class.
    This class is a container that stores a valid satellite code or phase bias entry from a Bias SINEX file.

    Attributes:
        epoch_start(Epoch): Start time for the bias estimate.
        epoch_end(Epoch): End time for the bias estimate.
        bias_type(str): Bias type identifier. Available types are DSB or OSB.
        sat(Satellite): Satellite for the bias estimate
        service_list(tuple[str,str]): services used for estimating the biases.
            The OBS2 field remains blank in case of absolute (OSB) estimates.
        datatype_list(tuple[DataType,DataType]): Observables used for estimating the biases.
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
        osb_data(OrderedDict): ... to be filled
        nav_data(src.data_mng.gnss.navigation_data.NavigationData): attribute that stores the broadcast navigation data
        bias_enum(EnumSatelliteBias): Enumeration that selects which type of bias is configured by the user
        log(logging.Logger): logger instance
        _cache(dict): internal cache to store intermediate data, to optimize consecutive calls
    """

    def __init__(self):
        """ Constructor of the `BiasManager` instance """
        self.dcb_data = OrderedDict()  # internal dict to store DCB data
        self.osb_data = OrderedDict()  # internal dict to store OSB data
        self.nav_data = None  # contains the broadcast navigation data
        self.bias_enum = None  # enumeration with the bias type configured by the user
        self.log = get_logger(IO_LOG)
        self._cache = {"warnings": set(), "dcb_if": dict(), "final_bias_for_datatype": dict()}

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

        if self.bias_enum in (EnumSatelliteBias.DCB, EnumSatelliteBias.OSB):
            self.log.info(f"Launching Satellite Code and Phase Bias Reader: mode {self.bias_enum}")
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
        """
        Main function to fetch the satellite hardware code/phase bias correction for the provided satellite, datatype
        and epoch

        Args:
            epoch(Epoch): epoch to fetch the bias correction
            sat(Satellite): satellite to fetch the bias correction
            datatype(DataType): datatype of the PR/CP observation to be reconstructed

        Returns:
             float: returns the computed bias correction in seconds, to be applied in the PR/CP
                observation reconstruction.
        """
        bias = 0
        if self.bias_enum == EnumSatelliteBias.BROADCAST:
            # no cache here -> always use the closest navigation message
            bias = self._bias_correction_broadcast(epoch, sat, datatype)
        elif self.bias_enum == EnumSatelliteBias.DCB:
            if (sat, datatype) not in self._cache["final_bias_for_datatype"]:
                _bias = -self._bias_correction_dcb(sat, datatype)
                self._cache["final_bias_for_datatype"][(sat, datatype)] = _bias
                self.log.info(f"Satellite Hardware Bias for {sat} and {datatype}: {_bias} seconds.")
            bias = self._cache["final_bias_for_datatype"][(sat, datatype)]
        elif self.bias_enum == EnumSatelliteBias.OSB:
            pass
        return bias

    def _bias_correction_broadcast(self, epoch, sat, datatype):
        """ Get the satellite hardware code bias correction for the provided satellite and datatype
        (broadcast BGDs/TGDs). """
        nav_message = self.get_nav_message(sat, epoch)
        return get_bgd_correction(datatype, nav_message)

    def _bias_correction_dcb(self, sat: Satellite, datatype: DataType):
        """ Get the satellite hardware code bias correction for the provided satellite and datatype
        (precise DCB products). """
        if sat not in self.dcb_data:
            if sat not in self._cache["warnings"]:
                self.log.warning(f"No DCB data available for satellite {sat}. Setting satellite bias to 0.")
                self._cache["warnings"].add(sat)
            return 0.0

        dcb_list = self.dcb_data[sat]

        constellation = sat.sat_system
        ionofree_service = config_dict.get_services()[constellation]['clock_product_service']

        base_if_obs = ionofree_service[0]

        if sat not in self._cache["dcb_if"]:
            dcb_if = self._find_dcb_in_list(ionofree_service, dcb_list)
            if dcb_if is not None:
                self._cache["dcb_if"][sat] = dcb_if
            else:
                self.log.warning(f"Could not find DCB data for satellite {sat} and "
                                 f"iono-free combination {ionofree_service}. Setting satellite bias to 0.")
                self._cache["warnings"].add(sat)
                return 0.0

        dcb_if = self._cache["dcb_if"][sat]

        # find the user service associated with the provided datatype
        service = self._get_service_for_datatype(datatype, constellation)

        if service in ionofree_service:
            self.log.info(f"User service {service} (datatype {datatype} for {sat}) is part of iono-free combination "
                          f"{ionofree_service}. Only the iono-free DCB is required.")
            # fetch required bias for this service
            out_bias = BiasManager._get_code_bias_for_service_dcb(dcb_if, service)
        else:
            correction_service = [base_if_obs, service]
            self.log.info(f"User service {service} (datatype {datatype} for {sat}) is not part of iono-free "
                          f"combination {ionofree_service}. Bias convertion is required from {ionofree_service} to "
                          f"{correction_service}.")

            correction_dcb = self._find_dcb_in_list(correction_service, dcb_list)
            if correction_dcb is None:
                self.log.warning(f"Could not find DCB data for satellite {sat} and "
                                 f"service {correction_service}. Setting satellite bias to 0.0.")
                self._cache["warnings"].add(sat)
                return 0.0

            correction_bias = BiasManager._correct_if_clock_for_dcb(dcb_if, correction_dcb)
            out_bias = correction_bias + BiasManager._get_code_bias_for_service_dcb(correction_dcb, service)

        return out_bias

    def _find_dcb_in_list(self, service_to_find, dcb_list):
        """ Find BiasEntry instance from the provided dcb_list that matches the required service. """
        for dcb in dcb_list:
            if self._check_service_lists(dcb.service_list, service_to_find):
                return dcb
        return None

    @staticmethod
    def _get_service_for_datatype(datatype, constellation):
        """ Finds the service (RINEX notation) that matches the provided datatype and constellation. """
        # TODO: Esta funcao vai portar-se mal com atributos diferentes...
        user_services = config_dict.get_services()[constellation]['user_service']
        for service in user_services:
            if str(datatype.freq_number) == service[0]:
                found_service = service
                break
        else:
            # This exception should never happen...
            raise ReconstructionError(
                f"The provided datatype {datatype} does not match any user service: {user_services}. ")

        return found_service

    @staticmethod
    def _check_service_lists(service_dcb, service_to_find):
        """ Returns True if the two observation lists match (also in code attribute). """
        if len(service_dcb) != 2 or len(service_to_find) != 2:
            raise ReconstructionError(f"DCB service {service_dcb} and service to find {service_to_find}"
                                      f" must have length 2.")
        for service in service_to_find:
            if f"C{service}" not in service_dcb:
                return False
        return True

    @staticmethod
    def _correct_if_clock_for_dcb(if_dcb, user_dcb):
        """ Obtains the bias correction to convert the iono-free the satellite clock to the user service using the
        provided user_dcb. """
        print(f"[PRINT]: correcting bias. IF dcb = {repr(if_dcb)}. User dcb = {repr(user_dcb)}")
        if_datatypes = if_dcb.datatype_list
        user_datatypes = user_dcb.datatype_list
        mu_base = (E1.freq_value / if_datatypes[0].freq.freq_value) ** 2  # this is always 1.
        mu_second = (E1.freq_value / if_datatypes[1].freq.freq_value) ** 2

        mu_user = (E1.freq_value / user_datatypes[1].freq.freq_value) ** 2

        return mu_base / (mu_second - mu_base) * if_dcb.value * 1E-9 - mu_base / (
                    mu_user - mu_base) * user_dcb.value * 1E-9

    @staticmethod
    def _get_code_bias_for_service_dcb(dcb, service):
        """ Gets the code bias for the provided service from the available DCB. """
        # TODO: check signal (+/-)
        # If we reach here, the provided service **must** be one of the services of the dcb
        print(f"[PRINT]: Fix datatype of service {service} with DCB {repr(dcb)}")
        datatypes = dcb.datatype_list
        mu1 = (E1.freq_value / datatypes[0].freq.freq_value) ** 2  # this is always 1.
        mu2 = (E1.freq_value / datatypes[1].freq.freq_value) ** 2
        _id = dcb.service_list.index(f"C{service}")
        if _id == 0:
            scale = mu1 / (mu2 - mu1)
        else:
            scale = mu2 / (mu2 - mu1)
        return scale * dcb.value * 1E-9
