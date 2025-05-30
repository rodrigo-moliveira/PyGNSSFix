"""
Satellite Code and Phase Biases Manager Module
This module implements the Satellite Biases manager, that manages the computation of satellite code and phase bias
corrections in the reconstruction of GNSS observations.
"""
from collections import OrderedDict

from src.data_types.date import Epoch
from src.data_types.gnss import Satellite, DataType, data_type_from_rinex
from src.data_mng import Container
from src.errors import ReconstructionError
from src.io.rinex_parser import SinexBiasReader
from src.common_log import IO_LOG, get_logger
from src.io.config import EnumSatelliteBias, config_dict

__all__ = ["BiasManager", "BiasEntry"]

from src.models.gnss_models import get_bgd_correction


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
        osb_data(OrderedDict): attribute that stores the OSB bias products. It is an :py:class:`OrderedDict`
            defined as:
                * keys -> :py:class:`Satellite` instances.
                * values -> :py:class:`list` list of :py:class:`BiasEntry` objects for this satellite.
        nav_data(src.data_mng.gnss.navigation_data.NavigationData): attribute that stores the broadcast navigation data
        bias_enum(EnumSatelliteBias): Enumeration that selects which type of bias is configured by the user
        log(logging.Logger): logger instance
        _cache(dict): cache to store internal data, to optimize consecutive calls
    """

    def __init__(self):
        """ Constructor of the `BiasManager` instance """
        self.dcb_data = OrderedDict()  # internal dict to store DCB data
        self.osb_data = OrderedDict()  # internal dict to store OSB data
        self.nav_data = None  # contains the broadcast navigation data
        self.bias_enum = None  # enumeration with the bias type configured by the user
        self.log = get_logger(IO_LOG)
        self._cache = {"dcb_if": dict(), "dcb_atr_fix": dict(), "final_bias_for_datatype": dict()}

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
        """ Print the bias data to a string, for debug purposes. """
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
        Method to set a satellite bias from the SINEX file (either DCB or OSB) for the given satellite.

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
        and epoch.

        Args:
            epoch(Epoch): epoch to fetch the bias correction
            sat(Satellite): satellite to fetch the bias correction
            datatype(DataType): datatype of the PR/CP observation to be reconstructed

        Returns:
             float: returns the computed satellite bias correction in seconds, to be applied in the PR/CP
                observation reconstruction.

        Raises:
            ReconstructionError: if the bias correction was not computed successfully, an exception is raised.
        """
        bias = 0
        if self.bias_enum == EnumSatelliteBias.BROADCAST:
            # if the datatype is carrier, convert it to code
            if DataType.is_carrier(datatype):
                datatype = DataType.carrier_to_code(datatype)

            # no cache here -> always use the closest navigation message
            bias = self._bias_correction_broadcast(epoch, sat, datatype)
        elif self.bias_enum == EnumSatelliteBias.DCB:
            # if the datatype is carrier, convert it to code
            if DataType.is_carrier(datatype):
                datatype = DataType.carrier_to_code(datatype)

            if (sat, datatype) not in self._cache["final_bias_for_datatype"]:
                _bias = -self._bias_correction_dcb(sat, datatype)
                self._cache["final_bias_for_datatype"][(sat, datatype)] = _bias
                self.log.info(f"Satellite Hardware Bias for {sat} and {datatype}: {_bias} seconds.")
            bias = self._cache["final_bias_for_datatype"][(sat, datatype)]
        elif self.bias_enum == EnumSatelliteBias.OSB:
            if (sat, datatype) not in self._cache["final_bias_for_datatype"]:
                _bias = self._bias_correction_osb(sat, datatype)
                self._cache["final_bias_for_datatype"][(sat, datatype)] = _bias
                self.log.info(f"Satellite Hardware Bias for {sat} and {datatype}: {_bias} seconds.")
            bias = self._cache["final_bias_for_datatype"][(sat, datatype)]
        return bias

    def _bias_correction_broadcast(self, epoch: Epoch, sat: Satellite, datatype: DataType) -> float:
        """ Get the satellite hardware code bias correction for the provided satellite and datatype
        (broadcast BGDs/TGDs). """
        nav_message = self.get_nav_message(sat, epoch)
        return get_bgd_correction(datatype, nav_message)

    def _bias_correction_dcb(self, sat: Satellite, datatype: DataType) -> float:
        """ Get the satellite hardware code bias correction for the provided satellite and datatype using the
        DCB products. Handles both combined (iono-free) and uncombined services.

        For GPS satellites, it also handles the attribute corrections. For example, if the user service is C1C and the
        clock product is C1W-C2W, the attribute of the clock product must be corrected to match the user service.
        """
        if sat not in self.dcb_data:
            self.log.warning(f"No DCB data available for satellite {sat}. Setting satellite bias to 0.")
            return 0.0

        dcb_list = self.dcb_data[sat]
        constellation = sat.sat_system
        ionofree_service = config_dict.get_services()[constellation]['clock_product_service']

        # Compute the bias correction for uncombined services
        if not DataType.is_iono_free(datatype):
            if sat not in self._cache["dcb_if"]:
                dcb_if = self._find_dcb_in_list(ionofree_service, dcb_list)
                if dcb_if is not None:
                    self._cache["dcb_if"][sat] = dcb_if
                    self._cache["dcb_atr_fix"][sat] = 0.0
                else:
                    if constellation != "GPS":
                        self.log.warning(f"Could not find DCB data for satellite {sat} and clock product "
                                         f"combination {ionofree_service}. Setting satellite bias to 0.")
                        return 0.0

                    # in case of GPS, attempt to correct the attribute of the clock product
                    success = self._correct_if_clock_attribute(ionofree_service, dcb_list, sat)
                    if not success:
                        self.log.warning(f"Could not fix iono-free clock combination {ionofree_service} "
                                         f"according to the available DCB data for satellite {sat}. "
                                         f"Setting satellite bias to 0.")
                        return 0.0

            fix_atr = self._cache["dcb_atr_fix"][sat]  # zero if not applicable
            return fix_atr + self._bias_correction_dcb_uncombined(dcb_list, sat, datatype)

        # Compute the bias correction for combined iono-free services
        else:
            return self._bias_correction_dcb_combined(dcb_list, ionofree_service, sat, datatype)

    def _bias_correction_dcb_uncombined(self, dcb_list: list[BiasEntry], sat: Satellite, datatype: DataType) -> float:
        """ Get the satellite hardware bias correction for the provided satellite and uncombined datatype. """
        dcb_if = self._cache["dcb_if"][sat]
        ionofree_service = dcb_if.service_list
        base_if_obs = ionofree_service[0]
        second_if_obs = ionofree_service[1]
        constellation = sat.sat_system

        # find the user service associated with the provided datatype
        service = self._get_service_for_datatype(datatype, constellation)

        if service in ionofree_service:
            self.log.info(f"User service {service} (datatype {datatype} for {sat}) is part of the clock product "
                          f"combination {ionofree_service}.")

            # fetch required bias for this service
            out_bias = BiasManager._get_code_bias_for_service_dcb(dcb_if, service)

        elif (service[1] == base_if_obs[1] and service[2] != base_if_obs[2]) or (
                service[1] == second_if_obs[1] and service[2] != second_if_obs[2]):
            if service[1] == base_if_obs[1]:
                correction_service = [base_if_obs, service]
                obs_index = 0
            else:
                correction_service = [second_if_obs, service]
                obs_index = 1

            self.log.info(f"User service {service} (datatype {datatype} for {sat}) does not match the clock product "
                          f"combination {ionofree_service}. Convertion of the code attribute of the "
                          f"{'first' if obs_index == 0 else 'second'} frequency is required using DCB "
                          f"{correction_service}.")

            correction_dcb = self._find_dcb_in_list(correction_service, dcb_list)
            if correction_dcb is None:
                self.log.warning(f"Could not find DCB data for satellite {sat} and combination {correction_service} "
                                 f"to fix the code attribute. Setting the attribute fix to 0.")
                atr_fix = 0.0
            else:
                if correction_dcb.service_list[0] == service:
                    atr_fix = -correction_dcb.value * 1E-9
                else:
                    atr_fix = correction_dcb.value * 1E-9

            out_bias = atr_fix + BiasManager._get_code_bias_for_service_dcb(
                dcb_if, base_if_obs if obs_index == 0 else second_if_obs)

        elif service[1] != second_if_obs[1]:
            # in this case, we need to change the 2nd frequency of the clock combination
            correction_service = [base_if_obs, service]
            self.log.info(f"User service {service} (datatype {datatype} for {sat}) does not match the clock product "
                          f"combination {ionofree_service}. Bias convertion is required from {ionofree_service} to "
                          f"{correction_service} (change of second frequency).")

            correction_dcb = self._find_dcb_in_list(correction_service, dcb_list)
            if correction_dcb is None:
                self.log.warning(f"Could not find DCB data for satellite {sat} and "
                                 f"service {correction_service}. Setting satellite bias to 0.")
                return 0.0

            correction_bias = BiasManager._correct_if_clock_for_dcb(dcb_if, correction_dcb)
            out_bias = correction_bias + BiasManager._get_code_bias_for_service_dcb(correction_dcb, service)
        else:
            raise NotImplementedError(f"The base observation for the clock service {ionofree_service} and the "
                                      f"combined user service {service} do not match. "
                                      f"This case is not yet implemented in the software.")
        return out_bias

    def _bias_correction_dcb_combined(self, dcb_list: list[BiasEntry], ionofree_service: tuple[str, str],
                                      sat: Satellite, datatype: DataType) -> float:
        """ Get the satellite hardware bias correction for the provided satellite and combined datatype. """
        constellation = sat.sat_system

        user_services = [f"C{service}" for service in config_dict.get_services()[constellation]['user_service']]
        if len(user_services) != 2 or len(ionofree_service) != 2:
            raise ReconstructionError(f"Error finding bias correction for datatype {datatype}: Clock product service "
                                      f"{ionofree_service} and user service {user_services} must have 2 observations.")

        if ionofree_service[0][1] != user_services[0][1]:
            raise ReconstructionError(f"Error finding bias correction for datatype {datatype}: Base frequency of the "
                                      f"clock product service {ionofree_service} does not match the base frequency "
                                      f"of the user service {user_services}.")

        # check if the user service is the same as the clock product
        if ionofree_service[0] == user_services[0] and ionofree_service[1] == user_services[1]:
            self.log.info(f"User combined service {user_services} matches the clock product (iono-free) "
                          f"{ionofree_service}. No bias correction is required for sat {sat} and datatype {datatype}.")
            return 0.0

        # Check if we need to change the second frequency (ignoring the code attribute)
        if ionofree_service[1][1] != user_services[1][1]:
            fix = 0.0
            # Need to fix the second frequency
            self.log.info(f"The combined user service {user_services} (datatype {datatype} for {sat}) does not match "
                          f"the clock service iono-free combination {ionofree_service}. Bias convertion is required "
                          f"from {ionofree_service} to {user_services}.")

            dcb_if = self._find_dcb_in_list(ionofree_service, dcb_list)
            if dcb_if is None:
                dcb_if = self._find_dcb_in_list(ionofree_service, dcb_list, ignore_attribute=True)
                if dcb_if.service_list[0][2] != ionofree_service[0][2]:
                    # fix the code attribute of the first frequency
                    fix = fix + self._fix_service_attribute(dcb_list, ionofree_service, dcb_if.service_list, 0,
                                                            constellation)
                if dcb_if.service_list[1][2] != ionofree_service[1][2]:
                    # fix the code attribute of the second frequency
                    fix = fix + self._fix_service_attribute(dcb_list, ionofree_service, dcb_if.service_list, 1,
                                                            constellation)

            dcb_user = self._find_dcb_in_list(user_services, dcb_list)
            if dcb_user is None:
                dcb_user = self._find_dcb_in_list([dcb_if.service_list[0], user_services[1]], dcb_list)
                if dcb_user is None:
                    dcb_user = self._find_dcb_in_list(user_services, dcb_list, ignore_attribute=True)

                if dcb_user.service_list[0][2] != user_services[0][2]:
                    # fix the code attribute of the first frequency
                    fix = fix + self._fix_service_attribute(dcb_list, dcb_user.service_list, user_services, 0,
                                                            constellation)
                if dcb_user.service_list[1][2] != user_services[1][2]:
                    # fix the code attribute of the second frequency
                    fix = fix + self._fix_service_attribute(dcb_list, dcb_user.service_list, user_services, 1,
                                                            constellation)

            if dcb_user is None or dcb_if is None:
                self.log.warning(f"Could not find DCB data for satellite {sat} and "
                                 f"services {dcb_user}, {ionofree_service}. Please check the provided "
                                 f"DCB file. Setting satellite bias to 0.")
                return 0.0

            # first apply DCB correction to the iono-free clock product
            out_bias = BiasManager._correct_if_clock_for_dcb(dcb_if, dcb_user) + fix

        # fix the code attribute only
        else:
            self.log.info(f"Fixing the code attribute from the clock product service {ionofree_service} to the user "
                          f"service {user_services} (datatype {datatype} for {sat}).")
            fix = 0.0
            if ionofree_service[0][2] != user_services[0][2]:
                # fix the code attribute of the first frequency
                fix = fix + self._fix_service_attribute(dcb_list, ionofree_service, user_services, 0, constellation)
            if ionofree_service[1][2] != user_services[1][2]:
                # fix the code attribute of the second frequency
                fix = fix + self._fix_service_attribute(dcb_list, ionofree_service, user_services, 1, constellation)

            out_bias = fix

        return out_bias

    def _find_dcb_in_list(self, service_to_find: tuple[str, str], dcb_list: list[BiasEntry],
                          ignore_attribute=False) -> BiasEntry:
        """ Find the `BiasEntry` instance from the provided dcb_list that matches the required service. """
        for dcb in dcb_list:
            if not ignore_attribute:
                if self._check_service_lists(dcb.service_list, service_to_find):
                    return dcb
            else:
                if self._check_service_lists_ignore_attribute(dcb.service_list, service_to_find):
                    return dcb
        return None

    @staticmethod
    def _get_service_for_datatype(datatype: DataType, constellation: str) -> str:
        """ Finds the service (RINEX notation) that matches the provided datatype and constellation. """
        user_services = config_dict.get_services()[constellation]['user_service']
        for service in user_services:
            if str(datatype.freq_number) == service[0]:
                found_service = service
                break
        else:
            # This exception should never happen...
            raise ReconstructionError(
                f"The provided datatype {datatype} does not match any user service: {user_services}. ")

        if DataType.is_code(datatype):
            return f"C{found_service}"
        elif DataType.is_carrier(datatype):
            return f"L{found_service}"
        else:
            raise ReconstructionError(f"The provided datatype {datatype} is not a code or carrier type. ")

    @staticmethod
    def _check_service_lists(service_dcb: tuple[str, str], service_to_find: tuple[str, str]) -> bool:
        """ Returns True if the two observation lists match (also in code attribute). """
        if len(service_dcb) != 2 or len(service_to_find) != 2:
            raise ReconstructionError(f"DCB service {service_dcb} and service to find {service_to_find}"
                                      f" must have length 2.")
        for service in service_to_find:
            if service not in service_dcb:
                return False
        return True

    @staticmethod
    def _check_service_lists_ignore_attribute(service1: tuple[str, str], service2: tuple[str, str]) -> bool:
        """ Returns True if the two observation lists match (ignoring the code attribute). """
        if len(service1) != 2 or len(service2) != 2:
            raise ReconstructionError(f"DCB service {service1} and service to find {service2} must have length 2.")

        freq1_1 = service1[0][1]
        freq1_2 = service1[1][1]
        freq2_1 = service2[0][1]
        freq2_2 = service2[1][1]

        if freq1_1 == freq2_1 and freq1_2 == freq2_2 or freq1_1 == freq2_2 and freq1_2 == freq2_1:
            return True
        return False

    def _correct_if_clock_attribute(self, ionofree_service: tuple[str, str], dcb_list: list[BiasEntry],
                                    sat: Satellite) -> bool:
        """ Corrects the attribute of the iono-free clock product using the available DCB data. """
        if sat.sat_system != "GPS":
            self.log.warning(f"Attribute correction is only available for GPS satellites. Satellite {sat} is not GPS.")
            return False

        dcb_if = self._find_dcb_in_list(ionofree_service, dcb_list, ignore_attribute=True)
        fix_atr = 0.0
        if dcb_if is None:
            self.log.warning(f"Could not find any DCB data for satellite {sat} and clock product combination "
                             f"{ionofree_service} (ignoring attributes). Setting satellite bias to 0.")
            return False
        if dcb_if.service_list[0] != ionofree_service[0]:
            fix_atr = fix_atr + self._fix_service_attribute(dcb_list, ionofree_service, dcb_if.service_list, 0,
                                                            sat.sat_system)
        if dcb_if.service_list[1] != ionofree_service[1]:
            fix_atr = fix_atr + self._fix_service_attribute(dcb_list, ionofree_service, dcb_if.service_list, 1,
                                                            sat.sat_system)

        self.log.info(
            f"Fixed iono free clock combination from {ionofree_service} to {dcb_if.service_list}.")
        self.log.debug(f"Attribute fix is {fix_atr} seconds.")

        self._cache["dcb_if"][sat] = dcb_if
        self._cache["dcb_atr_fix"][sat] = fix_atr
        return True

    def _fix_service_attribute(self, dcb_list: list[BiasEntry], ionofree_service: tuple[str, str],
                               user_service: tuple[str, str], index: int, constellation: str) -> float:
        """ Fix the attribute of the provided combined service using the available DCB data. """
        fix_service = [ionofree_service[index], user_service[index]]
        dcb_atr_fix = self._find_dcb_in_list(fix_service, dcb_list)
        if dcb_atr_fix is None:
            self.log.warning(f"Could not find DCB data for satellite to fix the code attribute of the clock product "
                             f"from {ionofree_service} to {user_service}. Setting the attribute fix to 0.")
            return 0.0
        else:
            self.log.info(f"Attribute correction for the clock product combination {ionofree_service} to "
                          f"{user_service}. Using DCB {dcb_atr_fix.service_list} to fix frequency index {index}.")

        first_obs = data_type_from_rinex(ionofree_service[0], constellation)
        second_obs = data_type_from_rinex(ionofree_service[1], constellation)
        first_freq = first_obs.freq.freq_value
        second_freq = second_obs.freq.freq_value

        mu_base = 1
        mu_second = (first_freq / second_freq) ** 2
        if index == 0:
            k = mu_second / (mu_second - mu_base)
        else:
            k = -mu_base / (mu_second - mu_base)

        if fix_service[0] == dcb_atr_fix.service_list[0] and fix_service[1] == dcb_atr_fix.service_list[1]:
            signal = 1
        else:
            signal = -1

        return signal * k * dcb_atr_fix.value * 1E-9

    @staticmethod
    def _correct_if_clock_for_dcb(if_dcb: BiasEntry, user_dcb: BiasEntry) -> float:
        """ Obtains the bias correction to convert the service (iono-free) satellite clock to the user service using
        the provided user_dcb. """
        if_datatypes = if_dcb.datatype_list
        user_datatypes = user_dcb.datatype_list
        mu_base = 1
        mu_second = (if_datatypes[0].freq.freq_value / if_datatypes[1].freq.freq_value) ** 2

        mu_user = (if_datatypes[0].freq.freq_value / user_datatypes[1].freq.freq_value) ** 2

        return mu_base / (mu_second - mu_base) * if_dcb.value * 1E-9 - mu_base / (
                mu_user - mu_base) * user_dcb.value * 1E-9

    @staticmethod
    def _get_code_bias_for_service_dcb(dcb: BiasEntry, service: str) -> float:
        """ Gets the code bias for the provided service from the available DCB.
        The provided service **must** be one of the services of the dcb. """
        datatypes = dcb.datatype_list
        mu1 = 1  # this is always 1.
        mu2 = (datatypes[0].freq.freq_value / datatypes[1].freq.freq_value) ** 2
        _id = dcb.service_list.index(service)
        if _id == 0:
            scale = mu1 / (mu2 - mu1)
        else:
            scale = mu2 / (mu2 - mu1)
        return scale * dcb.value * 1E-9

    def _bias_correction_osb(self, sat, datatype):
        """ Get the satellite hardware code bias correction for the provided satellite and datatype using the
        OSB products. Handles both combined (iono-free) and uncombined services.

        Typically, the OSB products provided by CODE do not contain all services. They contain:
            * GPS: corrections for L1 and L2
            * GAL: corrections for E1 and E5
        So this method is only appropriate for these datatypes.
        """
        if sat not in self.osb_data:
            self.log.warning(f"No OSB data available for satellite {sat}. Setting satellite bias to 0.")
            return 0.0

        if not DataType.is_iono_free(datatype):
            return self._bias_correction_osb_uncombined(sat, datatype)
        else:
            return self._bias_correction_osb_combined(sat)

    def _bias_correction_osb_uncombined(self, sat, datatype):
        """ Get the satellite hardware bias correction for the provided satellite and uncombined datatype. """
        bias = 0.0
        found = False
        user_service = self._get_service_for_datatype(datatype, sat.sat_system)

        for osb in self.osb_data[sat]:
            osb_datatype = osb.datatype_list[0]
            osb_service = osb.service_list[0]
            if datatype == osb_datatype:
                if user_service == osb_service:
                    bias = osb.value * 1E-9
                    self.log.info(f"Fetched bias OSB {bias} [s] for {sat} and {datatype}.")
                    found = True
                    break

        if not found:
            self.log.warning(f"Could not find OSB data for satellite {sat} and service {user_service}. "
                             f"Setting satellite bias correction to 0.")

        return bias

    def _bias_correction_osb_combined(self, sat):
        """ Get the satellite hardware bias correction for the provided satellite and combined datatype. """
        user_services = [f"C{service}" for service in config_dict.get_services()[sat.sat_system]['user_service']]

        if len(user_services) != 2:
            raise ReconstructionError(f"Error finding bias correction for satellite {sat}: User service "
                                      f"{user_services} must have 2 observations.")
        service1 = user_services[0]
        service2 = user_services[1]

        bias1 = bias2 = f1 = f2 = bias_if = 0.0
        found1 = found2 = False

        # find the bias for the user services
        for osb in self.osb_data[sat]:
            osb_service = osb.service_list[0]
            osb_datatype = osb.datatype_list[0]

            # matched service 1
            if service1 == osb_service:
                bias1 = osb.value * 1E-9
                f1 = osb_datatype.freq.freq_value
                found1 = True

            # matched service 2
            elif service2 == osb_service:
                bias2 = osb.value * 1E-9
                f2 = osb_datatype.freq.freq_value
                found2 = True

        # compute the combined iono-free bias
        if found1 and found2:
            gama1 = f1 * f1 / (f1 * f1 - f2 * f2)
            gama2 = f2 * f2 / (f1 * f1 - f2 * f2)
            bias_if = gama1 * bias1 - gama2 * bias2
        else:
            self.log.warning(f"Could not find OSB data for satellite {sat} and services {user_services}. "
                             f"Setting satellite bias correction to 0.")

        return bias_if
