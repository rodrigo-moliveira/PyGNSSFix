from src import constants
from src.errors import DataTypeError

__all__ = ["get_data_type"]


class DataType:
    """
    Class DataType
    Represents a GNSS observable datatype (Carrier Phase, Pseudo Range, Signal or Doppler)

    Attributes
        ----------
        freq_number : int
            The frequency index (1, 2, 5, ...)
        freq_value : float
            The frequency value in Hz
        data_type : str
            Short name for this type
        description : str
            Long description of this data type
    """

    # Force DataType objects to immutable
    __slots__ = ["_freq_number", "_freq_value", "_freq", "_data_type", "_description", "_constellation"]

    def __init__(self, freq_number: int = None, freq_value: float = None,
                 data_type: str = None, description: str = None, freq=None, constellation=None):

        if freq_number is not None:
            super(DataType, self).__setattr__('_freq_number', freq_number)
        if freq_value is not None:
            super(DataType, self).__setattr__('_freq_value', freq_value)
        if freq is not None:
            super(DataType, self).__setattr__('_freq', freq)
        if constellation is not None:
            super(DataType, self).__setattr__('_constellation', constellation)
        if data_type is not None:
            super(DataType, self).__setattr__('_data_type', data_type)
        else:
            super(DataType, self).__setattr__('_data_type', "UN")  # undefined
        if description is not None:
            super(DataType, self).__setattr__('_description', description)
        else:
            super(DataType, self).__setattr__('_description', "Unknown Data Type")

    # Objects of this class are immutable!
    def __setattr__(self, name, value):
        """Prevent modification of attributes."""
        raise AttributeError('DataType objects are immutable and cannot be modified')

    def __repr__(self):
        return f"DataType({self.data_type}, {self.description})"

    # conditional operations for use of ´in´ keyword in lists
    def __eq__(self, other):
        if isinstance(other, str):
            other = get_data_type(other)
        return self.data_type == other.data_type

    def __ne__(self, other):
        return not self.data_type == other.data_type

    def __lt__(self, other):
        return self.freq_number < other.freq_number

    def __str__(self):
        return self.data_type

    @property
    def data_type(self):
        if hasattr(self, "_data_type"):
            return self._data_type
        else:
            raise AttributeError(f"The provided datatype {str(self)} has no attribute '_data_type'")

    @property
    def freq(self):
        if hasattr(self, "_freq"):
            return self._freq
        else:
            raise AttributeError(f"The provided datatype {str(self)} has no attribute 'freq'")

    @property
    def freq_value(self):
        if hasattr(self, "_freq_value"):
            return self._freq_value
        else:
            raise AttributeError(f"The provided datatype {str(self)} has no attribute 'freq_value'")

    @property
    def freq_number(self):
        if hasattr(self, "_freq_number"):
            return self._freq_number
        else:
            raise AttributeError(f"The provided datatype {str(self)} has no attribute 'freq_number'")

    @property
    def description(self):
        if hasattr(self, "_description"):
            return self._description
        else:
            raise AttributeError(f"The provided datatype {str(self)} has no attribute 'description'")

    @staticmethod
    def is_code(data_type):
        return data_type in cAvailableCodes

    @staticmethod
    def is_iono_free_code(data_type):
        return data_type in cAvailableIonoFreeCodes

    @staticmethod
    def is_iono_free_smooth_code(data_type):
        return data_type in cAvailableIonoFreeSmoothCodes

    @staticmethod
    def is_carrier(data_type):
        return data_type in cAvailableCarriers

    @staticmethod
    def is_signal(data_type):
        return data_type in cAvailableSignals

    @staticmethod
    def get_code_datatypes(datatype_list):
        list_out = []
        for obs in datatype_list:
            if obs in cAvailableCodes + cAvailableSmoothCodes + cAvailableIonoFreeSmoothCodes + cAvailableIonoFreeCodes:
                list_out.append(obs)

        return list_out

    @staticmethod
    def get_signal_datatypes(datatype_list):
        list_out = []
        for obs in datatype_list:
            if obs in cAvailableSignals:
                list_out.append(obs)

        return list_out

    @staticmethod
    def get_iono_free_pseudorange(datatype1, datatype2):
        """
        mix data types to get the associated iono free datatype
        get_iono_free_datatype(PR1, PR1) => PR12
        """
        if DataType.is_code(datatype1) and DataType.is_code(datatype2):
            index1 = datatype1.freq_number
            index2 = datatype2.freq_number
            return get_data_type(f"PR{min(index1, index2)}{max(index1, index2)}")

        # if DataType.is_carrier(datatype1) and DataType.is_carrier(datatype2):
        #    index1 = datatype1.freq_number
        #    index2 = datatype2.freq_number
        #    return get_data_type(f"L{min(index1, index2)}{max(index1, index2)}")

        raise DataTypeError(f"Unable to obtain iono free pseudorange from the provided arguments. "
                            f"Datatype1 and datatype2 must be both code."
                            f"datatype1 = {datatype1.data_type}, datatype2 = {datatype2.data_type}")

    @staticmethod
    def get_smooth_pseudorange(datatype):
        """
        mix data types to get the associated smooth pseudo-range (SPR) or smooth iono free (IFSPR) datatype
        get_smooth_datatype(PR1) => SPR1
        """
        if DataType.is_code(datatype):
            # datatype = PR1 -> return SPR1
            index = datatype.freq_number
            return get_data_type(f"SPR{index}")

        if DataType.is_iono_free_code(datatype):
            # datatype = PR12 -> return SPR12
            index = datatype.freq_number
            return get_data_type(f"SPR{index}")

        raise DataTypeError(f"Unable to create smooth pseudorange from "
                            f"data type {str(datatype)}")


# Default Data Types

###############
# Frequencies #
###############

# GPS Frequencies
L1 = DataType(data_type="L1/E1", description="Frequency L1 for GPS", freq_value=constants.GPS_L1_FREQ,
              constellation="GPS")
L2 = DataType(data_type="L2", description="Frequency L2 for GPS", freq_value=constants.GPS_L2_FREQ,
              constellation="GPS")
L5 = DataType(data_type="L5", description="Frequency L5 for GPS", freq_value=constants.GPS_L5_FREQ,
              constellation="GPS")

# GAL Frequencies
E1 = DataType(data_type="E1", description="Frequency E1 (GAL)", freq_value=constants.GAL_E1_FREQ,
              constellation="GAL")
E5a = DataType(data_type="E5a", description="Frequency E5a (GAL)", freq_value=constants.GAL_E5A_FREQ,
               constellation="GAL")
E5b = DataType(data_type="E5b", description="Frequency E5b (GAL)", freq_value=constants.GAL_E5B_FREQ,
               constellation="GPS")
E5ALTBOC = DataType(data_type="E5AltBOC", description="Frequency E5AltBOC (GAL)",
                    freq_value=constants.GAL_E5ALTBOC_FREQ, constellation="GAL")
E6 = DataType(data_type="E6", description="Frequency E6 (GAL)", freq_value=constants.GAL_E6_FREQ,
              constellation="GPS")

###################
# Raw Observables #
###################

# GPS Code, Carrier, Signal and DOPs
PR1_GPS = DataType(data_type="PR1", description="PseudoRange in L1 (GPS)", freq=L1, freq_number=1)
PR2_GPS = DataType(data_type="PR2", description="PseudoRange in L2 (GPS)", freq=L2, freq_number=2)
PR5_GPS = DataType(data_type="PR5", description="PseudoRange in L5 (GPS)", freq=L5, freq_number=5)

CP1_GPS = DataType(data_type="CP1", description="Carrier Phase in L1 (GPS)", freq=L1, freq_number=1)
CP2_GPS = DataType(data_type="CP2", description="Carrier Phase in L2 (GPS)", freq=L2, freq_number=2)
CP5_GPS = DataType(data_type="CP5", description="Carrier Phase in L5 (GPS)", freq=L5, freq_number=5)

S1_GPS = DataType(data_type="S1", description="Signal Strength in L1 (GPS)", freq=L1, freq_number=1)
S2_GPS = DataType(data_type="S2", description="Signal Strength in L2 (GPS)", freq=L2, freq_number=2)
S5_GPS = DataType(data_type="S5", description="Signal Strength in L5 (GPS)", freq=L5, freq_number=5)

D1_GPS = DataType(data_type="D1", description="Doppler in L1 (GPS)", freq=L1, freq_number=1)
D2_GPS = DataType(data_type="D2", description="Doppler in L2 (GPS)", freq=L2, freq_number=2)
D5_GPS = DataType(data_type="D5", description="Doppler in L5 (GPS)", freq=L5, freq_number=5)

# GAL Code, Carrier, Signal and DOPs
PR1_GAL = DataType(data_type="PR1", description="PseudoRange in E1 (GAL)", freq=E1, freq_number=1)
PR5_GAL = DataType(data_type="PR5", description="PseudoRange in E5a (GAL)", freq=E5a, freq_number=5)
PR7_GAL = DataType(data_type="PR7", description="PseudoRange in E5b (GAL)", freq=E5b, freq_number=7)
PR8_GAL = DataType(data_type="PR8", description="PseudoRange in E5AltBOC (GAL)", freq=E5ALTBOC, freq_number=8)
PR6_GAL = DataType(data_type="PR6", description="PseudoRange in E6 (GAL)", freq=E6, freq_number=6)

CP1_GAL = DataType(data_type="CP1", description="Carrier Phase in E1 (GAL)", freq=E1, freq_number=1)
CP5_GAL = DataType(data_type="CP5", description="Carrier Phase in E5a (GAL)", freq=E5a, freq_number=5)
CP7_GAL = DataType(data_type="CP7", description="Carrier Phase in E5b (GAL)", freq=E5b, freq_number=7)
CP8_GAL = DataType(data_type="CP8", description="Carrier Phase in E5AltBOC (GAL)", freq=E5ALTBOC, freq_number=8)
CP6_GAL = DataType(data_type="CP6", description="Carrier Phase in E6 (GAL)", freq=E6, freq_number=6)

S1_GAL = DataType(data_type="S1", description="Signal Strength in E1 (GAL)", freq=E1, freq_number=1)
S5_GAL = DataType(data_type="S5", description="Signal Strength in E5a (GAL)", freq=E5a, freq_number=5)
S7_GAL = DataType(data_type="S7", description="Signal Strength in E5b (GAL)", freq=E5b, freq_number=7)
S8_GAL = DataType(data_type="S8", description="Signal Strength in E5AltBOC (GAL)", freq=E5ALTBOC, freq_number=8)
S6_GAL = DataType(data_type="S6", description="Signal Strength in E6 (GAL)", freq=E6, freq_number=6)

D1_GAL = DataType(data_type="D1", description="Doppler in E1 (GAL)", freq=E1, freq_number=1)
D5_GAL = DataType(data_type="D5", description="Doppler in E5a (GAL)", freq=E5a, freq_number=5)
D7_GAL = DataType(data_type="D7", description="Doppler in E5b (GAL)", freq=E5b, freq_number=7)
D8_GAL = DataType(data_type="D8", description="Doppler in E5AltBOC (GAL)", freq=E5ALTBOC, freq_number=8)
D6_GAL = DataType(data_type="D6", description="Doppler in E6 (GAL)", freq=E6, freq_number=6)

#####################################
# Iono Free PseudoRange Observables #
#####################################
# NOTE: currently, only allow for combinations between L1/E1 and secondary frequencies (L1/E1 must be present)
# GPS
PR12_GPS = DataType(data_type="PR12", description="L1-L2 Iono-Free PseudoRange (GPS)", freq_number=12)
PR15_GPS = DataType(data_type="PR15", description="L1-L5 Iono-Free PseudoRange (GPS)", freq_number=15)

# GAL
PR15_GAL = DataType(data_type="PR15", description="E1-E5a Iono-Free PseudoRange (GAL)", freq_number=15)
PR17_GAL = DataType(data_type="PR17", description="E1-E5b Iono-Free PseudoRange (GAL)", freq_number=17)
PR18_GAL = DataType(data_type="PR18", description="E1-E5AltBOC Iono-Free PseudoRange (GAL)", freq_number=18)
PR16_GAL = DataType(data_type="PR16", description="E1-E6 Iono-Free PseudoRange (GAL)", freq_number=16)

##################################
# Smooth PseudoRange Observables #
##################################
# GPS
SPR1_GPS = DataType(data_type="SPR1", description="Smooth PseudoRange in L1 (GPS)", freq=L1, freq_number=1)
SPR2_GPS = DataType(data_type="SPR2", description="Smooth PseudoRange in L2 (GPS)", freq=L2, freq_number=2)
SPR5_GPS = DataType(data_type="SPR5", description="Smooth PseudoRange in L5 (GPS)", freq=L5, freq_number=5)

# GAL
SPR1_GAL = DataType(data_type="SPR1", description="Smooth PseudoRange in E1 (GAL)", freq=E1, freq_number=1)
SPR5_GAL = DataType(data_type="SPR5", description="Smooth PseudoRange in E5a (GAL)", freq=E5a, freq_number=5)
SPR7_GAL = DataType(data_type="SPR7", description="Smooth PseudoRange in E5b (GAL)", freq=E5b, freq_number=7)
SPR8_GAL = DataType(data_type="SPR8", description="Smooth PseudoRange in E5AltBOC (GAL)", freq=E5ALTBOC, freq_number=8)
SPR6_GAL = DataType(data_type="SPR6", description="Smooth PseudoRange in E6 (GAL)", freq=E6, freq_number=6)

################################
# Smooth Iono-Free Observables #
################################
# GPS
SPR12_GPS = DataType(data_type="SPR12", description="L1-L2 Iono-Free Smooth PseudoRange (GPS)", freq_number=12)
SPR15_GPS = DataType(data_type="SPR15", description="L1-L5 Iono-Free Smooth PseudoRange (GPS)", freq_number=15)

# GAL
SPR15_GAL = DataType(data_type="SPR15", description="E1-E5a Iono-Free Smooth PseudoRange (GAL)", freq_number=15)
SPR17_GAL = DataType(data_type="SPR17", description="E1-E5b Iono-Free Smooth PseudoRange (GAL)", freq_number=17)
SPR18_GAL = DataType(data_type="SPR18", description="E1-E5AltBOC Iono-Free Smooth PseudoRange (GAL)", freq_number=18)
SPR16_GAL = DataType(data_type="SPR16", description="E1-E6 Iono-Free Smooth PseudoRange (GAL)", freq_number=16)

# Unknown
UN = DataType(data_type="UN", description="Unknown inputs type")

##############
# Containers #
##############
cAvailableCodes = [PR1_GPS, PR2_GPS, PR5_GPS, PR1_GAL, PR5_GAL, PR6_GAL, PR7_GAL, PR8_GAL]
cAvailableSignals = [S1_GPS, S2_GPS, S5_GPS, S1_GAL, S5_GAL, S6_GAL, S7_GAL, S8_GAL]
cAvailableCarriers = [CP1_GPS, CP2_GPS, CP5_GPS, CP1_GAL, CP5_GAL, CP6_GAL, CP7_GAL, CP8_GAL]
cAvailableFrequencies = [L1, L2, L5, E1, E5a, E5b, E5ALTBOC, E6]
cAvailableSmoothCodes = [SPR1_GPS, SPR2_GPS, SPR5_GPS, SPR1_GAL, SPR5_GAL, SPR6_GAL, PR7_GAL, PR8_GAL]
cAvailableIonoFreeSmoothCodes = [SPR12_GPS, SPR15_GPS, SPR15_GAL, SPR16_GAL, SPR17_GAL, SPR18_GAL]
cAvailableIonoFreeCodes = [PR12_GPS, PR15_GPS, PR15_GAL, PR16_GAL, PR17_GAL, PR18_GAL]

cGPSObsSignals = {"C": {"1": PR1_GPS, "2": PR2_GPS, "5": PR5_GPS},
                  "L": {"1": CP1_GPS, "2": CP2_GPS, "5": CP5_GPS},
                  "S": {"1": S1_GPS, "2": S2_GPS, "5": S5_GPS}}
cGALObsSignals = {"C": {"1": PR1_GAL, "5": PR5_GAL, "6": PR6_GAL, "7": PR7_GAL, "8": PR8_GAL},
                  "L": {"1": CP1_GAL, "5": CP5_GAL, "6": CP6_GAL, "7": CP7_GAL, "8": CP8_GAL},
                  "S": {"1": S1_GAL, "5": S5_GAL, "6": S6_GAL, "7": S7_GAL, "8": S8_GAL}}


def get_data_type(datatype: str):
    """
    Factory for inputs types. Receives a string representing the type of the datatype and returns the associated
    DataType instance

    Args:
        datatype (str): string with the short descriptor of the datatype to fetch (ex: C1, L5,...)
    Return:
         DataType : returns the corresponding DataType instance
    """
    for container in [cAvailableCodes, cAvailableSignals, cAvailableFrequencies, cAvailableCarriers,
                      cAvailableSmoothCodes, cAvailableIonoFreeCodes, cAvailableIonoFreeSmoothCodes]:
        for _type in container:
            if _type.data_type == datatype:
                return _type
    return UN


def data_type_from_rinex(data_type: str, constellation: str):
    if constellation == "GPS":
        return cGPSObsSignals[data_type[0]][data_type[1]]
    elif constellation == "GAL":
        return cGALObsSignals[data_type[0]][data_type[1]]
    return UN
