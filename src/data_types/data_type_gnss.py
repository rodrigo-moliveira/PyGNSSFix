from ...math_utils.Constants import Constant


class DataType:
    """
    Class DataType
    Represents a GNSS observable datatype (CarrierPhase, PseudoRange, or Signal)

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
    __slots__ = ["freq_number", "freq_value", "freq", "data_type", "description"]

    def __init__(self, freq_number: int = None, freq_value: float = None,
                 data_type: str = None, description: str = None, freq=None):

        if freq_number is not None:
            super(DataType, self).__setattr__('freq_number', freq_number)
        if freq_value is not None:
            super(DataType, self).__setattr__('freq_value', freq_value)
        if freq is not None:
            super(DataType, self).__setattr__('freq', freq)
        if data_type is not None:
            super(DataType, self).__setattr__('data_type', data_type)
        else:
            super(DataType, self).__setattr__('data_type', "UN")  # undefined
        if description is not None:
            super(DataType, self).__setattr__('description', description)
        else:
            super(DataType, self).__setattr__('description', "Unknown Data Type")

    # Objects of this class are immutable!
    def __setattr__(self, name, value):
        """Prevent modification of attributes."""
        raise AttributeError('DataType objects are immutable and cannot be modified')

    def __str__(self):
        return self.data_type

    def __repr__(self):
        return f"DataType({self.data_type}, {self.description})"

    # conditional operations for use of ´in´ keyword in lists
    def __eq__(self, other):
        if isinstance(other, str):
            other = DataTypeFactory(other)
        return self.data_type == other.data_type

    def __ne__(self, other):
        return not self.data_type == other.data_type

    def __lt__(self, other):
        return self.freq_number < other.freq_number

    def get_short_type(self):
        return str(self)

    def get_frequency(self):
        if hasattr(self, "freq"):
            return self.freq
        elif hasattr(self, "freq_value"):
            return self.freq_value
        else:
            raise ValueError(f"The provided datatype {str(self)} has no attribute 'frequency'")

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
    def is_iono_free_carrier(data_type):
        return data_type in cAvailableIonoFreeCarriers

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
    def get_iono_free_datatype(datatype1, datatype2):
        """
        mix data types to get the associated iono free datatype
        """
        if DataType.is_code(datatype1) and DataType.is_code(datatype2):
            index1 = datatype1.freq_number
            index2 = datatype2.freq_number
            return DataTypeFactory(f"C{min(index1, index2)}{max(index1, index2)}")

        if DataType.is_carrier(datatype1) and DataType.is_carrier(datatype2):
            index1 = datatype1.freq_number
            index2 = datatype2.freq_number
            return DataTypeFactory(f"L{min(index1, index2)}{max(index1, index2)}")

        raise TypeError(f"datatype1 and datatype2 must be both code or both carrier."
                        f"datatype1 = {str(datatype1)}, datatype2 = {str(datatype2)}")

    @staticmethod
    def get_smooth_datatype(datatype):
        """
        mix data types to get the associated smooth pseudorange (SPR) or smooth iono free (IFSPR) datatype
        get_smooth_datatype(C1) => SPR1
        get_smooth_datatype(C12) => SPR12
        """
        if DataType.is_code(datatype) or DataType.is_carrier(datatype):
            # datatype = C1 -> return SPR1
            index = datatype.freq_number
            return DataTypeFactory(f"SPR{index}")

        if DataType.is_iono_free_code(datatype) or DataType.is_iono_free_carrier(datatype):
            # datatype = C12 -> return SPR12
            index = datatype.freq_number
            return DataTypeFactory(f"SPR{index}")

        raise TypeError(f"Unable to create smooth pseudorange from "
                        f"data types {str(datatype)}")


# Default Data Types

# Frequencies
f1 = DataType(data_type="f1", description="Frequency L1 (GPS)", freq_value=Constant.L1_FREQ, freq_number=1)
f2 = DataType(data_type="f2", description="Frequency L2 (GPS)", freq_value=Constant.L2_FREQ, freq_number=2)
f5 = DataType(data_type="f5", description="Frequency L5 (GPS)", freq_value=Constant.L2_FREQ, freq_number=5)

###################
# Raw Observables #
###################

# Code
C1 = DataType(data_type="C1", description="PseudoRange in Frequency L1", freq=f1, freq_number=1)
C2 = DataType(data_type="C2", description="PseudoRange in Frequency L2", freq=f2, freq_number=2)
C5 = DataType(data_type="C5", description="PseudoRange in Frequency L5", freq=f5, freq_number=5)

# Carriers
L1 = DataType(data_type="L1", description="CarrierPhase in Frequency L1", freq=f1, freq_number=1)
L2 = DataType(data_type="L2", description="CarrierPhase in Frequency L2", freq=f2, freq_number=2)
L5 = DataType(data_type="L5", description="CarrierPhase in Frequency L5", freq=f5, freq_number=5)

# Signals
S1 = DataType(data_type="S1", description="Signal Strength in Frequency L1", freq=f1, freq_number=1)
S2 = DataType(data_type="S2", description="Signal Strength in Frequency L2", freq=f2, freq_number=2)
S5 = DataType(data_type="S5", description="Signal Strength in Frequency L5", freq=f5, freq_number=5)

#########################
# Iono Free Observables #
#########################
C12 = DataType(data_type="C12", description="L1-L2 Iono-Free PseudoRange", freq_number=12)
C15 = DataType(data_type="C15", description="L1-L5 Iono-Free PseudoRange", freq_number=15)
C25 = DataType(data_type="C25", description="L2-L5 Iono-Free PseudoRange", freq_number=25)
L12 = DataType(data_type="L12", description="L1-L2 Iono-Free Carrier Phase", freq_number=12)
L15 = DataType(data_type="L15", description="L1-L5 Iono-Free Carrier Phase", freq_number=15)
L25 = DataType(data_type="L25", description="L2-L5 Iono-Free Carrier Phase", freq_number=25)

######################
# Smooth Observables #
######################
SPR1 = DataType(data_type="SPR1", description="Smooth PseudoRange in Frequency L1", freq=f1, freq_number=1)
SPR2 = DataType(data_type="SPR2", description="Smooth PseudoRange in Frequency L2", freq=f2, freq_number=2)
SPR5 = DataType(data_type="SPR5", description="Smooth PseudoRange in Frequency L5", freq=f5, freq_number=5)

################################
# Smooth Iono-Free Observables #
################################
SPR12 = DataType(data_type="SPR12", description="L1-L2 Iono-Free Smooth PseudoRange", freq_number=12)
SPR15 = DataType(data_type="SPR15", description="L1-L5 Iono-Free Smooth PseudoRange", freq_number=15)
SPR25 = DataType(data_type="SPR25", description="L2-L5 Iono-Free Smooth PseudoRange", freq_number=25)

# Unknown
UN = DataType(data_type="UN", description="Unknown inputs type")

# Containers
cAvailableCodes = [C1, C2, C5]
cAvailableSignals = [S1, S2, S5]
cAvailableCarriers = [L1, L2, L5]
cAvailableFrequencies = [f1, f2, f5]
cAvailableSmoothCodes = [SPR1, SPR2, SPR5]
cAvailableIonoFreeSmoothCodes = [SPR12, SPR15, SPR25]
cAvailableIonoFreeCodes = [C12, C15, C25]
cAvailableIonoFreeCarriers = [L12, L15, L25]


def DataTypeFactory(datatype: str):
    """
    Factory for inputs types. Receives a string representing the type of the datatype and returns the associated
    DataType instance

    Args:
        datatype (str): string with the short descriptor of the datatype to fetch (ex: C1, L5,...)
    Return:
         DataType : returns the corresponding DataType instance
    """
    for container in [cAvailableCodes, cAvailableSignals, cAvailableFrequencies, cAvailableCarriers,
                      cAvailableSmoothCodes, cAvailableIonoFreeCodes, cAvailableIonoFreeCarriers,
                      cAvailableIonoFreeSmoothCodes]:
        for _type in container:
            if _type.data_type == datatype:
                return _type

    return UN
