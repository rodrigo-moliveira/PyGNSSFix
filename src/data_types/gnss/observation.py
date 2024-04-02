from .data_type import DataType


class Observation:
    """
    Class Observation to store GNSS observables
    The Observation objects contain a :py:class:`DataType` and the corresponding value

    Attributes:
        datatype (DataType): The DataType (gnss_models type) of this GNSS observation
        value (float): The numeric value for this gnss_models
    """

    def __init__(self, datatype: DataType, value: float):
        """
        Constructor of `Observation` objects

        Args:
            datatype (DataType): The DataType (gnss_models type) of this GNSS observation
            value (float): The numeric value for this gnss_models
        """
        if isinstance(datatype, DataType):
            self.datatype = datatype
        else:
            raise AttributeError(f'First argument of Observation constructor should be of type DataType, '
                                 f'{type(datatype)} was provided instead.')
        if isinstance(value, float) or isinstance(value, int):
            self.value = float(value)
        else:
            raise AttributeError(f'Second argument of Observation constructor should be of type float, {type(value)}'
                                 f' was provided instead')

    def __eq__(self, other):
        if not isinstance(other, Observation):
            raise ValueError(f"Cannot compare DataType object")
        return self.datatype == other.datatype and self.value == other.value

    def __repr__(self):
        return "Observation({} : {})".format(self.value, repr(self.datatype))

    # Algebraic methods to allow an Observation object to be treated as a float
    def __add__(self, other):
        # operator self + other
        if isinstance(other, float) or isinstance(other, int):
            return Observation(self.datatype, self.value + other)
        elif isinstance(other, Observation):
            return Observation(self.datatype, self.value + other.value)
        else:
            raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

    def __iadd__(self, other):
        # operator self += other
        if isinstance(other, float) or isinstance(other, int):
            self.value += other
        elif isinstance(other, Observation):
            self.value += other.value
        else:
            raise TypeError(f"unsupported operand type(s) for +=: '{type(self)}' and '{type(other)}'")
        return self

    def __sub__(self, other):
        # operator self - other
        if isinstance(other, float) or isinstance(other, int):
            return Observation(self.datatype, self.value - other)
        elif isinstance(other, Observation):
            return Observation(self.datatype, self.value - other.value)
        else:
            raise TypeError(f"unsupported operand type(s) for -: '{type(self)}' and '{type(other)}'")

    def __isub__(self, other):
        # operator self -= other
        if isinstance(other, float) or isinstance(other, int):
            self.value -= other
        elif isinstance(other, Observation):
            self.value -= other.value
        else:
            raise TypeError(f"unsupported operand type(s) for -=: '{type(self)}' and '{type(other)}'")
        return self

    def __lt__(self, other):
        # operator self < other
        if isinstance(other, float) or isinstance(other, int):
            return self.value < other
        elif isinstance(other, Observation):
            return self.value < other.value
        else:
            raise TypeError(f"unsupported operand type(s) for <: '{type(self)}' and '{type(other)}'")

    def __le__(self, other):
        # operator self <= other
        if isinstance(other, float) or isinstance(other, int):
            return self.value <= other
        elif isinstance(other, Observation):
            return self.value <= other.value
        else:
            raise TypeError(f"unsupported operand type(s) for <=: '{type(self)}' and '{type(other)}'")

    def __gt__(self, other):
        # operator self > other
        if isinstance(other, float) or isinstance(other, int):
            return self.value > other
        elif isinstance(other, Observation):
            return self.value > other.value
        else:
            raise TypeError(f"unsupported operand type(s) for >: '{type(self)}' and '{type(other)}'")

    def __ge__(self, other):
        # operator self >= other
        if isinstance(other, float) or isinstance(other, int):
            return self.value >= other
        elif isinstance(other, Observation):
            return self.value >= other.value
        else:
            raise TypeError(f"unsupported operand type(s) for >=: '{type(self)}' and '{type(other)}'")

    def __float__(self):
        return self.value

    def copy(self):
        """
        Return:
            Observation: returns a deep copy of this object
        """
        return Observation(self.datatype, float(self.value))
