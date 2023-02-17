from PositioningSolver.src.data_types.basics.DataType import DataType


class Observation:
    """
    Class Observation
    Stores a single observation
    Attributes
        ----------
        datatype : DataType
            The DataType (observation type) of this observation
        value : float or int
            The numeric value of this observation
    """

    def __init__(self, datatype: DataType, value: float):
        if isinstance(datatype, DataType):
            self.datatype = datatype
        else:
            raise TypeError(f'First argument of Observation constructor should be of type DataType, {type(datatype)}'
                            f' was provided instead.')
        if isinstance(value, float) or isinstance(value, int):
            self.value = float(value)
        else:
            raise TypeError(f'Second argument of Observation constructor should be of type float, {type(value)}'
                            f' was provided instead')

    def __eq__(self, other):
        return self.datatype == other.datatype

    def __repr__(self):
        return "Observation({} : {})".format(self.value, self.datatype)

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
