from enum import Enum


# Enum classes already perform input validation
# >>> EnumSolver(1) is okay
# >>> EnumSolver(10) will return ValueError

class EnumSolver(Enum):
    LS = 0
    WLS = 1

    @classmethod
    def show_options(cls):
        return f"[0 - LEAST SQUARES, 1 - WEIGHTED LEAST SQUARES]"


class EnumOnOff(Enum):
    DISABLED = 0
    ENABLED = 1

    @classmethod
    def show_options(cls):
        return f"[0 - DISABLED, 1 - ENABLED]"


class EnumIono(Enum):
    DISABLED = 0
    ENABLED = 1
    APRIORI = 2

    @classmethod
    def show_options(cls):
        return f"[0 - DISABLED, 1 - ENABLED, 2 - A PRIORI METHOD]"


class EnumTransmissionTime(Enum):
    GEOMETRIC = 0
    PSEUDORANGE = 1

    @classmethod
    def show_options(cls):
        return f"[0 - GEOMETRIC, 1 - PSEUDORANGE]"
