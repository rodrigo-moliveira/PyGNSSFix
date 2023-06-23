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
    APRIORI = 1
    IONO_FREE = 2

    @classmethod
    def show_options(cls):
        return f"[0 - DISABLED, 1 - A PRIORI METHOD, 2 - IONO FREE COMBINATIONS]"


class EnumTropo(Enum):
    DISABLED = 0
    SAASTAMOINEM = 1

    @classmethod
    def show_options(cls):
        return f"[0 - DISABLED, 1 - SAASTAMOINEM MODEL]"


class EnumTransmissionTime(Enum):
    GEOMETRIC = 0
    PSEUDORANGE = 1

    @classmethod
    def show_options(cls):
        return f"[0 - GEOMETRIC, 1 - PSEUDORANGE]"
