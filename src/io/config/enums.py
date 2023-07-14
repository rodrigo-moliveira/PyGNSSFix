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

    @classmethod
    def show_options(cls):
        return f"[0 - DISABLED, 1 - A PRIORI METHOD]"


class EnumCombined(Enum):
    UNCOMBINED_MODEL = 0
    COMBINED_MODEL = 1

    @classmethod
    def show_options(cls):
        return f"[0 - UNCOMBINED OBSERVATIONS MODEL, 1 - COMBINED OBSERVATIONS MODEL]"

    def __repr__(self):
        if self == EnumCombined.COMBINED_MODEL:
            return "Combined Model (Iono-Free Observation)"
        if self == EnumCombined.UNCOMBINED_MODEL:
            return "Uncombined Model"


class EnumModel(Enum):
    SINGLE_FREQ = 0
    DUAL_FREQ = 1

    @classmethod
    def show_options(cls):
        return f"[0 - SINGLE FREQUENCY MODEL, 1 - DUAL FREQUENCY MODEL]"


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
