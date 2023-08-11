from enum import Enum

from src.errors import EnumError


class EnumSolver(Enum):
    LS = 0
    WLS = 1

    @classmethod
    def show_options(cls):
        return f"[0 - LEAST SQUARES, 1 - WEIGHTED LEAST SQUARES]"


class EnumPositioningMode(Enum):
    SPS = 0
    SPS_IF = 1

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "sps":
            return EnumPositioningMode.SPS
        elif model_str.lower() == "sps_if":
            return EnumPositioningMode.SPS_IF
        else:
            raise EnumError(f"Unsupported model {model_str}. Available options are 'SPS', 'SPS_IF'")

    @classmethod
    def show_options(cls):
        return f"[0 - SPS, 1 - SPS (Iono-Free) ]"


class EnumModel(Enum):
    SINGLE_FREQ = 0
    DUAL_FREQ = 1

    @classmethod
    def show_options(cls):
        return f"[0 - SINGLE FREQUENCY MODEL, 1 - DUAL FREQUENCY MODEL]"


class EnumOnOff(Enum):
    DISABLED = 0
    ENABLED = 1

    @classmethod
    def show_options(cls):
        return f"[0 - DISABLED, 1 - ENABLED]"


class EnumIono(Enum):
    DISABLED = 0
    KLOBUCHAR = 1
    NEQUICK = 2

    @classmethod
    def show_options(cls):
        return f"[0 - DISABLED, 1 - Klobuchar (for GPS), 2 - NeQuick (for GAL)]"


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
