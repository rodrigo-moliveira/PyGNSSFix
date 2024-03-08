"""Definition of some useful python Enumeration objects"""

from enum import Enum
from src.errors import EnumError


class EnumSolver(Enum):
    LS = 0
    WLS = 1

    @classmethod
    def show_options(cls):
        return f"[ LS, WSL ]"

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "wsl":
            return EnumSolver.WLS
        elif model_str.lower() == "ls":
            return EnumSolver.LS
        else:
            raise EnumError(f"Unsupported Solver Model {model_str}. Available options are {cls.show_options()}")


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
        return f"[SPS, SPS_IF ]"


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
    NTCMG = 2

    @classmethod
    def show_options(cls):
        return f"[ NONE, Klobuchar, NTCM-G]"

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "none":
            return EnumIono.DISABLED
        elif model_str.lower() == "klobuchar":
            return EnumIono.KLOBUCHAR
        elif model_str.lower() == "ntcm-g" or model_str.lower() == "ntcmg":
            return EnumIono.NTCMG
        else:
            raise EnumError(f"Unsupported Ionospheric Model {model_str}. Available options are {cls.show_options()}")


class EnumTropoModel(Enum):
    DISABLED = 0
    SAASTAMOINEM = 1
    GPT3 = 2

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "none":
            return EnumTropoModel.DISABLED
        elif model_str.lower() == "saastamoinen":
            return EnumTropoModel.SAASTAMOINEM
        elif model_str.lower() == "gpt3":
            return EnumTropoModel.GPT3
        else:
            raise EnumError(f"Unsupported Tropospheric Model {model_str}. Available options are {cls.show_options()}")

    @classmethod
    def show_options(cls):
        return f"[ NONE, Saastamoinen, GPT3]"


class EnumTropoMask(Enum):
    SAASTAMOINEM = 0
    GMF = 1
    VMF1 = 2
    VMF3 = 3

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "saastamoinen":
            return EnumTropoMask.SAASTAMOINEM
        elif model_str.lower() == "gmf":
            return EnumTropoMask.GMF
        elif model_str.lower() == "vmf1":
            return EnumTropoMask.VMF1
        elif model_str.lower() == "vmf3":
            return EnumTropoMask.VMF3
        else:
            raise EnumError(f"Unsupported Tropospheric Mask Function {model_str}. "
                            f"Available options are {cls.show_options()}")

    @classmethod
    def show_options(cls):
        return f"[ Saastamoinen, GMF, VMF1, VMF3]"


class EnumTransmissionTime(Enum):
    GEOMETRIC = 0
    PSEUDORANGE = 1

    @classmethod
    def show_options(cls):
        return f"[0 - GEOMETRIC, 1 - PSEUDORANGE]"
