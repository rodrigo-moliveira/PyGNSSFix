""" Definition of some useful python Enumeration objects """

from enum import Enum
from src.errors import EnumError


class EnumSolver(Enum):
    """ Enumeration for the type of solver (Least-Squares or Weighted Least-Squares) """
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


class EnumAlgorithmPNT(Enum):
    """ Enumeration for the GNSS PNT algorithm (SPS, PR_PPP) """
    SPS = 0
    PR_PPP = 1
    CP_PPP = 2

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "sps":
            return EnumAlgorithmPNT.SPS
        elif model_str.lower() == "pr-ppp":
            return EnumAlgorithmPNT.PR_PPP
        elif model_str.lower() == "cp-ppp":
            return EnumAlgorithmPNT.CP_PPP
        else:
            raise EnumError(f"Unsupported algorithm {model_str}. Available options are 'SPS', 'PR-PPP', 'CP-PPP'.")

    @classmethod
    def show_options(cls):
        return f"[SPS, PR-PPP ]"


class EnumObservationModel(Enum):
    """
    Enumeration for the selected GNSS Observation model. Available models are:
        * uncombined model (single or dual frequency)
        * combined (iono-free combination)
    """
    UNCOMBINED = 0
    COMBINED = 1

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "uncombined":
            return EnumObservationModel.UNCOMBINED
        elif model_str.lower() == "combined":
            return EnumObservationModel.COMBINED
        else:
            raise EnumError(f"Unsupported observation model {model_str}. Available options are 'uncombined', "
                            f"'combined'")

    @classmethod
    def show_options(cls):
        return f"[0 - Uncombined Model (Single or Dual Frequency), 1 - Combined (Iono-Free Combination)]"


class EnumOnOff(Enum):
    """ Enumeration for the activation or deactivation of models (Disabled or Enabled) """
    DISABLED = 0
    ENABLED = 1

    @classmethod
    def show_options(cls):
        return f"[0 - DISABLED, 1 - ENABLED]"


class EnumFrequencyModel(Enum):
    """ Enumeration for the Frequency Model to be implemented in the estimation process (single or dual frequency) """
    SINGLE_FREQ = 0
    DUAL_FREQ = 1

    @classmethod
    def show_options(cls):
        return f"[ SINGLE_FREQ, DUAL_FREQ]"

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "single":
            return EnumFrequencyModel.SINGLE_FREQ
        elif model_str.lower() == "dual":
            return EnumFrequencyModel.DUAL_FREQ
        else:
            raise EnumError(f"Unsupported Frequency Model {model_str}. Available options are {cls.show_options()}")


class EnumIonoModel(Enum):
    """ Enumeration for the Ionosphere A-priori Model (Disabled, Klobuchar, NTCM-G or IONEX) """
    DISABLED = 0
    KLOBUCHAR = 1
    NTCMG = 2
    IONEX = 3

    @classmethod
    def show_options(cls):
        return f"[ NONE, Klobuchar, NTCM-G, IONEX]"

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "none":
            return EnumIonoModel.DISABLED
        elif model_str.lower() == "klobuchar":
            return EnumIonoModel.KLOBUCHAR
        elif model_str.lower() == "ntcm-g" or model_str.lower() == "ntcmg":
            return EnumIonoModel.NTCMG
        elif model_str.lower() == "ionex":
            return EnumIonoModel.IONEX
        else:
            raise EnumError(f"Unsupported Ionospheric Model {model_str}. Available options are {cls.show_options()}")


class EnumTropoModel(Enum):
    """ Enumeration for the Troposphere A-priori Model (Disabled, Saastamoinen, GPT3) """
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


class EnumTropoMap(Enum):
    """ Enumeration for the Troposphere Map function (Saastamoinen, GMF, VMF1, VMF3) """
    SAASTAMOINEM = 0
    GMF = 1
    VMF1 = 2
    VMF3 = 3

    @classmethod
    def init_model(cls, model_str: str):
        if model_str.lower() == "saastamoinen":
            return EnumTropoMap.SAASTAMOINEM
        elif model_str.lower() == "gmf":
            return EnumTropoMap.GMF
        elif model_str.lower() == "vmf1":
            return EnumTropoMap.VMF1
        elif model_str.lower() == "vmf3":
            return EnumTropoMap.VMF3
        else:
            raise EnumError(f"Unsupported Tropospheric Map Function {model_str}. "
                            f"Available options are {cls.show_options()}")

    @classmethod
    def show_options(cls):
        return f"[ Saastamoinen, GMF, VMF1, VMF3]"


class EnumTransmissionTime(Enum):
    """ Enumeration for the Computation of Transmission Time model (Geometric, Pseudorange, NAPEOS) """
    GEOMETRIC = 0
    PSEUDORANGE = 1
    NAPEOS = 2

    @classmethod
    def show_options(cls):
        return f"[0 - GEOMETRIC, 1 - PSEUDORANGE, 2 - NAPEOS]"


class EnumSatelliteBias(Enum):
    """
    Enumeration for the selected Satellite Bias type. The available types are:
        * broadcast code biases (BGD for GAL and TGD for GPS)
        * precise DCBs (differential code bias)
        * precise OSBs (observation specific bias)
    """
    BROADCAST = 0
    DCB = 1
    OSB = 2

    @classmethod
    def show_options(cls):
        return f"[0 - BGD/TGD, 1 - DCB, 2 - OSB]"


class EnumPCVModel(Enum):
    """
    Enumeration for the Phase Center Variation (PCV) Model:
        * 0: non-azimuth dependent
        * 1: azimuth dependent
    """
    NON_AZI_DEPENDENT = 0
    AZI_DEPENDENT = 1

    @classmethod
    def show_options(cls):
        return f"[0 - non-azimuth dependent, 1 - azimuth dependent]"
