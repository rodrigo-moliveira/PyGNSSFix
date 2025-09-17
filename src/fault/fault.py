""" Module with definition of Faults. """
import abc
from typing import Dict, Any

from src.data_types.date import Epoch


# TODO: class cleanup

class Fault(abc.ABC):
    """Abstract base class for all faults."""

    def __init__(self, params: Dict[str, Any]):
        """
        Initialize fault with arbitrary parameters.

        Args:
            params (dict): Dictionary of fault-specific parameters, parsed from the fault file.
        """
        self.params = params

    @abc.abstractmethod
    def apply(self, state_in):
        return state_in

    @abc.abstractmethod
    def check_fault(self, **params):
        # TODO: docstring
        pass

    def __str__(self):
        """Generic string representation for debugging."""
        return f"{self.__class__.__name__}({self.params})"


# ==============================
# Measurement Bias Fault Example
# ==============================
class MeasurementBias(Fault):
    """ Injects a measurement bias into the provided GNSS observation. """

    def __init__(self, params: Dict[str, Any]):
        super().__init__(params)

        # additional initializations
        self.epoch_start = Epoch.strptime(self.params["epoch_start"], scale=self.params["time_scale"])
        self.epoch_end = Epoch.strptime(self.params["epoch_end"], scale=self.params["time_scale"])
        self.bias = float(self.params["bias"])
        self.sat = self.params["sat_id"]
        self.rinex_obs = self.params["rinex_obs"]

    def apply(self, state_in):
        return state_in + self.bias

    def check_fault(self, **params):
        epoch = params["epoch"]
        sat = params["sat"]
        obs = params["obs"]

        # epoch check
        if self.epoch_start <= epoch <= self.epoch_end:

            # satellite check
            if str(sat) == self.sat:

                # observation type check
                if self.rinex_obs == obs:

                    return True
        return False

    def __str__(self):
        return (f"MeasurementBias(sat={self.sat}, "
                f"obs={self.rinex_obs}, "
                f"epoch=[{str(self.epoch_start)},{str(self.epoch_end)}], "
                f"bias={self.bias})")


# ================================
# Example: Filter Reset Fault
# ================================
class FilterResetFault(Fault):
    """Triggers a Kalman filter reset at a specific epoch."""

    def __init__(self, params: Dict[str, Any]):
        super().__init__(params)

        # additional initializations
        self.params["epoch"] = Epoch.strptime(self.params["epoch"], scale=self.params["time_scale"])

    def apply(self, state_in):
        return state_in

    def check_fault(self, **params):
        return False

    def __str__(self):
        return f"FilterResetFault(epoch={self.params['epoch']})"
