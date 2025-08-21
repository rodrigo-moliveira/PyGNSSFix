""" Module with definition of Faults. """
import abc
from typing import Dict, Any


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
    def apply(self, state: Dict[str, Any], epoch: int) -> Dict[str, Any]:
        """
        Apply the fault to the system state or measurements.

        Args:
            state(dict): System state or measurement set at current epoch.
            epoch (int): Current epoch index.

        Returns:
            dict: Modified state/measurements.
        """
        pass

    def __str__(self):
        """Generic string representation for debugging."""
        return f"{self.__class__.__name__}({self.params})"


# ================================
# Cycle Slip Fault Example
# ================================
class CycleSlipFault(Fault):
    """ Injects a cycle slip into carrier phase measurements. """

    def apply(self, state: Dict[str, Any], epoch: int) -> Dict[str, Any]:
        sat_id = self.params["sat_id"]
        freq = self.params["freq"]
        ep_start = int(self.params["epoch_start"])
        ep_end = int(self.params["epoch_end"])
        val = float(self.params["val"])

        if ep_start <= epoch <= ep_end:
            if "CP" in state and sat_id in state["CP"]:
                if freq in state["CP"][sat_id]:
                    state["CP"][sat_id][freq] += val
        return state

    def __str__(self):
        return (f"CycleSlipFault(sat={self.params['sat_id']}, "
                f"freq={self.params['freq']}, "
                f"epoch=[{self.params['epoch_start']},{self.params['epoch_end']}], "
                f"val={self.params['val']})")


# ================================
# Example: Filter Reset Fault
# ================================
class FilterResetFault(Fault):
    """Triggers a Kalman filter reset at a specific epoch."""

    def apply(self, state: Dict[str, Any], epoch: int) -> Dict[str, Any]:
        trigger_epoch = int(self.params["epoch"])
        if epoch == trigger_epoch:
            state["KF_reset"] = True
        return state

    def __str__(self):
        return f"FilterResetFault(epoch={self.params['epoch']})"
