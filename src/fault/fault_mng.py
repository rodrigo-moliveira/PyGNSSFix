""" Module Manager class for fault injections """

from typing import List, Dict, Any

from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger
from src.fault.fault import Fault, MeasurementBias, FilterResetFault


# TODO: class cleanup

class FaultInjector:
    """Handles parsing a fault file and applying faults."""

    # Available faults
    FAULT_CLASSES = {
        "MEAS_BIAS": MeasurementBias,
        "FILTER_RESET": FilterResetFault,
        # Future: "PR_BIAS": PRBiasFault, "OUTAGE": OutageFault, ...
    }

    def __init__(self, config_dict: dict):
        """
        Constructor of the Fault Manager.

        Args:
            config_dict (dict): dict instance with the user configurations

        Raises:
            ValueError: an exception is raised if the provided fault type is not correct
        """
        self.enabled = config_dict.get("fault_injector", "enabled")
        self.faults: List[Fault] = []

        if self.enabled:
            log = get_logger(IO_LOG)
            fault_path = config_dict.get("fault_injector", "fault_file")
            fault_file = WORKSPACE_PATH / fault_path
            log.info(f"Loading Fault Injector file {fault_file}...")
            self.load_faults(fault_file)

    def load_faults(self, fault_file: str):
        """ Reads the fault definition file and creates fault objects. """
        if self.enabled:
            with open(fault_file, "r") as f:
                for line in f:
                    if line.strip().startswith("#") or not line.strip():
                        continue
                    parts = [p.strip() for p in line.strip().split(",")]
                    fault_id = parts[0]

                    fault_cls = self.FAULT_CLASSES.get(fault_id)
                    if not fault_cls:
                        raise ValueError(f"Unknown fault type: {fault_id}")

                    # Convert remaining key=val pairs into dict
                    params = {}
                    for token in parts[1:]:
                        if "=" in token:
                            k, v = token.split("=")
                            params[k.strip()] = v.strip()

                    fault = fault_cls(params)
                    self.faults.append(fault)

    def check_faults(self, fault_type, state_in, **params):
        """
        Checks if there are any faults of `fault_type` to be applied at the provided `epoch`

        Args:
            epoch():
            fault_type():

        Returns:
            bool: True if there faults to be injected for this type and epoch, and false otherwise

        Raises:
            ValueError: an exception is raised if the provided fault type is not correct
        """
        state_out = state_in
        fault_cls = self.FAULT_CLASSES.get(fault_type)
        if not fault_cls:
            raise ValueError(f"Unknown fault type: {fault_type}")

        for fault in self.faults:
            if isinstance(fault, fault_cls):
                if fault.check_fault(**params):
                    state_out = fault.apply(state_out)

        return state_out

    def __str__(self):
        """Pretty print the list of loaded faults."""
        if not self.faults:
            return "FaultInjector(no faults loaded)"
        return f"FaultInjector (enabled={self.enabled}) with faults:\n  " + "\n  ".join(str(f) for f in self.faults)
