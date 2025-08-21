""" Module Manager class for fault injections """

from typing import List, Dict, Any

from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger
from src.fault.fault import Fault, CycleSlipFault, FilterResetFault


class FaultInjector:
    """Handles parsing a fault file and applying faults."""

    FAULT_CLASSES = {
        "CYCLE_SLIP": CycleSlipFault,
        "FILTER_RESET": FilterResetFault,
        # Future: "PR_BIAS": PRBiasFault, "OUTAGE": OutageFault, ...
    }

    def __init__(self, config_dict: dict):
        self.enabled = config_dict.get("fault_injector", "enabled")
        self.faults: List[Fault] = []

        if self.enabled:
            log = get_logger(IO_LOG)
            fault_path = config_dict.get("fault_injector", "fault_file")
            fault_file = WORKSPACE_PATH / fault_path
            log.info(f"Loading Fault Injector file {fault_file}...")
            self.load_faults(fault_file)

    def load_faults(self, fault_file: str):
        """Reads the fault definition file and creates fault objects."""
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

    def apply_faults(self, state: Dict[str, Any], epoch: int) -> Dict[str, Any]:
        """Applies all active faults for a given epoch."""
        if self.enabled:
            for fault in self.faults:
                state = fault.apply(state, epoch)
            return state

    def __str__(self):
        """Pretty print the list of loaded faults."""
        if not self.faults:
            return "FaultInjector(no faults loaded)"
        return f"FaultInjector (enabled={self.enabled}) with faults:\n  " + "\n  ".join(str(f) for f in self.faults)
