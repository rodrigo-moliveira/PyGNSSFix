""" Type Consistency Module """
from . import Filter


class TypeConsistencyFilter(Filter):

    def __init__(self, types, trace_path):
        """
        Constructor for the Type Consistency Filter
        The `TypeConsistencyFilter` cleans the Observation dataset by removing all epoch data for a given satellite
        that does not have all required observation types. The required observation types are provided as input to the
        filter

        Parameters:
            types(dict): Dictionary with constellation string as keys and list of all mandatory observation types as
                values
            trace_path(str or None): path to the trace file
        """
        super().__init__()
        self.mandatory = types
        self.write_header(trace_path)

    def write_header(self, trace_path):
        if trace_path is not None:
            self.fd = open(trace_path + "/TypeConsistentFilter.txt", "w")
            self.fd.write("Epoch, Satellite, Observable, To Remove\n")

    def is_applicable(self, sat, epoch, observation, **kwargs):
        # check if this datatype is part of the required observations
        to_remove = observation.datatype not in self.mandatory[sat.sat_system]

        if not to_remove:
            # check if all mandatory observations are available for this epoch
            for obs_required in self.mandatory[sat.sat_system]:
                found = False
                for obs in kwargs["obs_list"]:
                    if obs.datatype == obs_required:
                        found = True
                        break
                if not found:
                    to_remove = True  # delete this observable
                    break

        if self.fd is not None:
            self.fd.write(f"{epoch}, {sat}, {observation.datatype}, {to_remove}\n")
        return to_remove
