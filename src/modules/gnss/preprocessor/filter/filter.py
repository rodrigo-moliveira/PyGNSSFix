""" Module with the base Filter Class """

from src.data_types.gnss.observation import Observation
from src.data_types.gnss.satellite import Satellite
from src.data_types.date import Epoch


class Filter:
    """
    Base class for the Filters.

    Filter classes inherited from this must implement the following methods:
        * ``is_applicable``
        * ``write_header``
    """
    def __init__(self):
        self.fd = None

    def write_header(self, trace_path):
        """ Open the trace file and write the header content """
        pass

    def close_file(self):
        """ Close the trace file. """
        if self.fd is not None:
            self.fd.close()

    def is_applicable(self, sat: Satellite, epoch: Epoch, observation: Observation, **kwargs):
        """
        Returns:
            bool: return True to remove this `observation`, return False to keep it
        """
        return True

    def apply(self, sat: Satellite, epoch: Epoch, observation: Observation, v_removable: list):
        """
        If the observation triggers the filter, it is appended to the `v_removable` list (to be removed).
        """
        v_removable.append(observation)
