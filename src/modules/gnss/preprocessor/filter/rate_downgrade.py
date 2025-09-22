""" Rate Downgrade Filter Module """
from .filter import Filter


class RateDowngradeFilter(Filter):

    def __init__(self, rate_out: float, first_epoch, trace_path):
        """
        Constructor for the Rate Downgrade Filter.

        The `RateDowngradeFilter` is designed to reduce the sampling rate of a dataset by selectively removing data
        points according to a user-defined rate. This filter allows users to down-sample their data by specifying a
        new rate that is a multiple of the original sampling rate.

        Args:
            rate_out(float): desired output rate in seconds
            first_epoch(src.data_types.date.Epoch): epoch of first observation
            trace_path(str or None): path to the trace file
        """
        super().__init__()
        self.rate_out = rate_out
        self.first_epoch = first_epoch
        self.write_header(trace_path)

    def write_header(self, trace_path):
        if trace_path is not None:
            self.fd = open(trace_path + "/RateDowngradeFilter.txt", "w")
            self.fd.write("Epoch, Satellite, Observation, To Remove\n")

    def is_applicable(self, sat, epoch, observation, **kwargs):
        to_remove = (epoch - self.first_epoch).total_seconds() % self.rate_out != 0

        if self.fd is not None:
            self.fd.write(f"{epoch}, {epoch}, {observation.datatype}, {to_remove}\n")

        return to_remove
