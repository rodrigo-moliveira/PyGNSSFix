""" SNR Filter Module """
from src.data_types.gnss.data_type import DataType
from . import Filter


class SignalCheckFilter(Filter):

    def __init__(self, observation_data, snr_threshold, trace_path):
        """
        Constructor of the Signal Check Filter

        The `SignalCheckFilter` is designed to clean the quality of a dataset by removing data points that correspond
        to weak signals, based on the Carrier-to-Noise ratio (C/N0). This filter removes all observations whose C/N0
        (signal observable) is below the user-defined threshold.

        Parameters:
            observation_data(src.data_mng.gnss.ObservationData): input observation data set (required for the
                validation logic - read only)
            snr_threshold(float): threshold for the C/N0 noise
            trace_path(str or None): path to the trace file
        """
        super().__init__()
        self.observation_data = observation_data
        self.snr_threshold = snr_threshold
        self.write_header(trace_path)

    def write_header(self, trace_path):
        if trace_path is not None:
            self.fd = open(trace_path + "/SNRFilter.txt", "w")
            self.fd.write("Epoch, Satellite, Freq Band, SNR Observation, SNR Threshold, To Remove\n")

    def apply(self, sat, epoch, observation, v_removable):
        # only apply this filter to Signal Observables
        if not DataType.is_signal(observation.datatype):
            return

        to_remove = False
        if observation.value < self.snr_threshold:
            # remove all observables associated with this frequency
            to_remove = True
            v_observables = self.observation_data.get_observables_at_epoch(epoch, sat)

            # iterate over all observables and flag those with the same frequency as 'gnss_models'
            for this_obs in v_observables:
                if this_obs.datatype.freq == observation.datatype.freq:
                    if this_obs not in v_removable:
                        v_removable.append(this_obs)

        if self.fd is not None:
            self.fd.write(f"{str(epoch)}, {sat}, {observation.datatype.freq}, {observation.value}, {self.snr_threshold}, "
                          f"{to_remove}\n")
