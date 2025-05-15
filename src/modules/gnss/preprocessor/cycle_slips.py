from math import sqrt

from src.constants import SPEED_OF_LIGHT
from src.data_mng.gnss import ObservationData
from src.data_mng import TimeSeries
from src.errors import TimeSeriesError


class CycleSlipDetector:
    def __init__(self, mw_data: ObservationData):
        self.mw_data = mw_data

        self.detect_cycle_slip_mw()

    def detect_cycle_slip_mw(self):
        """
        Detect cycle slips in the Melbourne-Wubbena combination of carrier phase measurements using a mean and standard
        deviation based method.

        Reference:
            https://gssc.esa.int/navipedia/index.php?title=Detector_based_in_code_and_carrier_phase_data:_The_Melbourne-W%C3%BCbbena_combination

        """
        # Implement the MW cycle slip detection algorithm here
        vEpochs = self.mw_data.get_epochs()
        sat_list = self.mw_data.get_satellites()
        k_factor = 4.0
        max_gap = 5.0  # in seconds
        window_size = 100

        mean_mw = dict()
        cov_mw = dict()
        cycle_slips = dict()
        for sat in sat_list:
            mean_mw[sat] = TimeSeries()
            cov_mw[sat] = TimeSeries()
            cycle_slips[sat] = TimeSeries()

        for sat in sat_list:
            mean = 0.0
            cov = 0.0
            iWindow = 0

            for iEpoch, epoch in enumerate(vEpochs):
                k = iWindow + 1
                try:
                    mw_obs = self.mw_data.get_observables_at_epoch(epoch, sat)
                except TimeSeriesError:
                    # reset
                    iWindow = 0
                    mean = 0
                    cov = 0
                    continue

                obs = mw_obs[0].value
                # checks
                # 1 - time gap check
                # ...

                # 2 - comparison
                if iWindow > 0:
                    prev_epoch = vEpochs[iEpoch - 1]
                    if mean_mw[sat].has_epoch(prev_epoch):
                        prev_mean = mean_mw[sat].get_data_for_epoch(prev_epoch)
                        prev_cov = cov_mw[sat].get_data_for_epoch(prev_epoch)
                        if obs - prev_mean > k_factor * sqrt(prev_cov):
                            print(f"Cycle slip detected at epoch {epoch} for satellite {sat}.")
                            cycle_slips[sat].set_data(epoch, True)

                # Update mean and covariance
                mean = (k - 1) / k * mean + obs / k
                if iWindow == 0:
                    # initialize covariance for first epoch in window
                    cov = SPEED_OF_LIGHT / mw_obs[0].datatype.freq_value / 2
                else:
                    cov = (k - 1) / k * cov + (obs - mean) ** 2 / k
                mean_mw[sat].set_data(epoch, mean)
                cov_mw[sat].set_data(epoch, cov)

                iWindow += 1

        #exit()
        return cycle_slips
