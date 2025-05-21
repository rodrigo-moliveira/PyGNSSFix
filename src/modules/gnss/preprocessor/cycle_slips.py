""" Module for detecting cycle slips in GNSS data """
from math import sqrt
import numpy as np

from src.common_log import PREPROCESSOR_LOG, get_logger
from src.constants import SPEED_OF_LIGHT
from src.data_mng.gnss import ObservationData
from src.data_mng import TimeSeries
from src.errors import TimeSeriesError
from src.io.config import config_dict


class CycleSlipDetector:
    """
    Class for detecting cycle slips in GNSS data using Melbourne-Wubbena and Geometry-Free methods.

    Attributes:
        mw_data (ObservationData): Melbourne-Wubbena observation data.
        gf_data (ObservationData): Geometry-Free observation data.
        log (logging.Logger): Logger for the class.
        write_trace (bool): Flag to write trace data to file.
        trace_path (str): Path to write trace data.
        mw_k_factor (float): K-factor for Melbourne-Wubbena method.
        max_gap (int): Maximum gap in number of measurements.
        gf_threshold (float): Threshold for Geometry-Free method.
        cycle_slips (dict): Dictionary to store cycle slips for each satellite (output)

    The detected cycle slips are stored in the `cycle_slips` attribute, which is a dictionary where
    the keys are satellite identifiers and the values are TimeSeries objects containing the cycle slip data.
    """
    def __init__(self, mw_data: ObservationData, gf_data, write_trace, trace_path):
        self.mw_data = mw_data
        self.gf_data = gf_data
        self.log = get_logger(PREPROCESSOR_LOG)
        self.write_trace = write_trace
        self.trace_path = trace_path

        # configurations
        self.mw_k_factor = config_dict.get("preprocessor", "cycle_slips", "mw_k_factor")
        self.max_gap = config_dict.get("preprocessor", "cycle_slips", "max_gap")  # in number of measurements
        self.gf_threshold = config_dict.get("preprocessor", "cycle_slips", "gf_threshold")

        self.cycle_slips = None

    def compute(self):
        """
        Main method to compute cycle slips using Melbourne-Wubbena and Geometry-Free methods.
        """
        cycle_slip_mw, cycle_slip_gf = dict(), dict()
        if not self.mw_data.is_empty():
            f = None
            if self.write_trace:
                f = open(self.trace_path + "/MW_CycleSlipDetection.txt", "w")
            try:
                self.detect_cycle_slip_mw(cycle_slip_mw, f)
            except Exception as e:
                self.log.error(f"Error in Melbourne Wubbena cycle slip detection: {e}")
            if self.write_trace:
                f.close()
        if not self.gf_data.is_empty():
            f = None
            if self.write_trace:
                f = open(self.trace_path + "/GF_CycleSlipDetection.txt", "w")
            try:
                self.detect_cycle_slip_gf(cycle_slip_gf, f)
            except Exception as e:
                self.log.error(f"Error in Geometry Free cycle slip detection: {e}")
            if self.write_trace:
                f.close()

        # merged cycle slips
        self._merge_results(cycle_slip_mw, cycle_slip_gf)

    def detect_cycle_slip_mw(self, cycle_slips, f_handle=None):
        """
        Detect cycle slips in the Melbourne-Wubbena combination of carrier phase measurements using a mean and standard
        deviation based method.

        Reference:
            https://gssc.esa.int/navipedia/index.php?title=Detector_based_in_code_and_carrier_phase_data:_The_Melbourne-W%C3%BCbbena_combination

        """
        vEpochs = self.mw_data.get_epochs()
        vSatellites = self.mw_data.get_satellites()
        sampling_rate = self.mw_data.get_rate()

        mean_mw = dict()
        cov_mw = dict()
        for sat in vSatellites:
            mean_mw[sat] = TimeSeries()
            cov_mw[sat] = TimeSeries()
            cycle_slips[sat] = TimeSeries()

        for sat in vSatellites:
            mean = 0.0
            cov = 0.0
            iWindow = 0
            prev_epoch = None

            for iEpoch, epoch in enumerate(vEpochs):
                n_elements = iWindow + 1

                # fetch observation for this epoch and satellite
                try:
                    mw_obs = self.mw_data.get_observables_at_epoch(epoch, sat)
                    obs = mw_obs[0].value
                except TimeSeriesError:
                    continue

                # comparisons
                if iWindow > 0:
                    if prev_epoch is not None:
                        detect = False
                        if (epoch - prev_epoch).total_seconds() > self.max_gap * sampling_rate:
                            if f_handle:
                                f_handle.write(f"Cycle slip detected at epoch {epoch} for satellite {sat}: maximum gap "
                                               f"{self.max_gap * sampling_rate} seconds exceeded.\n")
                            detect = True

                        prev_mean = mean_mw[sat].get_data_for_epoch(prev_epoch)
                        prev_cov = cov_mw[sat].get_data_for_epoch(prev_epoch)
                        if obs - prev_mean > self.mw_k_factor * sqrt(prev_cov):
                            if f_handle:
                                f_handle.write(f"Cycle slip detected at epoch {epoch} for satellite {sat}: "
                                               f"difference obs - prev_mean = {obs - prev_mean} exceeds threshold "
                                               f"{self.mw_k_factor * sqrt(prev_cov)}\n")
                            detect = True

                        if f_handle:
                            f_handle.write(f"{epoch}: {sat}: MW observable = {obs}. Previous mean: {prev_mean}. " 
                                            f"Previous covariance: {prev_cov}. Difference: {obs - prev_mean}. "
                                            f"Threshold: {self.mw_k_factor * sqrt(prev_cov)}\n")

                        if detect:
                            # cycle slip detected: reset algorithm
                            cycle_slips[sat].set_data(epoch, True)
                            mean = 0.0
                            cov = 0.0
                            iWindow = 0
                            prev_epoch = None
                            continue

                # Update mean and covariance
                mean = (n_elements - 1) / n_elements * mean + obs / n_elements
                if iWindow == 0:
                    cov = SPEED_OF_LIGHT / mw_obs[0].datatype.freq_value / 2
                else:
                    cov = (n_elements - 1) / n_elements * cov + (obs - mean) ** 2 / n_elements
                mean_mw[sat].set_data(epoch, mean)
                cov_mw[sat].set_data(epoch, cov)
                prev_epoch = epoch
                iWindow += 1

    def detect_cycle_slip_gf(self, cycle_slips, f_handle=None):
        """
        Detect cycle slips in the geometry-free combination of carrier phase measurements using a polynomial
        fitting method.

        Reference:
            https://gssc.esa.int/navipedia/index.php/Detector_based_in_carrier_phase_data:_The_geometry-free_combination

        """
        vEpochs = self.gf_data.get_epochs()
        vSatellites = self.gf_data.get_satellites()
        sampling_rate = self.gf_data.get_rate()
        threshold = self.gf_threshold
        window_size = 10

        for sat in vSatellites:
            cycle_slips[sat] = TimeSeries()

        for sat in vSatellites:

            if threshold < 0:
                threshold = self._get_gf_threshold(sat, sampling_rate)

            iLast_slip = -1
            prev_epoch = None
            for iEpoch, epoch in enumerate(vEpochs):

                # fetch GF observable for this epoch
                try:
                    obs_list = self.gf_data.get_observables_at_epoch(epoch, sat)
                    for obs in obs_list:
                        if "CP" in obs.datatype.data_type:
                            gf_obs = obs
                            break
                    else:
                        continue
                except TimeSeriesError:
                    continue

                # declare cycle slip when data gap greater than tolerance
                if prev_epoch is not None and (epoch - prev_epoch).total_seconds() > self.max_gap * sampling_rate:
                    if f_handle:
                        f_handle.write(f"Cycle slip detected at epoch {epoch} for satellite {sat}: maximum gap "
                                       f"{self.max_gap * sampling_rate} seconds exceeded.\n")
                    cycle_slips[sat].set_data(epoch, True)
                    iLast_slip = iEpoch
                    prev_epoch = None
                    continue

                # fetch the observations for the last window_size epochs (moving window), after the previous cycle-slip
                indexes = [iEpoch - 1 - i for i in range(window_size) if iEpoch - 1 - i >= 0]
                indexes = [i for i in indexes if i > iLast_slip]  # remove previous cycle slip
                data = []
                time = []
                indexes.reverse()
                for i in indexes:
                    try:
                        obs_list = self.gf_data.get_observables_at_epoch(vEpochs[i], sat)
                        for obs in obs_list:
                            if "CP" in obs.datatype.data_type:
                                data.append(obs.value)
                                time.append((vEpochs[i] - epoch).total_seconds())
                                break
                    except TimeSeriesError:
                        continue

                # compute second-degree polynomial
                if len(time) > 2:
                    poly_coeff = np.polyfit(time, data, deg=2)
                    p = np.poly1d(poly_coeff)
                    y_fit = p(0)

                    if f_handle:
                        f_handle.write(f"{epoch}: {sat}: GF observable = {gf_obs}. Indexes = {indexes}. "
                                       f"Difference: {gf_obs.value - y_fit} m. Threshold: {threshold} m.\n")

                    if gf_obs.value - y_fit > threshold:
                        # cycle slip detected
                        if f_handle:
                            f_handle.write(f"Cycle slip detected at epoch {epoch} for satellite {sat}: "
                                           f"difference gf_obs.value - y_fit {gf_obs.value - y_fit} exceeds "
                                           f"threshold {threshold}\n")
                        cycle_slips[sat].set_data(epoch, True)
                        iLast_slip = iEpoch
                        prev_epoch = None
                        continue
                prev_epoch = epoch

    def _merge_results(self, cycle_slip_mw, cycle_slip_gf):
        merged_cycle_slip = dict()
        sat_list = []
        if cycle_slip_mw is not None:
            sat_list = list(cycle_slip_mw.keys())
        if cycle_slip_gf is not None:
            _list = cycle_slip_gf.keys()
            for sat in _list:
                if sat not in sat_list:
                    sat_list.append(sat)

        for sat in sat_list:
            merged_cycle_slip[sat] = TimeSeries()

            # get TimeSeries
            mw_data = cycle_slip_mw.get(sat, None)
            gf_data = cycle_slip_gf.get(sat, None)

            for cycle_slip_data in (mw_data, gf_data):
                if isinstance(cycle_slip_data, TimeSeries):
                    if not cycle_slip_data.is_empty():
                        vEpochs = cycle_slip_data.get_all_epochs()
                        for epoch in vEpochs:
                            if not merged_cycle_slip[sat].has_epoch(epoch):
                                merged_cycle_slip[sat].set_data(epoch, cycle_slip_data.get_data_for_epoch(epoch))

        for sat in sat_list:
            if merged_cycle_slip[sat].is_empty():
                merged_cycle_slip.pop(sat)

        self.cycle_slips = merged_cycle_slip

    def _get_gf_threshold(self, sat, sampling_rate):
        vTypes = self.gf_data.get_types(sat.sat_system)
        wavelength = SPEED_OF_LIGHT / vTypes[0].freq_value
        a0 = 1.5 * wavelength  # Maximum threshold
        a1 = a0 / 2
        T0 = 60  # seconds
        Delta_t = 1 / sampling_rate  # seconds
        threshold = a0 - a1 * np.exp(-Delta_t / T0)
        return threshold
