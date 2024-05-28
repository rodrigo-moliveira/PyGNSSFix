import os

from src.common_log import get_logger, PREPROCESSOR_LOG
from src.data_types.gnss import data_type_from_rinex
from src.data_types.gnss.service_utils import get_freq_from_service
from src.io.config import config_dict, EnumPositioningMode
from src.errors import PreprocessorError
from .filter import *
from .functor import *


class PreprocessorManager:
    """
    PreprocessorManager class: performs the preprocessing operations to the GNSS measurements
    to prepare them for the GNSS Solver

    The operations involved are:
        * Filtering: clean and delete the observation dataset according to some checks:
            - SNR Filter: remove measurements with low SNR (signal-to-noise ratio)
            - Satellite Health Filter: remove measurements from unhealthy satellites (check GPS URA or GAL SISA)
            - Type Consistency Filter: clean satellite epoch data for missing observables
            - Downgrade Filter: downgrade the observation dataset for the provided user-defined rate

        * Creating Processed ObservationData:
                - Compute IonoFree observation dataset (only if defined in the user configurations)
                - Compute Smooth and IonoFree Smooth observation dataset (only if defined in the user configurations)

        * Write intermediate trace files for each operation
    """

    def __init__(self, trace_path, data_manager):
        """
        Constructor of the PreprocessorManager instance

        Parameters:
            trace_path(str): path to store the Preprocessor trace files
            data_manager(src.data_mng.gnss.gnss_data_mng.GnssDataManager): GNSS data manager instance
        """
        self.log = get_logger(PREPROCESSOR_LOG)
        self.data_manager = data_manager
        self.services = config_dict.get_services()

        # set up trace path
        self.write_trace = config_dict.get("preprocessor", "trace_files")
        self.trace_path = None
        if self.write_trace:
            self.trace_path = f"{trace_path}\\preprocessor"
            try:
                os.makedirs(self.trace_path)
            except:
                raise IOError(f"Cannot create dir: {self.trace_path}")

    def compute(self):
        """
        Launch the Preprocessor algorithm

        Raises:
            PreprocessorError: if an error is detected during the Preprocessor algorithm, an error is raised.
        """
        self.log.info("Starting Preprocessor...")

        # SNR Filter
        obs_data = self.data_manager.get_data("obs_data")
        try:
            self.snr_filter(obs_data)
        except Exception as e:
            raise PreprocessorError(f"PreprocessorManager -> Error performing SNR filter: {e}")

        # Satellite Health Check (SiS reports)
        try:
            self.sv_ura_health_filter(obs_data)
        except Exception as e:
            raise PreprocessorError(f"PreprocessorManager -> Error performing SV URA and Health filter: {e}")

        # Type Consistency Filter
        try:
            self.consistency_filter(obs_data)
        except Exception as e:
            raise PreprocessorError(f"PreprocessorManager -> Error performing Consistency Type filter: {e}")

        # Prepare ObservationData for output (downgrade output rate)
        try:
            self.downgrade(obs_data)
        except Exception as e:
            raise PreprocessorError(f"PreprocessorManager -> Error performing Downgrade Rate filter: {e}")

        # check to compute or not iono free dataset from raw observables
        model = config_dict.get("model", "mode")
        compute_iono_free = (model == EnumPositioningMode.SPS_IF)
        if compute_iono_free:
            iono_free_data = self.data_manager.get_data("iono_free_obs_data")
            for constellation in self.services.keys():
                obs_list = config_dict.get("model", constellation, "observations")
                n_obs = len(obs_list)
                if n_obs != 2:
                    raise PreprocessorError(f"Unable to compute iono free data for constellation {constellation} "
                                            f"due to lack of data. Need 2 observations, but have {n_obs} ({obs_list})")
                else:
                    self.log.info(f"Computing iono free data for constellation {constellation} with observations "
                                  f"{obs_list}")
                    try:
                        self.iono_free(obs_data, iono_free_data, constellation)
                    except Exception as e:
                        raise PreprocessorError(
                            f"PreprocessorManager -> Error computing Iono Free Observation Data: {e}")
            if self.write_trace:
                self.log.debug(
                    "Writing IonoFree Observation Data to trace file {}".format("IonoFreeObservationData.txt"))
                f = open(self.trace_path + "/IonoFreeObservationData.txt", "w")
                f.write(str(iono_free_data))
                f.close()
        else:
            self.log.info(f"Iono free data not computed due to user choice")
            # select raw observables for this constellation
            # constellation_filter = ConstellationFilter(self.services.keys())
            # mapper = FilterMapper(constellation_filter)
            # mapper.apply(obs_data)

        # Get Smooth Observation Data
        compute_smooth = config_dict.get("preprocessor", "compute_smooth")
        if compute_smooth:
            try:
                data = self.data_manager.get_data("iono_free_obs_data") if compute_iono_free else obs_data
                self.smooth(data)
            except Exception as e:
                raise PreprocessorError(f"Error computing Smooth Observation Data: {e}")

        # Saving Output Observation Data to File
        if self.write_trace:
            self.log.debug(
                "Writing Preprocessor output Observation Data to trace file {}".format(
                    "PreprocessedObservationData.txt"))
            f = open(self.trace_path + "/PreprocessedObservationData.txt", "w")
            f.write(str(obs_data))
            f.close()

        self.log.info("End of module Preprocessor")

    def consistency_filter(self, observation_data):
        """ Perform Consistency Filter
        """
        self.log.info("Applying consistency filter to remove unnecessary datatypes and data-less satellites")
        required_datatypes = {}

        keep_doppler = config_dict.get("model", "estimate_velocity")  # doppler only needed when velocity is estimated
        keep_carrier = config_dict.get("preprocessor", "compute_smooth")  # carrier only needed for smoothing

        for constellation, services in self.services.items():
            required_datatypes[constellation] = []
            for service in services:
                pr = data_type_from_rinex(f"C{service}", constellation)
                required_datatypes[constellation].append(pr)
                if keep_carrier:
                    cp = data_type_from_rinex(f"L{service}", constellation)
                    required_datatypes[constellation].append(cp)
                if keep_doppler:
                    doppler = data_type_from_rinex(f"D{service}", constellation)
                    required_datatypes[constellation].append(doppler)

        self.log.info(f"Type Consistency Filter: Required datatypes are {required_datatypes}")
        type_filter = TypeConsistencyFilter(required_datatypes, self.trace_path)
        mapper = FilterMapper(type_filter)
        mapper.apply(observation_data)

        # Write report to log
        self.log.info(f"Type Consistency Filter Report: {mapper.report}")

        # Saving Consistent data to file
        # if self.write_trace:
        #    self.log.debug(
        #       "Writing Type Consistent Observation Data to trace file {}".format("TypeConsistentObservationData.txt"))
        #    f = open(self.trace_path + "/TypeConsistentObservationData.txt", "w")
        #    f.write(str(observation_data))
        #    f.close()

    def snr_filter(self, observation_data):
        """ Perform SNR Filter """
        snr_threshold = config_dict.get("preprocessor", "snr_filter")
        self.log.info(f"Applying SNR check filter. SNR Threshold is {snr_threshold}")

        snr_filter = SignalCheckFilter(observation_data, snr_threshold, self.trace_path)
        mapper = FilterMapper(snr_filter)
        mapper.apply(observation_data)

        # Write report to log
        self.log.info(f"SNR Filter Report: {mapper.report}")

        # Saving debug data to file
        # if self.write_trace:
        #    self.log.debug(
        #        "Writing SNR Checked Observation Data to trace file {}".format("SNRCheckObservationData.txt"))
        #    f = open(self.trace_path + "/SNRCheckObservationData.txt", "w")
        #    f.write(str(observation_data))
        #    f.close()

    def sv_ura_health_filter(self, observation_data):
        """ Perform Satellite Health Check (GPS URA and GAL SISA test) """
        gps_ura_check = config_dict.get("preprocessor", "satellite_status", "GPS", "URA")
        gps_ura_val = config_dict.get("preprocessor", "satellite_status", "GPS", "max_URA")
        gps_health = config_dict.get("preprocessor", "satellite_status", "GPS", "health")

        gal_sisa_check = config_dict.get("preprocessor", "satellite_status", "GAL", "SISA")
        gal_sisa_val = config_dict.get("preprocessor", "satellite_status", "GAL", "max_SISA")
        gal_health = config_dict.get("preprocessor", "satellite_status", "GAL", "health")

        self.log.info(f"Applying GPS URA, GAL SISA and Health Status checks. "
                      f"GPS URA filter = {gps_ura_check} URA threshold = {gps_ura_val}m, health status check is "
                      f"{gps_health}. "
                      f"GAL SISA filter = {gal_sisa_check} SISA threshold = {gal_sisa_val}m, health status check is "
                      f"{gal_health}")

        nav_data = self.data_manager.get_data("nav_data")
        ura_filter = SatFilterHealthURA(nav_data, gps_ura_check, gps_ura_val, gps_health,
                                        gal_sisa_check, gal_sisa_val, gal_health, self.log, self.trace_path)
        mapper = FilterMapper(ura_filter)
        mapper.apply(observation_data)

        # Write report to log
        self.log.info(f"GPS URA, GAL SISA and Health Status Filter Report: {mapper.report}")

        # Saving debug data to file
        # if self.write_trace:
        #    self.log.debug(
        #        "Writing SV URA and Health Check Observation Data to trace file {}".
        #        format("SvURAHealthObservationData.txt"))
        #    f = open(self.trace_path + "/SvURAHealthObservationData.txt", "w")
        #    f.write(str(observation_data))
        #    f.close()

    def iono_free(self, raw_data, data_out, constellation):
        """ Compute IonoFree data """
        services = self.services[constellation]
        if len(services) != 2:
            raise AttributeError(f"Problem getting base and second frequencies in Iono Free Computation for "
                                 f"constellation {constellation}, observations provided are {services}. "
                                 f"There should be 2 observation types defined, and not {len(services)}")
        base_freq = get_freq_from_service(services[0], constellation)
        second_freq = get_freq_from_service(services[1], constellation)
        if base_freq is not None and second_freq is not None:
            functor = IonoFreeFunctor(constellation, base_freq, second_freq)
            mapper = FunctorMapper(functor)
            mapper.apply(raw_data, data_out)
        else:
            raise AttributeError(f"Unable to fetch base ({base_freq}) and second ({second_freq}) frequencies "
                                 f"for constellation {constellation} and services {services}. Please review"
                                 f" the input observation data and the configurations.")

    def smooth(self, data):
        """ Compute Smooth data """
        time_constant = config_dict.get("preprocessor", "smooth_time_constant_secs")
        rate = data.get_rate()
        self.log.info(f"Computing smooth data with time constant set to {time_constant}[s]. Data rate is {rate}")

        smooth_functor = SmoothFunctor(time_constant, rate)
        mapper = FunctorMapper(smooth_functor)
        smooth_data = self.data_manager.get_data("smooth_obs_data")
        mapper.apply(data, smooth_data)

        if self.write_trace:
            self.log.debug("Writing Smooth Observation Data to trace file {}".format("SmoothObservationData.txt"))
            f = open(self.trace_path + "/SmoothObservationData.txt", "w")
            f.write(str(smooth_data))
            f.close()

    def downgrade(self, data):
        """ Perform Downgrading of the GNSS measurements input rate """
        # Downgrade gnss_models data rate
        rate_out = config_dict.get("inputs", "rate")
        rate_in = data.get_rate()

        if rate_in != rate_out:
            self.log.info(f"Downgrading gnss_models data from input rate {rate_in} [s] to output rate "
                          f"{rate_out} [s]")

            if rate_out % rate_in != 0:
                self.log.warning(f"It is not possible to downgrade to the selected rate. Output rate is not divisible "
                                 f"by input rate. Keeping input rate of {rate_in}...")

            else:
                epoch_list = data.get_epochs()
                downgrade_filter = RateDowngradeFilter(rate_out, epoch_list[0], self.trace_path)
                mapper = FilterMapper(downgrade_filter)
                mapper.apply(data)

                # Write report to log
                self.log.info(f"Rate Downgrade Filter Report: {mapper.report}")
