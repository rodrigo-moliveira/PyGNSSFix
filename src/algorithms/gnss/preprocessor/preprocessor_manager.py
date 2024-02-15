from src.common_log import get_logger
from src.data_types.gnss.data_type import data_type_from_rinex
from src.io.config import config_dict
from src.errors import PreprocessorError
from src.io.config.enums import EnumPositioningMode
from .filter import FilterMapper, TypeConsistencyFilter, RateDowngradeFilter, SignalCheckFilter
from .filter.ura_health_check import SatFilterHealthURA
from .functor import FunctorMapper, IonoFreeFunctor, SmoothFunctor


class PreprocessorManager:

    def __init__(self, trace_path, data_manager):
        self.log = get_logger("PREPROCESSOR")
        self.trace_path = trace_path
        self.data_manager = data_manager
        self.services = config_dict.get_services()

    def compute(self):
        """
        Algorithms to apply:
            Acting on Raw ObservationData:
                * Apply SNR Check Filter -> remove observables with low SNR
                * Type Consistency Filter -> remove all unnecessary observations and empty satellites

            Creating Processed ObservationData:
                * Compute IonoFree Observation Data -> Compute iono-free observables from the raw gnss_obs data,
                                                        creating a new dataset
                * Compute (IonoFree) Smooth Observation Data -> Compute smooth code observables

        """
        self.log.info("Starting Preprocessor...")

        # SNR Check Filter
        # TODO: write a filter report (x% of data was removed, etc.)
        # TODO: add log messages to all filters... (removed obs for epoch...)
        obs_data = self.data_manager.get_data("obs_data")
        try:
            self.snr_filter(obs_data)
        except Exception as e:
            raise PreprocessorError(f"PreprocessorManager -> Error performing SNR filter: {e}")

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
        compute_iono_free = (config_dict.get_model() == EnumPositioningMode.SPS_IF)
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
        self.log.debug(
            "Writing Preprocessor output Observation Data to trace file {}".format("PreprocessedObservationData.txt"))
        f = open(self.trace_path + "/PreprocessedObservationData.txt", "w")
        f.write(str(obs_data))
        f.close()

        self.log.info("End of module Preprocessor")

    def consistency_filter(self, observation_data):
        self.log.info("Applying consistency filter to remove unnecessary datatypes and data-less satellites")
        datatypes = {}

        # NOTE: currently we filter out Signal and Doppler Observables -> only need Pseudorange and Carrier Phase
        for constellation, services in self.services.items():
            datatypes[constellation] = []
            for service in services:
                pr = data_type_from_rinex(f"C{service}", constellation)
                if pr is not None:
                    datatypes[constellation].append(pr)
                cp = data_type_from_rinex(f"L{service}", constellation)
                if cp is not None:
                    datatypes[constellation].append(cp)

        type_filter = TypeConsistencyFilter(datatypes)
        mapper = FilterMapper(type_filter)
        mapper.apply(observation_data)

        # Saving Consistent data to file
        self.log.debug(
            "Writing Type Consistent Observation Data to trace file {}".format("TypeConsistentObservationData.txt"))
        f = open(self.trace_path + "/TypeConsistentObservationData.txt", "w")
        f.write(str(observation_data))
        f.close()

    def snr_filter(self, observation_data):
        snr_threshold = config_dict.get("preprocessor", "snr_filter")
        self.log.info(f"Applying SNR check filter. SNR Threshold is {snr_threshold}")

        snr_filter = SignalCheckFilter(observation_data, snr_threshold)
        mapper = FilterMapper(snr_filter)
        mapper.apply(observation_data)

        # Saving debug data to file
        self.log.debug(
            "Writing SNR Checked Observation Data to trace file {}".format("SNRCheckObservationData.txt"))
        f = open(self.trace_path + "/SNRCheckObservationData.txt", "w")
        f.write(str(observation_data))
        f.close()

    def sv_ura_health_filter(self, observation_data):
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
                                        gal_sisa_check, gal_sisa_val, gal_health, self.log)
        mapper = FilterMapper(ura_filter)
        mapper.apply(observation_data)

        # Saving debug data to file
        self.log.debug(
            "Writing SV URA and Health Check Observation Data to trace file {}".
            format("SvURAHealthObservationData.txt"))
        f = open(self.trace_path + "/SvURAHealthObservationData.txt", "w")
        f.write(str(observation_data))
        f.close()

    def iono_free(self, raw_data, data_out, constellation):
        functor = IonoFreeFunctor(constellation, self.services[constellation])
        mapper = FunctorMapper(functor)
        mapper.apply(raw_data, data_out)

    def smooth(self, data):
        time_constant = config_dict.get("preprocessor", "smooth_time_constant_secs")
        rate = data.get_rate()
        self.log.info(f"Computing smooth data with time constant set to {time_constant}[s]. Data rate is {rate}")

        smooth_functor = SmoothFunctor(time_constant, data.get_rate())
        mapper = FunctorMapper(smooth_functor)
        smooth_data = self.data_manager.get_data("smooth_obs_data")
        mapper.apply(data, smooth_data)

        self.log.debug("Writing Smooth Observation Data to trace file {}".format("SmoothObservationData.txt"))
        f = open(self.trace_path + "/SmoothObservationData.txt", "w")
        f.write(str(smooth_data))
        f.close()

    def downgrade(self, data):
        # Downgrade gnss_obs data rate
        rate_out = config_dict.get("inputs", "rate")
        rate_in = data.get_rate()

        if rate_in != rate_out:
            self.log.info(f"Downgrading gnss_obs data from input rate {rate_in} [s] to output rate "
                          f"{rate_out} [s]")

            if rate_out % rate_in != 0:
                self.log.warning(f"It is not possible to downgrade to the selected rate. Output rate is not divisible "
                                 f"by input rate. Keeping input rate of {rate_in}...")

            else:
                epoch_list = data.get_epochs()
                downgrade_filter = RateDowngradeFilter(rate_out, epoch_list[0])
                mapper = FilterMapper(downgrade_filter)
                mapper.apply(data)

        return data
