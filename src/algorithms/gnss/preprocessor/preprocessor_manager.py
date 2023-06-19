from src.common_log import get_logger
from src.data_types.gnss.data_type import data_type_from_rinex
from src.io.config import config_dict
from src.data_types.gnss.observation_data import ObservationData
from src.errors import PreprocessorError
from .filter import FilterMapper, TypeConsistencyFilter, RateDowngradeFilter, SignalCheckFilter
from src.algorithms.gnss.preprocessor.functor.constellation_filter import ConstellationFunctor
from .functor import FunctorMapper, IonoFreeFunctor


class PreprocessorManager:

    def __init__(self, trace_path, raw_data):
        self.log = get_logger("PREPROCESSOR")
        self.trace_path = trace_path
        self.raw_data = raw_data
        self.services = config_dict.get_services()

    def compute(self):
        """
        Algorithms to apply:
            Acting on Raw ObservationData:
                * Apply SNR Check Filter -> remove observables with low SNR
                * Type Consistency Filter -> remove all unnecessary observations and empty satellites

            Creating Processed ObservationData:
                * Compute IonoFree Observation Data -> Compute iono-free observables from the raw observation data,
                                                        creating a new dataset
                * Compute (IonoFree) Smooth Observation Data -> Compute smooth code observables

        """
        self.log.info("Starting Preprocessor...")

        # SNR Check Filter
        # TODO: write a filter report (x% of data was removed, etc.)
        try:
            self.snr_filter(self.raw_data)
        except Exception as e:
            raise PreprocessorError(f"PreprocessorManager -> Error performing SNR filter: {e}")

        # Type Consistency Filter
        try:
            self.consistency_filter(self.raw_data)
        except Exception as e:
            raise PreprocessorError(f"PreprocessorManager -> Error performing Consistency Type filter: {e}")

        obs_data_out = ObservationData()

        # check to compute or not iono free dataset from raw observables
        for constellation in self.services.keys():
            compute_iono_free = config_dict.is_iono_free(constellation)
            if compute_iono_free:
                obs_list = config_dict.get("model", constellation, "observations")
                n_obs = len(obs_list)
                if n_obs != 2:
                    raise PreprocessorError(f"Unable to compute iono free data for constellation {constellation} "
                                            f"due to lack of data. Need 2 observations, but have {n_obs} ({obs_list})")
                else:
                    self.log.info(f"Computing iono free data for constellation {constellation} with observations "
                                  f"{obs_list}")
                    try:
                        self.iono_free(self.raw_data, obs_data_out, constellation)
                    except Exception as e:
                        raise PreprocessorError(
                            f"PreprocessorManager -> Error computing Iono Free Observation Data: {e}")
            else:
                self.log.info(f"Iono free data not computed for {constellation} due to user choice")
                # select raw observables for this constellation
                functor = ConstellationFunctor(constellation)
                mapper = FunctorMapper(functor)
                mapper.apply(self.raw_data, obs_data_out)

        """Currently, the Smooth Algorithm is turned off. Further investigation is needed to fix it"""
        # Get Smooth Observation Data
        # try:
        #    _data_out = self.smooth(_data_out)
        # except Exception as e:
        #    raise PreprocessorError(f"Error computing Smooth Observation Data: {e}")

        # Prepare ObservationData for output (downgrade output rate)
        try:
            self.downgrade(obs_data_out)
        except Exception as e:
            raise PreprocessorError(f"PreprocessorManager -> Error performing Downgrade Rate filter: {e}")

        # Saving Output Observation Data to File
        self.log.debug(
            "Writing Preprocessor output Observation Data to trace file {}".format("PreprocessedObservationData.txt"))
        f = open(self.trace_path + "/PreprocessedObservationData.txt", "w")
        f.write(str(obs_data_out))
        f.close()

        self.log.info("End of module Preprocessor")
        return obs_data_out

    def consistency_filter(self, observation_data):
        self.log.info("Applying consistency filter to remove unnecessary datatypes and data-less satellites")
        datatypes = {}

        # NOTE: currently we filter out Signal and Carrier Observables -> only need Pseudorange
        for constellation, services in self.services.items():
            datatypes[constellation] = []
            for service in services:
                datatype = data_type_from_rinex(f"C{service}", constellation)
                if datatype is not None:
                    datatypes[constellation].append(datatype)

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

    def iono_free(self, raw_data, data_out, constellation):
        functor = IonoFreeFunctor(constellation, self.services[constellation])
        mapper = FunctorMapper(functor)
        mapper.apply(raw_data, data_out)

    """def smooth(self, data):
        self.log.info("Computing smooth data")

        smooth_functor = SmoothFunctor(config["model"]["smooth"]["time_constant"],
                                       data.get_rate())
        mapper = FunctorMapper(smooth_functor)
        smooth_data = ObservationData()
        mapper.apply(data, smooth_data)

        self.log.debug("Writing Smooth Observation Data to trace file {}".format("SmoothObservationData.txt"))
        f = open(self.trace_path + "/SmoothObservationData.txt", "w")
        f.write(str(smooth_data))
        f.close()

        # change pointer of _data_out
        data = smooth_data
        return data"""

    def downgrade(self, data):
        # Downgrade observation data rate
        rate_out = config_dict.get("inputs", "rate")
        rate_in = data.get_rate()

        if rate_in != rate_out:
            self.log.info(f"Downgrading observation data from input rate {rate_in} [s] to output rate "
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
