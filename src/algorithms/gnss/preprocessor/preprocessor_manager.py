from .... import get_logger
from ....config import config
from PositioningSolver.src.gnss.data_types.ObservationData import ObservationData
from ....gnss.data_types.ServicesUtils import get_code_type_from_service
from ....utils.errors import PreprocessorError
from .filter import FilterMapper, TypeConsistencyFilter, RateDowngradeFilter, SignalCheckFilter
from .functor import FunctorMapper, IonoFreeFunctor, SmoothFunctor


class Preprocessor:

    def __init__(self, trace_path, service_manager, constellation, raw_data, compute_iono_free, output_rate):

        log = get_logger("preprocessor")
        log.info("###############################################################")
        log.info("###### Starting module 'Processing Observation Data' ... ######")

        self.log = log
        self.trace_path = trace_path
        self.service_manager = service_manager
        self.constellation = constellation
        self.raw_data = raw_data

        self.compute_iono_free = compute_iono_free
        self.output_rate = output_rate

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
        # SNR Check Filter
        try:
            self.snr_filter(self.raw_data)
        except Exception as e:
            raise PreprocessorError(f"Error performing SNR filter: {e}")

        # Type Consistency Filter
        try:
            self.consistency_filter(self.raw_data)
        except Exception as e:
            raise PreprocessorError(f"Error performing Consistency Type filter: {e}")

        # pointer to output dataset (processed by Iono, Smooth and RateDowngrade functors)
        _data_out = self.raw_data

        # Get Iono Free Observation Data
        try:
            _data_out = self.iono_free(_data_out)
        except Exception as e:
            raise PreprocessorError(f"Error computing Iono Free Observation Data: {e}")

        """Currently, the Smooth Algorithm is turned off. Further investigation is needed to fix it"""
        # Get Smooth Observation Data
        # try:
        #    _data_out = self.smooth(_data_out)
        # except Exception as e:
        #    raise PreprocessorError(f"Error computing Smooth Observation Data: {e}")

        # Prepare ObservationData for output (downgrade output rate)
        try:
            _data_out = self.downgrade(_data_out)
        except Exception as e:
            raise PreprocessorError(f"Error performing Downgrade Rate filter: {e}")

        self.log.info("####### End of module 'Process Observation Data' ... #######\n")

        return _data_out

    def consistency_filter(self, observation_data):
        self.log.info("Applying consistency filter to remove unnecessary datatypes and data-less satellites")
        types = get_code_type_from_service(self.service_manager.services[self.constellation], self.constellation)

        # if self.compute_smooth:
        #    types += get_carrier_type_from_service(self.service_manager.services[self.constellation],
        #                                           self.constellation)

        type_functor = TypeConsistencyFilter(types)
        mapper = FilterMapper(type_functor)
        mapper.apply(observation_data)

        # Saving Consistent data to file
        self.log.debug(
            "Writing Type Consistent Observation Data to trace file {}".format("TypeConsistentObservationData.txt"))
        f = open(self.trace_path + "/TypeConsistentObservationData.txt", "w")
        f.write(str(observation_data))
        f.close()

    def snr_filter(self, observation_data):
        self.log.info("Applying SNR check filter to remove data with low signal to noise ratio (SNR)")

        snr_functor = SignalCheckFilter(self.raw_data)
        mapper = FilterMapper(snr_functor)
        mapper.apply(observation_data)

        # Saving Consistent data to file
        self.log.debug(
            "Writing SNR Checked Observation Data to trace file {}".format("SNRCheckObservationData.txt"))
        f = open(self.trace_path + "/SNRCheckObservationData.txt", "w")
        f.write(str(observation_data))
        f.close()

    def iono_free(self, data):
        if self.compute_iono_free:
            self.log.info("Computing iono free data")
            iono_free_functor = IonoFreeFunctor(self.service_manager[self.constellation])
            mapper = FunctorMapper(iono_free_functor)
            iono_free_data = ObservationData()
            mapper.apply(data, iono_free_data)

            # Saving Iono Free data to file
            self.log.debug("Writing Iono Free Observation Data to trace file {}".format("IonoFreeObservationData.txt"))
            f = open(self.trace_path + "/IonoFreeObservationData.txt", "w")
            f.write(str(iono_free_data))
            f.close()

            # change pointer of _data_out
            data = iono_free_data

        return data

    def smooth(self, data):
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
        return data

    def downgrade(self, data):
        # Downgrade observation data rate
        if self.output_rate:
            rate_in = data.get_rate()
            self.log.info(f"Downgrading observation data from input rate {rate_in} [s] to rate "
                          f"{self.output_rate} [s]")

            if self.output_rate % rate_in != 0:
                self.log.warning(f"It is not possible to downgrade to the selected rate. Output rate is not divisible "
                                 f"by input rate! Keeping input rate of {rate_in}...")

            else:
                epoch_list = data.get_epochs()
                downgrade_filter = RateDowngradeFilter(self.output_rate, epoch_list[0])
                mapper = FilterMapper(downgrade_filter)
                mapper.apply(data)

                self.log.debug(
                    "Writing Downgraded Observation Data to trace file {}".format("DowngradedObservationData.txt"))
                f = open(self.trace_path + "/DowngradedObservationData.txt", "w")
                f.write(str(data))
                f.close()

        return data
