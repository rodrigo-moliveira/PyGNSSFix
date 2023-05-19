import os
import time

from src import RUNS_PATH, WORKSPACE_PATH
from src.data_mng.gnss_data_mng import GnssDataManager
from src.common_log import set_logs, get_logger
from src.data_types.gnss.observation_data import ObservationData
from src.errors import PyGNSSFixError
from src.io.rinex.obs_reader import RinexObsReader


class GnssAlgorithmManager:

    def __init__(self, algorithm, config):
        """
        Add data to available.
        Args:
            algorithm (src.algorithms.algorithm.Algorithm) : name..
            config (src.io.config.gnss_config.ConfigGNSS) : name..
        """

        self.data_manager = GnssDataManager()
        self.algorithm = algorithm
        self.config = config

        # create output folder
        data_dir = self.config.performance_evaluation.output_path
        self.data_dir = self._check_data_dir(data_dir)

        # creating logger object
        # initialize logger objects
        set_logs(config.log.log_level, f"{self.data_dir}\\log.txt")

    def _read_inputs(self, logger):

        try:
            # read navigation data
            nav_data_dir = self.config.inputs.rinex_nav
            nav_data = None
            # self.data_manager.add_data("navigation_data", nav_data)

            # get observation data variables
            obs_file = self.config.inputs.rinex_obs[0]
            services = self.config.get_services()
            first_epoch = self.config.inputs.first_epoch
            last_epoch = self.config.inputs.last_epoch
            snr_check = self.config.inputs.snr_control

            handler = RinexObsReader(obs_file, services, logger, first_epoch, last_epoch, snr_check)

            # self.data_manager.add_data("services", services)
            self.data_manager.add_data("obs_data", handler.obs)

            # ... add more here
        except PyGNSSFixError as e:
            logger.error(f"{str(e)}")
        except Exception as e:
            logger.error(f"Exception caught -> {str(e)}")

    def _run(self):
        pass

        """# fetch input variables to this algorithm
                input_names = self.algorithm.inputs
                inputs = []
    
                for _in_name in input_names:
                    _in = self.data_manager.get_data(_in_name)
                    inputs.append(_in)
    
                self.algorithm.compute(*inputs)
    
                # get results and add them to the data manager
                _results = self.algorithm.get_results()
                for _data, _name in zip(_results, self.algorithm.outputs):
                    if _data is not None:
                        self.data_manager.add_data(_name, _data)
                """

    def run(self):
        # get loggers
        main_log = get_logger("GNSS_ALG")
        io_log = get_logger("IO")
        preprocessor_log = get_logger("PREPROCESSOR")

        # start algorithm
        main_log.info(f"Running {str(self.algorithm)}")

        # read inputs
        self._read_inputs(io_log)

        # run algorithm

        # process results
        self._results()

    def _results(self, trace=True, performance=False, save=False):

        # save data files
        # if save:
        #    self.data_manager.save_data(data_dir)

        if trace:
            # store trace files
            trace_dir = f"{self.data_dir}\\trace\\observation_data.txt"
            with open(trace_dir, "w") as file:
                file.write(str(self.data_manager.get_data("obs_data")))

        if performance:
            pass
            # GnssQualityManager.process(self.data_manager, data_dir, self.algorithm.name, performance, plot,
            # separate_axis)

    def _check_data_dir(self, data_dir):
        """
        check if data_dir is a valid dir. If not, use the default dir.
        check if the data_dir exists. If not, create it.
        Args:
            data_dir: all generated files are saved in data_dir
        Returns:
            data_dir: valid data dir.
        """
        # check data dir
        # data_dir is not specified, automatically create one
        if data_dir is None or data_dir == '':
            data_dir = str(RUNS_PATH)
            if data_dir[-1] != '//':
                data_dir = data_dir + '//'
            data_dir = data_dir + time.strftime('%Y-%m-%dT%HH%MM%SS', time.localtime()) + '//'
            data_dir = os.path.abspath(data_dir)

        # try to create data dir
        if not os.path.exists(data_dir):
            try:
                data_dir = os.path.abspath(data_dir)
                os.makedirs(data_dir)
                os.makedirs(f"{data_dir}\\trace")
            except:
                raise IOError(f"Cannot create dir: {data_dir}")
        return data_dir
