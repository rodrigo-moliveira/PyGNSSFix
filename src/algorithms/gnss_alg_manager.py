import os
import time

from src import RUNS_PATH
from src.data_mng.gnss_data_mng import GnssDataManager
from src.common_log import set_logs, get_logger
from src.data_types.gnss.observation_data import ObservationData
from src.errors import PyGNSSFixError
from src.io.rinex.nav_reader import RinexNavReader
from src.io.rinex.obs_reader import RinexObsReader
from src.data_types.gnss.navigation_data import NavigationData
from src.io.config import config_dict


class GnssAlgorithmManager:

    def __init__(self, algorithm):
        """
        Add data to available.
        Args:
            algorithm (src.algorithms.algorithm.Algorithm) : name..
        """
        # create data members
        self.data_manager = GnssDataManager()
        self.algorithm = algorithm

        # create output folder
        data_dir = config_dict.get("performance_evaluation", "output_path")
        self.data_dir = self._check_data_dir(data_dir)

        # initialize logger objects
        set_logs(config_dict.get("log", "minimum_level"), f"{self.data_dir}\\log.txt")

    def _read_inputs(self):

        logger = get_logger("IO")
        # TODO: add log messages
        try:
            # read navigation data
            # limpar aqui
            nav_file = config.get_config().inputs.rinex_nav[0]
            obs_file = config.get_config().inputs.rinex_obs[0]
            services = config.get_config().get_services()
            first_epoch = config.get_config().inputs.first_epoch
            last_epoch = config.get_config().inputs.last_epoch
            snr_check = config.get_config().inputs.snr_control

            nav = NavigationData()
            obs = ObservationData()

            RinexNavReader(nav_file, nav)
            RinexObsReader(obs, obs_file, services, logger, first_epoch, last_epoch, snr_check)

            # self.data_manager.add_data("services", services)
            self.data_manager.add_data("obs_data", obs)
            self.data_manager.add_data("nav_data", nav)

            # ... add more here
        except PyGNSSFixError as e:
            logger.error(f"{str(e)}")
            return False
        except Exception as e:
            logger.error(f"Exception caught -> {str(e)}")
            return False

        return True

    def run(self):
        main_log = get_logger("GNSS_ALG")

        # start algorithm
        main_log.info(f"Running {str(self.algorithm)}")

        # read inputs
        main_log.info(f"Starting IO Module...")
        success = self._read_inputs()
        if not success:
            main_log.warn("Stopping execution of program due to error encountered in execution of IO Module")
            return

        # computing algorithm
        main_log.info(f"Starting Main Algorithm Module...")
        self.algorithm.compute(self.data_manager)

        # process results
        main_log.info(f"Starting Performance Module...")
        self._results()

        main_log.info(f"Successfully executed algorithm {str(self.algorithm)}")

    def _results(self, trace=True, performance=False, save=False):

        # save data files
        # if save:
        #    self.data_manager.save_data(data_dir)

        if trace:
            # store trace files
            trace_dir = f"{self.data_dir}\\trace"
            with open(f"{trace_dir}\\observation_data.txt", "w") as file:
                file.write(str(self.data_manager.get_data("obs_data")))
            with open(f"{trace_dir}\\navigation_data.txt", "w") as file:
                file.write(str(self.data_manager.get_data("nav_data")))

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
