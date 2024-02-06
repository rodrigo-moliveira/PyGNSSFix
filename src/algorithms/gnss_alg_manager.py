import os
import time
import traceback

from src import RUNS_PATH
from src.data_mng.gnss_data_mng import GnssDataManager
from src.common_log import set_logs, get_logger
from src.data_types.gnss.observation_data import ObservationData
from src.io.rinex.nav_reader import RinexNavReader
from src.io.rinex.obs_reader import RinexObsReader
from src.data_types.gnss.navigation_data import NavigationData
from src.io.config import config_dict


class GnssAlgorithmManager:

    def __init__(self, algorithm):
        # create data members
        self.data_manager = GnssDataManager()
        self.algorithm = algorithm

        # create output folder
        data_dir = config_dict.get("output", "output_path")
        self.data_dir = self._check_data_dir(data_dir)

        # initialize logger objects
        set_logs(config_dict.get("log", "minimum_level"), f"{self.data_dir}\\log.txt")

        main_log = get_logger("MAIN_LOG")
        main_log.info("Starting GNSS Algorithm Manager")

    def _read_inputs(self):

        logger = get_logger("IO")
        #TODO: add log messages
        #TODO: add possibility of multiple files

        # read navigation data
        nav_file = config_dict.get("inputs", "nav_files")[0]
        obs_file = config_dict.get("inputs", "obs_files")[0]
        services = config_dict.get_services()
        first_epoch = config_dict.get("inputs", "arc", "first_epoch")
        last_epoch = config_dict.get("inputs", "arc", "last_epoch")
        snr_check = config_dict.get("inputs", "snr_control")

        RinexNavReader(nav_file, self.data_manager.get_data("nav_data"))
        RinexObsReader(self.data_manager.get_data("obs_data"), obs_file, services, logger, first_epoch, last_epoch,
                       snr_check)

        # trace data files
        trace_dir = f"{self.data_dir}\\trace"
        with open(f"{trace_dir}\\RawObservationData.txt", "w") as file:
            file.write(str(self.data_manager.get_data("obs_data")))
        with open(f"{trace_dir}\\RawNavigationData.txt", "w") as file:
            file.write(str(self.data_manager.get_data("nav_data")))

    def run(self):
        main_log = get_logger("MAIN_LOG")
        main_log.info(f"Running algorithm {str(self.algorithm)}")

        # Input Reader Module
        try:
            main_log.info(f"Starting Input Reader Module...")
            self._read_inputs()
        except Exception as e:
            main_log.error(f"Stopping execution of program due to error in execution of Input Reader Module: {e}")
            print(traceback.format_exc())
            exit(-1)

        # Main Algorithm Module
        try:
            main_log.info(f"Starting Main Algorithm Module...")
            self.algorithm.compute(self.data_manager, f"{self.data_dir}\\trace")
        except Exception as e:
            main_log.error(f"Stopping execution of program due to error in execution of Main Algorithm Module: {e}")
            print(traceback.format_exc())
            exit(-1)

        # Output Writer Module
        try:
            main_log.info(f"Starting Output Writer Module...")
            self._save_run()
        except Exception as e:
            main_log.error(f"Stopping execution of program due to error in execution of Output Writer Module: {str(e)}")
            print(traceback.format_exc())
            exit(-1)

        main_log.info(f"Successfully executed algorithm {str(self.algorithm)}")

    def _save_run(self):
        log = get_logger("IO")
        self.data_manager.save_data(f"{self.data_dir}\\output", log)

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
        if data_dir is None or data_dir == '' or data_dir == "default":
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
                os.makedirs(f"{data_dir}\\output")
            except:
                raise IOError(f"Cannot create dir: {data_dir}")
        return data_dir
