""" Main scrip to execute the GNSS Post Processing / Performance Manager algorithm.

The post_processing_gnss.py` script launches the execution of the GNSS Post Processing algorithm.

This algorithm is responsible for evaluating the performance of the PNT estimation:
    * Computes the estimation positioning errors
    * Performs statistical tests on the Normal Equations post-fit residuals (residual analysis)
    * Plots the estimated states, errors and additional data

To run `post_processing_gnss.py` please do:
    $ python post_processing_gnss.py <path_to_config>

where:
    <path_to_config>: Path to the configuration file. Examples of configuration files are stored in the project folder
        :file:`../configs`

Usage:
    To execute the script, pass the path to the configuration file:

        $ python post_processing_gnss.py configs/sample_config.xml
"""
import json
import sys
import os

from src import RUNS_PATH
from src.common_log import set_logs, get_logger, POST_PROC_LOG
from src.data_mng.csv.csv_data_mng import GnssRunStorageManager
from src.io.config import config_dict
from src.io.io_utils import clean_text_file, get_last_run_folder
from src.modules.performance.performance_manager import PerformanceManager


def main():
    # get json config file from command line args
    try:
        config_filename = sys.argv[1]
    except IndexError:
        print("ERROR: No configuration json file was provided as command argument")
        print("To run `post_processing_gnss.py` please do:\n\t$ python post_processing_gnss.py <path_to_config>")
        exit(-1)

    # read config file
    try:
        print(f"Reading config file {config_filename}")
        with open(config_filename) as json_file:
            data = json.load(json_file)
        config_dict.init(data, alg="post_processing")
    except Exception as e:
        print(f"Error Reading Configuration File: {e}\nfilename = {config_filename}")
        exit()

    # Fetching run path folder
    try:
        run_path = config_dict.get("performance_evaluation", "run_folder")
        if run_path == "last":
            run_path = get_last_run_folder()

        run_path = RUNS_PATH / run_path
    except IndexError:
        print("ERROR: No valid run path was provided in the configurations")
        print("Please check the field 'performance_evaluation.run' of the json file")
        exit(-1)

    # validate path
    if not os.path.exists(run_path):
        print(f"ERROR: Provided path {run_path} does not exist")
        exit(-1)

    # initialize logger object
    set_logs(config_dict.get("log", "minimum_level"), f"{run_path}\\log_post_proc.txt")
    clean_text_file(f"{run_path}\\log_post_proc.txt")
    log = get_logger(POST_PROC_LOG)

    # Launching post processing algorithm
    try:
        log.info(f"Starting Post Processing Script for GNSS Run in {run_path}")

        # load run
        data_manager = GnssRunStorageManager(log)
        data_manager.read_data(run_path / 'output')

        # run Performance Evaluation Module
        log.info("Executing the GNSS Performance Evaluation Manager...")
        eval_manager = PerformanceManager(data_manager, log)
        eval_manager.process(run_path)

        log.info("Exiting Post Processing Script!")
        eval_manager.show_plots()

    except Exception as e:
        log.error(f"Unexpected error running while running program: {e}")
        exit()


print("#--------------------------------------------------#")
print("#           Welcome to GNSSNavPy Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
