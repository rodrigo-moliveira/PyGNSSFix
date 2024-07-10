import json
import sys
import os

from src import RUNS_PATH
from src.data_mng.csv.csv_data_mng import GnssRunStorageManager
from src.io.config import config_dict
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
        run_path = config_dict.get("performance_evaluation", "run")
        run_path = RUNS_PATH / run_path
        print(f"Running Post-Processing Performance Evaluation Module for run {run_path}")
    except IndexError:
        print("ERROR: No valid run path was provided in the configurations")
        print("Please check the field 'performance_evaluation.run' of the json file")
        exit(-1)

    # validate path
    if not os.path.exists(run_path):
        print(f"ERROR: Provided path {run_path} does not exist")
        exit(-1)

    # Launching post processing algorithm
    try:
        # load run
        data_manager = GnssRunStorageManager()
        data_manager.read_data(run_path / 'output')

        # run Performance Evaluation Module
        print("Executing Performance Evaluation Manager for GNSS run...")
        eval_manager = PerformanceManager(data_manager, config_dict)
        eval_manager.process(run_path)

    except Exception as e:
        print(f"Unexpected error running while running program: {e}")
        exit()


print("#--------------------------------------------------#")
print("#           Welcome to GNSSNavPy Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
