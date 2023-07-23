import json
import sys
import os

from src import RUNS_PATH
from src.io.config import config_dict


def main():
    try:
        run_path = sys.argv[1]
        run_path = RUNS_PATH / run_path
        print(f"Running Performance Evaluation Module for run {run_path}")
    except IndexError as e:
        print("ERROR: No run path was provided as command argument")
        print("To run `performance_eval.py` please do:\n\t$ python performance_eval.py <path_to_run>")
        exit(-1)

    # validate path
    if not os.path.exists(run_path):
        print(f"ERROR: Provided path {run_path} does not exist")
        exit(-1)

    # read config file
    try:
        print(f"Reading config file {run_path/'config.json'}")
        with open(run_path/"config.json") as json_file:
            data = json.load(json_file)
        config_dict.init(data)
    except Exception as e:
        print(f"Error Reading Configuration File: {e}\nfilename = {run_path/'config.json'}")
        exit()

    try:
        # load run


        # run Performance Evaluation Module
        alg_mng.run()

        
    except Exception as e:
        print(f"Unexpected error running while running program: {e}")
        exit()



print("#--------------------------------------------------#")
print("#           Welcome to PyGNSSFix Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
