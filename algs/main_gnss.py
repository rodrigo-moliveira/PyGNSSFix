import sys
import json
import shutil

from src.algorithms.gnss_alg_manager import GnssAlgorithmManager
from src.io.config import config_dict
from src.errors import ConfigError, ConfigTypeError
from src.algorithms.gnss_sps import GnssSinglePointSolution


def main():
    # get json config file from command line args
    try:
        config_filename = sys.argv[1]
    except IndexError as e:
        print("ERROR: No configuration json file was provided as command argument")
        print("To run `main_gnss.py` please do:\n\t$ python main_gnss.py <path_to_config>")
        exit()

    # read config file
    try:
        with open(config_filename) as json_file:
            data = json.load(json_file)
        config_dict.init(data)
    except Exception as e:
        print(f"Error Reading Configuration File: {e}\nfilename = {config_filename}")
        exit()

    # create algorithm and algorithm manager
    try:
        alg = GnssSinglePointSolution()
        alg_mng = GnssAlgorithmManager(alg)

        # run algorithm
        alg_mng.run()

        # copy config json file to output dir
        shutil.copyfile(config_filename, f"{alg_mng.data_dir}\\config.json")
    except Exception as e:
        print(f"Unexpected error running while running algorithm: {e}")
        exit()


print("#--------------------------------------------------#")
print("#           Welcome to PyGNSSFix Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
