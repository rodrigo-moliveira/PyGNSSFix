import sys
import json
import shutil

from src.algorithms.gnss_alg_manager import GnssAlgorithmManager
from src.io.config.gnss_config import ConfigGNSS
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
        config = ConfigGNSS(**data)
    except (ConfigError, ConfigTypeError) as e:
        print(f"[ERROR READING CONFIGURATION FILE]: {e}\nfilename = {config_filename}")
        exit()

    # create algorithm and algorithm manager
    alg = GnssSinglePointSolution()
    alg_mng = GnssAlgorithmManager(alg, config)

    # run algorithm
    alg_mng.run()

    # copy config json file to output dir
    shutil.copyfile(config_filename, f"{alg_mng.data_dir}\\config.json")


print("#--------------------------------------------------#")
print("#           Welcome to PyGNSSFix Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
