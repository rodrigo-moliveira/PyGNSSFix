import sys
import json
import shutil

from src.modules.gnss import GnssAlgorithmManager
from src.io.config import config_dict
from src.common_log import set_logs


def main():
    # get json config file from command line args
    try:
        config_filename = sys.argv[1]
    except IndexError:
        print("ERROR: No configuration json file was provided as command argument")
        print("To run `main_gnss.py` please do:\n\t$ python main_gnss.py <path_to_config>")
        exit(-1)

    # read config file
    try:
        with open(config_filename) as json_file:
            data = json.load(json_file)
        config_dict.init(data, alg="GNSS")
    except Exception as e:
        print(f"Error Reading Configuration File: {e}\nfilename = {config_filename}")
        exit()

    # create algorithm and algorithm manager
    try:
        alg_mng = GnssAlgorithmManager()

        # initialize logger objects
        set_logs(config_dict.get("log", "minimum_level"), f"{alg_mng.data_dir}\\log.txt")

        # run algorithm
        alg_mng.run()

        # copy config json file to output dir
        shutil.copyfile(config_filename, f"{alg_mng.data_dir}\\config.json")
    except Exception as e:
        print(f"Unexpected error running while running program: {e}")
        exit()


print("#--------------------------------------------------#")
print("#           Welcome to GNSSNavPy Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
