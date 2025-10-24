""" Parser of Ocean Loading Parameters in BLQ Format.
Although I could not find any official documentation for these kind of products, an example is provided here:
http://holt.oso.chalmers.se/loading/example_blq.html
"""

from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger
from src.errors import FileError


class BLQReader:
    """
    TO BE WRITTEN
    """

    def __init__(self, file, ocean_loading_manager, station):
        self.ocean_loading_manager = ocean_loading_manager
        f_handler = open(f"{WORKSPACE_PATH}/{file}", "r")

        self.log = get_logger(IO_LOG)
        self.log.info(f"Reading Ocean Loading Parameters file {WORKSPACE_PATH}/{file}...")
        self.log.info(f"Selected Station for Ocean Loading: {station}")

        self._read_data(f_handler, station)

        f_handler.close()

    def _read_data(self, file, station):
        # NOTE: currently the logic of this function is designed to only accept one station entry.
        # Undefined/untested behaviour if several station entries are provided
        line = " "
        start = False
        line_control = 0
        while line:
            line = file.readline()

            # ignore everything until $$ END HEADER
            if line.strip() == "$$ END HEADER":
                start = True
                continue

            if line.strip() == "$$ END TABLE":
                # end of file
                break

            if len(line) == 0 or (len(line) > 2 and line[0:2] == "$$"):
                # ignore comments
                continue

            if start:
                print(line[0:-1])
                line_control = self.process_line(line, line_control, station)

                if line_control == 7:
                    print("here")
                    break


        exit()

    def process_line(self, line, line_control, station):
        tokens = line.split()

        if line_control == 0:
            # read station
            line_control = 1
            if station != tokens[0]:
                raise FileError(f"The station in the ocean loading file ({tokens[0]}) does not match the station "
                                f"defined in the user configuration: {station}")

        elif line_control == 1:
            # read Radial Amplitude parameters
            line_control = 2

        elif line_control == 2:
            # read West Amplitude parameters
            line_control = 3

        elif line_control == 3:
            # read South Amplitude parameters
            line_control = 4

        elif line_control == 4:
            # read Radial Phase parameters
            line_control = 5

        elif line_control == 5:
            # read West Phase parameters
            line_control = 6

        elif line_control == 6:
            # read South Phase parameters
            line_control = 7

        return line_control
