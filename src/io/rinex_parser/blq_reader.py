""" Parser of Ocean Loading Parameters in BLQ Format.
Although I could not find any official documentation for these kind of products, an example is provided here:
http://holt.oso.chalmers.se/loading/example_blq.html
"""

from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger
from src.errors import FileError


class BLQReader:
    """
    Parser of Ocean Loading Parameter files in BLQ Format.

    The following elements are parsed from these files:
        * Radial Amplitudes: 11 parameters, in meters
        * West Amplitudes: 11 parameters, in meters
        * South Amplitudes: 11 parameters, in meters
        * Radial Phases: 11 parameters, in degrees
        * West Phases: 11 parameters, in degrees
        * South Phases: 11 parameters, in degrees


    Limitations:
        * This parser currently only allows to read data for one station. If a file with multiple stations is provided
        an error will be raised.
    """

    def __init__(self, file, ocean_loading_manager, station):
        """
        Reads the provided Ocean Loading file and stores its content in the provided `ocean_loading_manager` instance.

        Args:
            file(str): path to the input Ocean Loading file
            ocean_loading_manager(src.data_mng.gnss.ocean_loading_data.OceanLoadingData): the `OceanLoadingData` object
                to store the parameters read.
            station(str): The station to be read

        Raises:
            FileError: an exception is raised if the data is not properly parsed
        """
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
                line_control = self._process_line(line, line_control, station)

                if line_control == 7:
                    break

    def _process_line(self, line, line_control, station):
        tokens = line.split()

        if line_control == 0:
            # read station
            line_control = 1
            if station != tokens[0]:
                raise FileError(f"The station in the ocean loading file ({tokens[0]}) does not match the station "
                                f"defined in the user configuration: {station}")
        else:
            data_array = self._read_data_line(line, tokens)
            if line_control == 1:
                # read Radial Amplitude parameters
                self.ocean_loading_manager.radial_amplitude = data_array

            elif line_control == 2:
                # read West Amplitude parameters
                self.ocean_loading_manager.west_amplitude = data_array

            elif line_control == 3:
                # read South Amplitude parameters
                self.ocean_loading_manager.south_amplitude = data_array

            elif line_control == 4:
                # read Radial Phase parameters
                self.ocean_loading_manager.radial_phase = data_array

            elif line_control == 5:
                # read West Phase parameters
                self.ocean_loading_manager.west_phase = data_array

            elif line_control == 6:
                # read South Phase parameters
                self.ocean_loading_manager.south_phase = data_array

            line_control = line_control + 1
        return line_control

    @staticmethod
    def _read_data_line(line, tokens):
        if len(tokens) != 11:
            raise FileError(f"The provided data line does not have 11 elements for reading:\n{line}")
        return [float(x) for x in tokens]
