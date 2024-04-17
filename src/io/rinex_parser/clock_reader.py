"""Parser of RINEX Clock files (version 3.XX)
The official documentation of these files may be found in https://files.igs.org/pub/data/format/rinex_clock304.txt
for v3.04
"""

from src.data_mng.gnss.sat_clock_data import SatelliteClocks
from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger
from src.data_types.gnss import get_satellite, get_constellation
from src.data_types.date import Epoch
from src.errors import FileError
from src.io.config.config import config_dict

from . import utils


class RinexClockReader:
    """
    Parser of Rinex Clock files
    """

    def __init__(self, file, clocks: SatelliteClocks, first_obs_epoch, last_obs_epoch):
        """
        Args:
            file(str): path to the input RINEX Clock file
            clocks(SatelliteClocks): the `SatelliteClocks` object to store the satellite clocks (precise)
            first_obs_epoch(str): First observation epoch
            last_obs_epoch(str): Last observation epoch
        """

        # instance variables
        self.clocks = clocks
        self._time_sys = None
        self._prn_list = []
        self._last_epoch = (None, None)
        self._active_constellations = config_dict.get("model", "constellations")

        f_handler = open(f"{WORKSPACE_PATH}/{file}", "r")
        self.log = get_logger(IO_LOG)
        self.log.info(f"Reading clock file {WORKSPACE_PATH}/{file}...")

        # read header
        self._read_header(f_handler)

        # TODO: self._first_epoch = Epoch()

        # read inputs
        self._read_data(f_handler)

        f_handler.close()

        test_epoch = Epoch(2019, 1, 14, 6, 0, 5, scale=str(self._time_sys))
        test_sat = get_satellite("E01")
        clock = self.clocks.get_clock(test_sat, test_epoch)
        print(test_epoch, test_sat, clock)
        exit()

    def _read_header(self, file):
        """
        Method to read header data

        Header lines to look for (others are ignored):
            * RINEX VERSION / TYPE  -> rinex_version and satellite_system
            * TIME SYSTEM ID        -> Time system used for time tags
            * # / TYPES OF DATA     -> List of clock data types (required AS)
            * PRN LIST              -> List of all satellites reported in this file
            * END OF HEADER         -> end of header section
        """
        line = " "

        while line:
            line = file.readline()

            if "RINEX VERSION / TYPE" in line:
                rinex_version = utils.to_float(line[0:4])
                if rinex_version < 3 or rinex_version >= 4:
                    raise FileError("The provided rinex_parser file {} is of version {}. Only versions 3.XX "
                                    "are supported. Error!".format(file, rinex_version))

                rinex_type = line[21]
                if rinex_type != 'C':
                    raise FileError("Rinex File {} should be a RINEX Clock File. Instead, "
                                    "a {} was provided (code {})".format(file,
                                                                         utils.RINEX_FILE_TYPES.get
                                                                         (rinex_type, "Unknown Data File"), rinex_type))

            elif "TIME SYSTEM ID" in line:
                tokens = line.split()
                self._time_sys = get_constellation(tokens[0])

            elif "# / TYPES OF DATA" in line:
                tokens = line.split()
                if "AS" not in tokens:
                    raise FileError(f"Missing # / TYPES OF DATA 'AS'")

            elif "PRN LIST" in line:
                tokens = line[0:65].split()
                self._prn_list += list(tokens)

            elif "END OF HEADER" in line:
                break

    def _read_data(self, file):
        line = " "

        while line:
            line = file.readline()

            tokens = line.split()
            if len(tokens) != 0 and tokens[0] == "AS":
                try:
                    # process new satellite clock
                    sat = get_satellite(tokens[1])

                    if sat.sat_system not in self._active_constellations:
                        continue

                    year = int(tokens[2])
                    month = int(tokens[3])
                    day = int(tokens[4])
                    hour = int(tokens[5])
                    minute = int(tokens[6])
                    second = int(float(tokens[7]))  # note: milliseconds are ignored

                    if self._last_epoch[0] == [year, month, day, hour, minute, second]:
                        date = self._last_epoch[1]
                    else:
                        date = Epoch(year, month, day, hour, minute, second, scale=str(self._time_sys))
                        self._last_epoch = ([year, month, day, hour, minute, second], date)

                    clock_bias = utils.to_float(tokens[9])  # clock bias in seconds

                    self.clocks.set_data(date, sat, clock_bias)
                except Exception as e:
                    self.log.warn(f"Unable to process line: {line} due to error: {e}")
