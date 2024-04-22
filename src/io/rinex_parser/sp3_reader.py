"""Parser of RINEX SP3-d files
The official documentation of these files may be found in http://epncb.eu/ftp/data/format/sp3d.pdf
(compatible with SP3-d version)
"""
import numpy as np

from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger
from src.data_types.gnss import get_satellite, get_constellation
from src.data_types.date import Epoch
from src.errors import FileError
from src.io.config.config import config_dict

from . import utils


class SP3OrbitReader:
    """
    Parser of SP3 files
    """

    def __init__(self, file, orbits, first_epoch=None, last_epoch=None):
        """
        Args:
            file(str): path to the input RINEX Clock file
            orbits(src.data_mng.gnss.sat_orbit_data.SatelliteOrbits): the `SatelliteOrbits` object to store the
                satellite orbits (precise)
            first_epoch(str or None): First epoch to read (ignore if None)
            last_epoch(str or None): Last epoch to read (ignore if None)
        """

        # instance variables
        self.orbits = orbits
        self._time_sys = None
        self._n_of_epochs = 0
        self._interval = 0
        self._active_constellations = config_dict.get("model", "constellations", fallback=["GPS", "GAL"])

        f_handler = open(f"{WORKSPACE_PATH}/{file}", "r")
        self.log = get_logger(IO_LOG)
        self.log.info(f"Reading SP3 file {WORKSPACE_PATH}/{file}...")

        # read header
        epoch = self._read_header(f_handler)

        #self._first_epoch = Epoch.strptime(first_epoch, scale=str(self._time_sys)) if first_epoch is not None else None
        #self._last_epoch = Epoch.strptime(last_epoch, scale=str(self._time_sys)) if last_epoch is not None else None

        # read inputs
        #self._read_data(f_handler)

        #f_handler.close()
        exit()

    def _read_header(self, file):
        """
        Method to read header data

        TODO: to update
        Header lines to look for (others are ignored):
            * RINEX VERSION / TYPE  -> rinex_version and satellite_system
            * TIME SYSTEM ID        -> Time system used for time tags
            * # / TYPES OF DATA     -> List of clock data types (required AS)
            * PRN LIST              -> List of all satellites reported in this file
            * END OF HEADER         -> end of header section
        """
        # read first line
        line = file.readline()
        version = line[1]
        self._n_of_epochs = int(line[32:39])
        if version.lower() not in ["d", "c"]:
            raise FileError(f"The provided SP3 file is not a valid SP3c or SP3d file. Version is {version}")

        # read second line
        line = file.readline()
        self._interval = utils.to_float(line[24:38])

        # discard lines starting with + and ++
        line = file.readline()
        while line[0] == '+':
            line = file.readline()

        # read first %c line
        self._time_sys = line[9:12]

        # ignore the rest of the header lines
        line = file.readline()
        while line[0] != '*':
            line = file.readline()

        return self.get_epoch(line)

    def get_epoch(self, line_str):
        tokens = line_str.split()

        year = int(tokens[1])
        month = int(tokens[2])
        day = int(tokens[3])
        hour = int(tokens[4])
        minute = int(tokens[5])
        second = int(float(tokens[6]))
        return Epoch(year, month, day, hour, minute, second, scale=str(self._time_sys))

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

                    if self._prev_epoch[0] == [year, month, day, hour, minute, second]:
                        date = self._prev_epoch[1]
                    else:
                        date = Epoch(year, month, day, hour, minute, second, scale=str(self._time_sys))
                        self._prev_epoch = ([year, month, day, hour, minute, second], date)

                    # check if this entry is inside the valid interval or not
                    if self._first_epoch:
                        if self._first_epoch > date:
                            continue
                    if self._last_epoch:
                        if self._last_epoch < date:
                            continue

                    clock_bias = utils.to_float(tokens[9])  # clock bias in seconds

                    self.clocks.set_data(date, sat, clock_bias)
                except Exception as e:
                    self.log.warn(f"Unable to process line: {line} due to error: {e}")
