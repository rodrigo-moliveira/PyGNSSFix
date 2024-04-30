"""Parser of RINEX SP3-c and SP3-d files
The official documentation of these files may be found in http://epncb.eu/ftp/data/format/sp3d.pdf
(compatible with SP3-d version)
"""
import numpy as np
import datetime

from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger
from src.data_types.gnss import get_satellite
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
            first_epoch(str or None): first observation epoch
            last_epoch(str or None): last observation epoch

        Note that for interpolation purposes, this reader starts saving data from 5 epochs before the `first_epoch`
        until 5 epochs after the `last_epoch` (if possible).
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

        self._first_epoch = Epoch.strptime(first_epoch, scale=str(self._time_sys)) - datetime.timedelta(
            seconds=5 * self._interval) if first_epoch is not None else None
        self._last_epoch = Epoch.strptime(last_epoch, scale=str(self._time_sys)) + datetime.timedelta(
            seconds=5 * self._interval) if first_epoch is not None else None

        # read inputs
        self._read_data(f_handler, epoch)

        f_handler.close()

    def _read_header(self, file):
        """
        Method to read header data
        Values read:
            * version (check if it is SP3c or SP3d)
            * number of epochs in the file
            * time interval between epochs
            * time system for the reported epochs
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
        """Build the epoch object for the provided SP3 line"""
        tokens = line_str.split()

        year = int(tokens[1])
        month = int(tokens[2])
        day = int(tokens[3])
        hour = int(tokens[4])
        minute = int(tokens[5])
        second = int(float(tokens[6]))
        return Epoch(year, month, day, hour, minute, second, scale=str(self._time_sys))

    def _read_data(self, file, first_epoch):
        """
        Method to read the orbit data (only position is read, velocities and clocks are ignored)
        """
        line = " "
        curr_epoch = first_epoch

        while line:
            line = file.readline()

            # End of File detected
            if "EOF" in line:
                break

            # New epoch detected
            if line[0] == '*':
                curr_epoch = self.get_epoch(line)

            # check if this entry is inside the valid interval or not
            if self._first_epoch:
                if self._first_epoch > curr_epoch:
                    continue
            if self._last_epoch:
                if self._last_epoch < curr_epoch:
                    break  # after the last epoch has been reached, we can stop reading

            tokens = line.split()
            if len(tokens) != 0 and tokens[0][0] == "P":
                try:
                    # process new satellite orbit
                    sat = get_satellite(tokens[0][1:])

                    if sat.sat_system not in self._active_constellations:
                        continue

                    orbit = [utils.to_float(tokens[i]) * 1000.0 for i in [1, 2, 3]]  # get x,y,z components, in m (ECEF)
                    self.orbits.set_data(curr_epoch, sat, np.array(orbit))
                except Exception as e:
                    self.log.warn(f"Unable to process line: {line} due to error: {e}")
