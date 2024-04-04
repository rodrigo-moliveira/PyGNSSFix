"""Parser of RINEX OBS files (version 3.XX)
The official documentation of these files may be found in https://igs.org/wg/rinex/#documents-formats
"""

from math import floor

from src.data_types.date import Epoch
from src.data_types.gnss import get_satellite, data_type_from_rinex, DataType, get_constellation, service_utils
from src.data_mng.gnss.observation_data import ObservationData
from src.errors import ConfigError, FileError
from src.constants import SPEED_OF_LIGHT
from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger

from .utils import *


class RinexObsReader:
    """
    Parser of Rinex OBS files
    """

    def __init__(self, obs: ObservationData, file: str, services,
                 first_arc_epoch=None, last_arc_epoch=None, snr_control_check=0):

        # instance variables
        self._services = services
        self._map = {}
        self._obs = obs
        self._snr_control_check = snr_control_check
        self._log = get_logger(IO_LOG)

        f_handler = open(f"{WORKSPACE_PATH}/{file}", "r")
        self._log.info(f"Reading observation file {WORKSPACE_PATH}/{file}...")

        # read header
        self._read_header(f_handler)

        # first and final arc epochs to read, if not None
        self._first_arc_epoch = Epoch.strptime(first_arc_epoch, scale=self.time_system)
        self._last_arc_epoch = Epoch.strptime(last_arc_epoch, scale=self.time_system)

        self._validate_requested_observations()

        # read gnss_models data
        self._read_obs(f_handler)

        f_handler.close()

    def _read_header(self, file):
        """
        Method to read header data and store it in container ObservationHeader object

        Tags to look for:
            * RINEX VERSION / TYPE  -> rinex_version and satellite_system
            * APPROX POSITION XYZ   -> receiver_position
            * LEAP SECONDS          -> leap_seconds
            * TIME OF FIRST OBS     -> first_epoch
            * TIME OF LAST OBS      -> last_epoch
            * SYS / # / OBS TYPES   -> fill attribute self._map (map between data types and columns in file)
            * END OF HEADER         -> end of header section
        """
        line = " "

        while line:
            line = file.readline()

            if "RINEX VERSION / TYPE" in line:
                rinex_version = float(line[5:10])
                if rinex_version < 3 or rinex_version >= 4:
                    raise FileError(f"Provided rinex_parser file {file} is of version {rinex_version}"
                                    f". Only versions 3.XX are supported.")
                self.satellite_system = line[40]

                rinex_type = line[20]
                if rinex_type != 'O':
                    raise FileError(f"Rinex File {file} should be a GNSS Observation Data File. Instead, a "
                                    f"{RINEX_FILE_TYPES.get(rinex_type, 'Unknown Data File')} file was provided")

            # elif "APPROX POSITION XYZ" in line:
            #    data = line[:RINEX_OBS_END_OF_DATA_HEADER]
            #    # static antenna position in ECEF frame (static receiver)
            #    self._obs.header.receiver_position = tuple(map(float, data.split()))

            # elif "LEAP SECONDS" in line:
            #    data = line[:RINEX_OBS_END_OF_DATA_HEADER].split()
            #    self._obs.header.leap_seconds = int(data[0])

            elif "TIME OF FIRST OBS" in line:
                data = line[:RINEX_OBS_END_OF_DATA_HEADER]
                self._set_time_system(data[48:51].upper())

                data = data.split()
                self.first_epoch = Epoch(int(data[0]), int(data[1]), int(data[2]),
                                         int(data[3]), int(data[4]), floor(float(data[5])),
                                         scale=self.time_system)

            elif "TIME OF LAST OBS" in line:
                data = line[:RINEX_OBS_END_OF_DATA_HEADER].split()
                self.last_epoch = Epoch(int(data[0]), int(data[1]), int(data[2]),
                                        int(data[3]), int(data[4]), floor(float(data[5])),
                                        scale=self.time_system)

            elif "SYS / # / OBS TYPES" in line:
                data = line[:RINEX_OBS_END_OF_DATA_HEADER].split()
                const_code = line[0]
                nmb_of_obs = int(data[1])

                if nmb_of_obs > len(data[2:]):
                    line2 = file.readline()[:RINEX_OBS_END_OF_DATA_HEADER].split()
                    data += line2

                # iterate over all constellations
                for constellation, services in self._services.items():
                    if RINEX_SATELLITE_SYSTEM.get(const_code, "") == constellation:
                        _map = {}

                        # iterate over all available gnss_models codes for this constellation
                        for this_obsCode in data[2:]:
                            if len(this_obsCode) == 3:
                                this_service = this_obsCode[1:]
                                if this_service in services:
                                    this_type = this_obsCode[0]
                                    if this_type in RINEX_OBS_TYPES_TO_READ:
                                        _map[this_obsCode] = data[2:].index(this_obsCode)
                        self._map[constellation] = _map

            elif "END OF HEADER" in line:
                break

    def _validate_requested_observations(self):
        for constellation, services in self._services.items():
            if constellation in self._map:
                services_read = set([x[1:] for x in self._map[constellation]])
                for service in services:
                    if service not in services_read:
                        raise ConfigError(f"User-selected service "
                                          f"'{service}' {constellation} does not exist in provided Observation File")

    def _set_time_system(self, time_system: str):

        # set default time system (according to the provided satellite system)
        if not time_system or time_system.isspace():
            if self.satellite_system == 'M':
                raise FileError(f"No valid time system provided in header line TIME OF FIRST OBS")
            self.time_system = RINEX_SATELLITE_SYSTEM.get(self.satellite_system, "GPS")
        else:
            self.time_system = time_system

    def _read_obs(self, file):
        """
        Read gnss_models data

        Each gnss_models is stored in a fixed field of 14 characters and
        with three trailing digits after the decimal point. Carrier phase
        and pseudorange observations are furthermore complemented by an
        optional loss-of-lock indicator (0, blank, or 1) and/or a
        single-digit signal-strength indicator in the two cells adjacent
        to the actual measurement value. Although the RINEX format supports
        the loss-of lock indicator for each gnss_models type, it is common
        practice to only indicate it on the phase observations. Likewise,
        the single-digit signal-strength field is often omitted if the
        signal strength (in dB-Hz) is explicitly provided as an gnss_models
        type.

        Epoch (example):
            > 2019 01 14 06 15  0.0000000  0 38
        if the epoch flag is different from zero, this epoch is ignored

        Data line:
        G01  21169746.382   111247790.951       -2243.909          50.312    21169747.722    86686606.781       -1748.49

        NOTE: if signal strength bit is < 5, then the gnss_models is discarded
        """

        line = " "
        this_epoch = None
        ignoring = False

        while line:
            line = file.readline()

            if not line.strip():
                # look for empty lines
                break

            # restart epoch
            if line[0] == ">":
                # Reading new epoch
                ignoring = False
                data = line[1:].split()
                this_epoch = Epoch(int(data[0]), int(data[1]), int(data[2]),
                                   int(data[3]), int(data[4]), floor(float(data[5])),
                                   scale=self.time_system)

                epochFlag = int(data[6])
                # nmbOfSats = data[7]

                if epochFlag != 0:
                    ignoring = True
                    self._log.debug(f"Discarding all data for {this_epoch} due to bad epoch flag {epochFlag}")
                    continue

                # check initial and final arc intervals
                if self._first_arc_epoch:
                    if self._first_arc_epoch > this_epoch:
                        ignoring = True
                        continue
                if self._last_arc_epoch:
                    if this_epoch > self._last_arc_epoch:
                        ignoring = True
                        continue

            # reading gnss_models-data line for the current epoch
            else:
                if ignoring:
                    continue
                # Reading observations
                # Each gnss_models word is 16 characters: 14 (gnss_models) + 1 (loss-of-lock indicator) + 1
                # (signal strength).
                # the last two are optional -> may be blank. Here they will be discarded.

                try:
                    constellation = get_constellation(service_utils.CodeToConstellationMap[line[0]])
                except (ConfigError, KeyError):
                    continue

                if constellation in self._map:
                    # get satellite
                    this_sat = get_satellite(line[0:3])
                    this_map = self._map[constellation]

                    # get available observations in the current line of the file
                    obs_str = line[3:]
                    observations = [obs_str[i:i + 16] for i in range(0, len(obs_str), 16)]
                    for this_obsCode, this_index in this_map.items():
                        this_type = data_type_from_rinex(this_obsCode[0:2], constellation)

                        # get observable
                        try:
                            this_str = observations[this_index]
                            if this_str.isspace():
                                continue

                            observable = float(this_str[0:-2])
                            signal_strength = this_str[-1]

                            try:
                                signal_strength = int(signal_strength)
                                if signal_strength < self._snr_control_check:
                                    self._log.debug(f"Discarding observable {this_type} at {this_epoch} "
                                                    f"for {this_sat} due to low signal strength: {signal_strength} < "
                                                    f"{self._snr_control_check}")
                                    continue
                            except ValueError:
                                pass

                            # convert carrier phase from cycles to meters
                            if DataType.is_carrier(this_type):
                                wavelength = SPEED_OF_LIGHT / this_type.freq_value
                                observable = observable * wavelength

                            # set observable
                            self._obs.set_observable(this_epoch, this_sat, this_type, observable)
                            # print("Setting observable", this_epoch, constellation, this_sat, this_type, observable)

                        except (ValueError, IndexError):
                            pass
                            # self._log.debug("problem parsing line {} for index {}: {}".format(line, this_index, e))
