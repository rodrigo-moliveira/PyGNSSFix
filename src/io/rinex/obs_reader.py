from math import floor

from src.data_types.date.date import Epoch
from src.data_types.gnss.satellite import get_satellite
from src.data_types.gnss.data_type import get_data_type
from src.data_types.gnss.observation_data import ObservationData, ObservationHeader
from src.data_types.gnss.service_utils import Services
from src.errors import ConfigError, FileError
from src import WORKSPACE_PATH

from .utils import *


class RinexObsReader:
    """
    Class RinexObsReader

    Attributes
        ----------
        obs: ObservationData
        header : ObservationHeader

    """

    def __init__(self, file: str, services, log=None,
                 first_arc_epoch=None, last_arc_epoch=None, snr_control_check=0):

        # instance variables
        self._file = file
        self._services = services
        self._map = {}
        self.obs = ObservationData()
        self.header = ObservationHeader()
        self._log = log
        self._snr_control_check = snr_control_check

        # first and final arc epochs to read, if not None
        self._first_arc_epoch = first_arc_epoch
        self._last_arc_epoch = last_arc_epoch

        f_handler = open(f"{WORKSPACE_PATH}/{self._file}", "r")

        # read header
        self._read_header(f_handler)

        self._validate_requested_observations()

        # read observation data
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
            print(line)

            if "RINEX VERSION / TYPE" in line:
                self.header.rinex_version = float(line[5:10])
                if self.header.rinex_version < 3:
                    raise FileError(f"The provided rinex file {self._file} is of version {self.header.rinex_version}. "
                                    f"Only version 3.00 or higher is supported.")

                rinex_type = line[20]
                if rinex_type != 'O':
                    raise FileError(f"Rinex File {self._file} should be a GNSS Observation Data File. Instead, a "
                                    f"{RINEX_FILE_TYPES.get(rinex_type, 'Unknown Data File')} file was provided")

                self.header.satellite_system = line[40]

            elif "APPROX POSITION XYZ" in line:
                data = line[:RINEX_OBS_END_OF_DATA_HEADER]
                # static antenna position in ECEF frame (static receiver)
                self.header.receiver_position = tuple(map(float, data.split()))

            elif "LEAP SECONDS" in line:
                data = line[:RINEX_OBS_END_OF_DATA_HEADER].split()
                self.header.leap_seconds = int(data[0])

            elif "TIME OF FIRST OBS" in line:
                data = line[:RINEX_OBS_END_OF_DATA_HEADER]
                self._set_time_system(data[48:51].upper())

                data = data.split()
                self.header.first_epoch = Epoch(int(data[0]), int(data[1]), int(data[2]),
                                                int(data[3]), int(data[4]), floor(float(data[5])),
                                                scale='GPS')  # TODO change scale for the appropriate one

            elif "TIME OF LAST OBS" in line:
                data = line[:RINEX_OBS_END_OF_DATA_HEADER].split()
                self.header.last_epoch = Epoch(int(data[0]), int(data[1]), int(data[2]),
                                                int(data[3]), int(data[4]), floor(float(data[5])),
                                               scale='GPS')  # TODO change scale for the appropriate one

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

                        # iterate over all available observation codes for this constellation
                        for this_obsCode in data[2:]:
                            if len(this_obsCode) == 3:
                                this_service = this_obsCode[1:]
                                if this_service in services:
                                    this_type = this_obsCode[0]
                                    if this_type in RINEX_OBS_TYPES_TO_READ:
                                        _map[this_obsCode] = data[2:].index(this_obsCode)
                        self._map[constellation] = _map

            elif "END OF HEADER" in line:
                # TODO: este mapa não está certo..
                print(self._map)
                exit()
                break

    def _validate_requested_observations(self):
        for constellation, services in self._services.items():
            if constellation in self._map:
                services_read = [x[1:] for x in self._map[constellation]]
                services_read = list(dict.fromkeys(services_read))

                for service in services:
                    if service not in services_read:
                        raise ConfigError(f"User-selected service "
                                          f"'{service}' does not exist in provided Observation File")

    def _set_time_system(self, time_system: str):

        # set default time system (according to the provided satellite system)
        if not time_system or time_system.isspace():
            self.header.time_system = RINEX_SATELLITE_SYSTEM.get(self.header.satellite_system, "GPS")
        else:
            self.header.time_system = time_system

        # NOTE: currently, this software only interprets GPS time system!!!
        if self.header.time_system != "GPS":
            self._log.warn(f"Time system `{self.header.time_system}´ will be interpreted internally as GPS")
            self.header.time_system = "GPS"

    def _read_obs(self, cFile):
        """
        Read observation data

        Each observation is stored in a fixed field of 14 characters and
        with three trailing digits after the decimal point. Carrier phase
        and pseudorange observations are furthermore complemented by an
        optional loss-of-lock indicator (0, blank, or 1) and/or a
        single-digit signal-strength indicator in the two cells adjacent
        to the actual measurement value. Although the RINEX format supports
        the loss-of lock indicator for each observation type, it is common
        practice to only indicate it on the phase observations. Likewise,
        the single-digit signal-strength field is often omitted if the
        signal strength (in dB-Hz) is explicitly provided as an observation
        type.

        Epoch (example):
            > 2019 01 14 06 15  0.0000000  0 38
        if the epoch flag is different from zero, this epoch is ignored

        Data line:
        G01  21169746.382   111247790.951       -2243.909          50.312    21169747.722    86686606.781       -1748.49

        NOTE: if signal strength bit is < 5, then the observation is discarded
        """

        line = " "
        this_epoch = None
        ignoring = False

        """
        while line:
            line = cFile.readline()

            if not line.strip():
                # look for empty lines
                break

            if line[0] == ">":
                # Reading new epoch
                ignoring = False
                data = line[1:].split()
                time_dict = {"year": int(data[0]),
                             "month": int(data[1]),
                             "day": int(data[2]),
                             "hour": int(data[3]),
                             "minute": int(data[4]),
                             "second": floor(float(data[5]))}
                epochFlag = int(data[6])
                # nmbOfSats = data[7]

                if epochFlag != 0:
                    ignoring = True
                    self.log.debug(f"Discarding all data for {this_epoch.to_time_stamp()} due to bad epoch flag")
                    continue

                this_epoch = Epoch(time_dict, time_system=self.header.time_system)

                # check initial and final arc intervals
                if self.first_arc_epoch:
                    if self.first_arc_epoch > this_epoch:
                        ignoring = True
                        continue
                if self.last_arc_epoch:
                    if this_epoch > self.last_arc_epoch:
                        ignoring = True
                        continue

            else:
                if ignoring:
                    continue
                # Reading observations
                # Each observation word is 16 characters: 14 (observation) + 1 (loss-of-lock indicator) + 1
                # (signal strength).
                # the last two are optional -> may be blank. Here they will be discarded.
                constellation_code = line[0]

                constellation = RinexUtils.RINEX_SATELLITE_SYSTEM[constellation_code]
                if constellation in self._map:

                    # get satellite
                    this_sat = SatelliteFactory(line[0:3])
                    this_map = self._map[constellation]

                    # get available observations in the current line of the file
                    obs_str = line[3:]
                    observations = [obs_str[i:i + 16] for i in range(0, len(obs_str), 16)]
                    for this_obsCode, this_index in this_map.items():
                        this_type = DataTypeFactory(this_obsCode[0:2])

                        # get observable
                        try:
                            this_str = observations[this_index]
                            if this_str.isspace():
                                continue

                            observable = float(this_str[0:-2])
                            signal_strength = this_str[-1]

                            try:
                                signal_strength = int(signal_strength)
                                if signal_strength < self.snr_control_check:
                                    self.log.debug(f"Discarding observable {this_type} at {this_epoch.to_time_stamp()} "
                                                   f"for {this_sat} due to low signal strength: {signal_strength} < "
                                                   f"{self.snr_control_check}")
                                    continue
                            except ValueError:
                                pass

                            # set observable
                            self.cObsData.set_observable(this_epoch, this_sat, this_type, observable)
                            # print("Setting observable", this_epoch.to_time_stamp(), this_sat, this_type, observable)

                        except (ValueError, IndexError) as e:
                            self.log.debug("problem parsing line {} for index {}: {}".format(line, this_index, e))
        """