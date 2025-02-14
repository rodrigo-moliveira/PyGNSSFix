import numpy as np
from math import floor

from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger
from src.data_mng.gnss.phase_center_mng import PhaseCenterManager
from src.data_types.date import Epoch
from src.data_types.gnss import get_data_type, data_type_from_rinex, SatelliteAntenna, PhaseCenter
from src.data_types.gnss import get_satellite
from src.errors import FileError
from src.io.config import config_dict
from src.io.rinex_parser.utils import RINEX_SATELLITE_SYSTEM, GPS_ANTENNAS, GAL_ANTENNAS


# TODO: add documentation


class AntexReader:
    def __init__(self, file, phase_center: PhaseCenterManager):
        self.file = file
        self.phase_center = phase_center
        self._active_constellations = config_dict.get("model", "constellations", fallback=["GPS", "GAL"])
        self._user_services = dict()
        for constellation in self._active_constellations:
            service_list = [f"C{service}" for service in config_dict.get_services()[constellation]['user_service']]
            self._user_services[constellation] = [data_type_from_rinex(service, constellation) for service in
                                                  service_list]
        self._epoch = Epoch.strptime(config_dict.get("inputs", "arc", "first_epoch"), scale=str("GPS"))

        f_handler = open(f"{WORKSPACE_PATH}/{file}", "r")

        self.log = get_logger(IO_LOG)
        self.log.info(f"Reading Antex file {WORKSPACE_PATH}/{file}...")
        self.log.debug(f"Receiver Antenna selected: {self.phase_center.get_receiver_antenna().ant_type}")
        self.log.debug(f"Active satellite constellations: {self._active_constellations}")
        self.log.debug(f"User services selected: {self._user_services}")

        # read header
        self._read_header(f_handler)

        # read inputs
        self._read_data(f_handler)

        f_handler.close()

    def _read_header(self, f_handler):
        line = f_handler.readline()
        version = float(line[5:10])

        if "ANTEX VERSION / SYST" not in line:
            raise FileError("The provided Antex file {} is not valid")

        if version > 1.4:
            raise FileError("The provided Antex file {} is not valid. "
                            "Only bias Antex files of version up to 1.4 are currently supported. A file with version {}"
                            " was provided instead.".format(f_handler.name, version))

        while line:
            line = f_handler.readline()

            if "PCV TYPE / REFANT" in line:
                mode = line[0]
                if mode != 'A':
                    raise FileError("The provided bias Antex file {} is not valid. "
                                    "Currently only Absolute (A) PCV's are supported. PCV Type {} was provided.".
                                    format(f_handler.name, mode))

            elif "END OF HEADER" in line:
                break

    def _read_data(self, f_handler):
        receiver_antenna = self.phase_center.get_receiver_antenna()
        line = " "
        while line:
            if "START OF ANTENNA" in line:
                line = f_handler.readline()
                antenna_type = line[0:20]
                if antenna_type == receiver_antenna.ant_type:
                    # matched antenna to process
                    self.log.debug(f"Storing receiver antenna {antenna_type}")
                    self._read_antenna(f_handler, receiver_antenna)
                if antenna_type.strip() in GPS_ANTENNAS + GAL_ANTENNAS:
                    sat_str = line[20:23]
                    sat = get_satellite(sat_str)

                    sat_antenna = SatelliteAntenna(sat)
                    sat_antenna.ant_type = antenna_type
                    if self._read_antenna(f_handler, sat_antenna):
                        self.log.debug(f"Storing {sat.sat_system} satellite antenna {antenna_type} for sat {sat_str}")
                        self.phase_center.add_satellite_antenna(sat, sat_antenna)
            line = f_handler.readline()

    def _read_antenna(self, f_handler, antenna):
        f_handler.readline()  # METH / BY / # / DATE
        line = f_handler.readline()  # DAZI
        d_az = float(line.split()[0])  # d_az may be zero
        if d_az > 0:
            n_az = int(360 // d_az) + 1
            azimuth_vec = np.linspace(0, 360, n_az)
        else:
            azimuth_vec = np.array([])
        antenna.azimuth_vec = azimuth_vec

        line = f_handler.readline()  # ZEN1 / ZEN2 / DZEN
        tokens = line.split()[0:3]
        d_zenith = [float(x) for x in tokens]
        n_zenith = int((d_zenith[1] - d_zenith[0]) // d_zenith[2]) + 1
        antenna.zenith_vec = np.linspace(d_zenith[0], d_zenith[1], n_zenith)

        while line:
            line = f_handler.readline()
            if "START OF FREQUENCY" in line:
                freq_name = line.split()[0]
                freq_datatype = self.datatype_from_freq(freq_name)

                if freq_datatype is not None:
                    for service in self._user_services[freq_datatype.constellation]:
                        if service.freq == freq_datatype:
                            # read this frequency data
                            phase_center = self._read_frequency(f_handler, antenna)
                            antenna.set_freq_data(freq_datatype, phase_center)

            if "VALID FROM" in line:
                first_epoch = self.get_epoch(line)
                if first_epoch > self._epoch:
                    # ignore this antenna due to out of date
                    return False

            elif "VALID UNTIL" in line:
                last_epoch = self.get_epoch(line)
                if last_epoch < self._epoch:
                    # ignore this antenna due to out of date
                    return False

            elif "END OF ANTENNA" in line:
                break
        return True

    def _read_frequency(self, f_handler, antenna):
        line = " "
        phase_center = PhaseCenter()

        n_az = len(antenna.azimuth_vec)
        n_zenith = len(antenna.zenith_vec)
        if n_az > 0:
            pcv_matrix = np.zeros((n_az, n_zenith))
        else:
            pcv_matrix = None

        while line:
            line = f_handler.readline()
            if "END OF FREQUENCY" in line:
                phase_center.pcv_azi = pcv_matrix
                break

            elif "NORTH / EAST / UP" in line:
                tokens = line.split()[0:3]
                phase_center.pco = [float(x) for x in tokens]

            else:
                tokens = line.split()
                this_pcv = [float(x) for x in tokens[1:]]
                if tokens[0] == 'NOAZI':
                    phase_center.pcv_noazi = this_pcv
                else:
                    this_az = float(tokens[0])
                    az_idx = np.where(antenna.azimuth_vec == this_az)[0]
                    if len(az_idx) != 1:
                        self.log.warning(f"Error reading azimuth {this_az} for antenna {antenna.ant_type}. Skipping...")
                        continue
                    az_idx = az_idx[0]
                    pcv_matrix[az_idx, :] = this_pcv

        return phase_center

    def datatype_from_freq(self, freq):
        constellation = RINEX_SATELLITE_SYSTEM.get(freq[0], "UNKNOWN")
        if constellation in self._active_constellations:
            freq_no = freq[1:]
            if constellation == "GPS":
                if freq_no == "01":
                    return get_data_type("L1", constellation)
                elif freq_no == "02":
                    return get_data_type("L2", constellation)
                elif freq_no == "05":
                    return get_data_type("L5", constellation)
            elif constellation == "GAL":
                if freq_no == "01":
                    return get_data_type("E1", constellation)
                elif freq_no == "05":
                    return get_data_type("E5a", constellation)
                elif freq_no == "07":
                    return get_data_type("E5b", constellation)
                elif freq_no == "08":
                    return get_data_type("E5AltBOC", constellation)
                elif freq_no == "06":
                    return get_data_type("E6", constellation)
        return None

    @staticmethod
    def get_epoch(line):
        tokens = line.split()
        year = int(tokens[0])
        month = int(tokens[1])
        day = int(tokens[2])
        hour = int(tokens[3])
        minute = int(tokens[4])
        second = floor(float(tokens[5]))
        return Epoch(year, month, day, hour, minute, second, scale=str("GPS"))
