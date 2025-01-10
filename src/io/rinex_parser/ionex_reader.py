""" Parser of Ionospheric Global Map (GIM) files IONEX (version 1.0 or 1.1)
The official documentation of these files may be found in
http://ftp.aiub.unibe.ch/ionex/draft/ionex11.pdf for version v1.1
"""
from math import ceil
import numpy as np

from src import WORKSPACE_PATH
from src.common_log import get_logger, IO_LOG
from src.errors import FileError
from src.io.config import config_dict

from . import utils
from ...data_types.date import Epoch


class IONEXReader:
    """
    Parser of Rinex IONEX files
    """

    def __init__(self, file, gim):
        """
        Reads the provided RINEX IONEX file and stores its content in the provided `GlobalIonoMap` instance.

        Args:
            file(str): path to the input RINEX IONEX file
            gim(src.data_mng.gnss.global_iono_map.GlobalIonoMap): the `GlobalIonoMap` object to store the GIM maps
        """

        first_obs_str = config_dict.get("inputs", "arc", "first_epoch")
        last_obs_str = config_dict.get("inputs", "arc", "last_epoch")
        first_obs_epoch = Epoch.strptime(first_obs_str, scale="UTC")
        last_obs_epoch = Epoch.strptime(last_obs_str, scale="UTC")

        f_handler = open(f"{WORKSPACE_PATH}/{file}", "r")
        self.log = get_logger(IO_LOG)
        self.log.info(f"Reading IONEX file {WORKSPACE_PATH}/{file}...")

        # read header
        self._read_header(f_handler, gim)

        # perform validations
        if gim.header.first_epoch is not None and first_obs_epoch < gim.header.first_epoch:
            self.log.warning(f"The first observation epoch {first_obs_epoch} "
                             f"is earlier than the first epoch in the IONEX file {gim.header.first_epoch}.")
        if gim.header.last_epoch is not None and last_obs_epoch > gim.header.last_epoch:
            self.log.warning(f"The last observation epoch {last_obs_epoch} "
                             f"is later than the last epoch in the IONEX file {gim.header.last_epoch}.")

        # read body
        gim.init_data()
        self._read_data(f_handler, gim)

        f_handler.close()

    @classmethod
    def _read_header(cls, f_handler, gim):
        """
        Method to read header data

        Header lines to look for (others are ignored):
            * IONEX VERSION / TYPE
            * EPOCH OF FIRST MAP
            * EPOCH OF LAST MAP
            * # OF MAPS IN FILE
            * MAPPING FUNCTION
            * BASE RADIUS
            * MAP DIMENSION
            * HGT1 / HGT2 / DHGT
            * LAT1 / LAT2 / DLAT
            * LON1 / LON2 / DLON
            * EXPONENT
            * END OF HEADER
        """
        line = " "

        while line:
            line = f_handler.readline()

            if "IONEX VERSION / TYPE" in line:
                version = utils.to_float(line.split()[0])
                rinex_type = line.split()[1][0]
                if rinex_type != 'I':
                    raise FileError("The provided file should be a RINEX IONEX File. Instead, "
                                    "a {} was provided (code {})".format(utils.RINEX_FILE_TYPES.get
                                                                         (rinex_type, "Unknown Data File"), rinex_type))
                if version > 1.1:
                    raise FileError(f"The provided IONEX file is of version {version}. "
                                    f"Only versions 1.0 and 1.1 are supported.")

            elif "EPOCH OF FIRST MAP" in line:
                gim.header.first_epoch = IONEXReader._get_epoch(line)

            elif "EPOCH OF LAST MAP" in line:
                gim.header.last_epoch = IONEXReader._get_epoch(line)

            elif "# OF MAPS IN FILE" in line:
                gim.header.n_maps = int(line.split()[0])

            elif "MAPPING FUNCTION" in line:
                gim.header.mapping_function = line.split()[0]

            elif "BASE RADIUS" in line:
                gim.header.base_radius = utils.to_float(line.split()[0])

            elif "MAP DIMENSION" in line:
                map_dimension = int(line.split()[0])
                if map_dimension != 2:
                    raise FileError(f"Only 2D maps are supported. The provided map has {map_dimension} dimensions.")

            elif "HGT1 / HGT2 / DHGT" in line:
                tokens = line.split()
                gim.header.hgt = [utils.to_float(tokens[0]), utils.to_float(tokens[1]), utils.to_float(tokens[2])]

            elif "LAT1 / LAT2 / DLAT" in line:
                tokens = line.split()
                gim.header.lat = [utils.to_float(tokens[0]), utils.to_float(tokens[1]), utils.to_float(tokens[2])]

            elif "LON1 / LON2 / DLON" in line:
                tokens = line.split()
                gim.header.lon = [utils.to_float(tokens[0]), utils.to_float(tokens[1]), utils.to_float(tokens[2])]

            elif "EXPONENT" in line:
                gim.header.exponent = int(line.split()[0])

            elif "END OF HEADER" in line:
                break

    def _read_data(self, f_handler, gim):
        """ Method to read the body of the IONEX file

        Only TEC Maps are read. RMS and Height maps are ignored.
        """
        line = " "
        epoch = None
        ignoring = False
        tec_matrix = None

        while line:
            line = f_handler.readline()
            if "END OF RMS MAP" in line or "END OF HEIGHT MAP" in line:
                ignoring = False
                continue

            if ignoring:
                continue

            if "END OF TEC MAP" in line:
                gim.add_data(epoch, tec_matrix)
                epoch = None
                tec_matrix = None

            elif "START OF TEC MAP" in line:
                tec_matrix = np.zeros((len(gim.lat_array), len(gim.lon_array)))

            elif "START OF RMS MAP" in line or "START OF HEIGHT MAP" in line:
                ignoring = True

            elif "EPOCH OF CURRENT MAP" in line:
                epoch = IONEXReader._get_epoch(line)

            elif "LAT/LON1/LON2/DLON/H" in line:
                lat = utils.to_float(line[2:8])
                lon1 = utils.to_float(line[8:14])
                # lon2 = utils.to_float(line[14:20])
                dlon = utils.to_float(line[20:26])
                # h = utils.to_float(line[26:32])
                # lon_interval = abs(lon2-lon1)
                # n_entries = int(lon_interval/dlon)+1

                n_lines = ceil(len(gim.lon_array)/16)
                lat_idx = np.where(gim.lat_array == lat)[0]
                if len(lat_idx) != 1:
                    self.log.warning(f"Latitude {lat} not found in the latitude array. Ignoring this TEC map.")
                    continue
                lat_idx = lat_idx[0]

                this_lon = lon1
                for i in range(n_lines):
                    line = f_handler.readline()
                    for j in range(16):
                        if i*16+j >= len(gim.lon_array):
                            break

                        lon_idx = np.where(gim.lon_array == this_lon)[0]
                        if len(lon_idx) != 1:
                            self.log.warning(f"Longitude {this_lon} not found in the longitude array. "
                                             f"Ignoring this TEC map.")
                            break

                        lon_idx = lon_idx[0]

                        # update TEC Matrix entry
                        tec_matrix[lat_idx, lon_idx] = utils.to_float(line[j*5:j*5+5])
                        this_lon = this_lon + dlon

    @classmethod
    def _get_epoch(cls, line):
        tokens = line.split()
        year = int(tokens[0])
        month = int(tokens[1])
        day = int(tokens[2])
        hour = int(tokens[3])
        minute = int(tokens[4])
        second = int(tokens[5])
        return Epoch(year, month, day, hour, minute, second, scale="UTC")
