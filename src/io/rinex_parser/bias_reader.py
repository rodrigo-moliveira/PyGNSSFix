""" Parser of SINEX BIAS—Solution (Software/technique) INdependent EXchange (SINEX) Format for GNSS Biases Version 1.00
The official documentation of these products may be found in
https://files.igs.org/pub/data/format/sinex_bias_100.pdf
"""

from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger
from src.data_types.gnss import get_satellite
from src.data_types.date import Epoch
from src.errors import FileError
from src.io.config import config_dict, EnumSatelliteBias

from . import utils


class SinexBiasReader:
    """
    Parser of Bias-SINEX files for GNSS Code and Phase Biases.
    Station biases are ignored.

    From the official documentation, the available bias types in SINEX files are:
        * DSB: Differential Signal Biases, also known as DCB.
        * ISB: Ionosphere-Free Signal Biases
        * OSB: Observable-Specific Signal Biases

    In the context of this software, only DSB and OSB products are considered (ISBs will be ignored).

    Each Bias file may only contain one of the bias products detailed (either OSBs or DSBs).
    The bias product is detailed in the ``BIAS_MODE`` field of the ``BIAS/DESCRIPTION`` Block (header).

    The ``BIAS_MODE`` describes how the included GNSS bias values have to be interpreted and applied, respectively.
    Possible modes are:
        * ``RELATIVE``
        * ``ABSOLUTE``

    Obviously, this implies that inclusion of either relative (differential) or absolute (pseudo-absolute)
    GNSS bias values is allowed in a SINEX_BIAS file.
    The Relative mode refers to DSBs and the Absolute mode refers to OSBs.

    Note: The bias mode is part of the SINEX_BIAS header line (encoded with ``R`` or ``A``).
    """

    def __init__(self, file, bias_manager, bias_enum):
        """
        Reads the provided Bias SINEX file and stores its content in the provided `bias_manager` instance.

        Args:
            file(str): path to the input RINEX Clock file
            bias_manager(src.data_mng.gnss.bias_manager.BiasManager): the `BiasManager` object to store the
                satellite biases (precise)
            bias_enum(EnumSatelliteBias): enumeration with the bias mode selected to be read (OSBs or DCBs)
        """
        first_epoch = config_dict.get("inputs", "arc", "first_epoch")
        # last_epoch = config_dict.get("inputs", "arc", "last_epoch")

        # instance variables
        self.bias_manager = bias_manager
        self._bias_mode = None
        self._time_sys = None
        self._active_constellations = config_dict.get("model", "constellations", fallback=["GPS", "GAL"])

        f_handler = open(f"{WORKSPACE_PATH}/{file}", "r")

        self.log = get_logger(IO_LOG)
        self.log.info(f"Reading Bias SINEX file {WORKSPACE_PATH}/{file}...")

        # read header
        self._read_header(f_handler, bias_enum)

        self._epoch = Epoch.strptime(first_epoch, scale=str(self._time_sys))

        # read inputs
        self._read_data(f_handler)

        f_handler.close()

    def _read_header(self, file, bias_enum: EnumSatelliteBias):
        """
        Method to read header data.

        This function performs the following tasks:
            * checks that the provided file is a valid bias SINEX file of version 1.00
            * checks the bias mode agrees with the user selected `bias_enum`
                * BIAS_MODE set to Absolute for `EnumSatelliteBias.OSB`
                * BIAS_MODE set to Relative for `EnumSatelliteBias.DCB`
            * reads the time system of the time-tags in TIME_SYSTEM
            * ignores all the other header information
        """
        line = file.readline()
        sinex = line[0:5]
        version = float(line[6:10])
        self._bias_mode = line[64]

        if sinex != "%=BIA":
            raise FileError("The provided bias SINEX file {} is not valid. Missing identifier {} in the "
                            "first line".format(file, "%=BIA"))

        if version != 1.0:
            raise FileError("The provided bias SINEX file {} is not valid. "
                            "Only bias SINEX files of version 1.00 are currently supported. A file with version {} "
                            "was provided instead.".format(file, version))

        if self._bias_mode not in ("A", "R"):
            raise FileError("The provided bias SINEX file {} is not valid. "
                            "Unknown bias mode {}. Available modes are absolute (A) or relative (R).".
                            format(file, self._bias_mode))

        if (self._bias_mode == "A" and bias_enum != EnumSatelliteBias.OSB) or \
                (self._bias_mode == "R" and bias_enum != EnumSatelliteBias.DCB):
            raise FileError("The provided bias SINEX file {} is not valid. "
                            "The file bias mode {} does not match the user choice {}.".format(file, self._bias_mode,
                                                                                              bias_enum))

        while line:
            line = file.readline()

            if "TIME_SYSTEM" in line:
                time_sys = line.split()[1]
                self._time_sys = utils.RINEX_SATELLITE_SYSTEM[time_sys]

            elif "+BIAS/SOLUTION" in line:
                # end of header
                break

    def _read_data(self, file):
        """
        Method to read the Satellite Bias Data.
        Station biases are ignored.
        """
        line = " "

        while line:
            line = file.readline()
            tokens = line.split()

            # check for end of block
            if "-BIAS/SOLUTION" in line:
                break
            if tokens[0] == "*BIAS":
                continue  # ignore this line

            bias_type = tokens[0]
            sat_str = tokens[2]

            if bias_type not in ["DSB", "OSB"]:
                # ignore ISB entries
                continue
            if bias_type == "DSB" and self._bias_mode == "A":
                self.log.warning(f"Illegal entry of DSB bias found in file with bias mode set to Absolute. "
                                 f"Ignoring the faulty line: {line}")
                continue
            if bias_type == "OSB" and self._bias_mode == "R":
                self.log.warning(f"Illegal entry of OSB bias found in file with bias mode set to Relative. "
                                 f"Ignoring the faulty line: {line}")
                continue
            if len(sat_str) != 3:
                # station entry -> ignore this line
                continue
            if sat_str[0] not in ("G", "E"):
                # only process GPS and GAL satellites
                continue

            try:
                sat = get_satellite(sat_str)
                if sat.sat_system not in self._active_constellations:
                    # ignore this satellite/constellation due to user configuration choice
                    continue
                obs_list = None
                if bias_type == "DSB":
                    obs_list = (tokens[3], tokens[4])
                elif bias_type == "OSB":
                    obs_list = (tokens[3],)

                epoch_start = self.get_epoch(line[35:49])
                epoch_end = self.get_epoch(line[50:64])

                if epoch_start <= self._epoch < epoch_end:
                    units = line[65:68].strip()
                    value = utils.to_float(line[70:91])

                    # add this bias to the `bias_manager`
                    self.bias_manager.add_bias(epoch_start, epoch_end, bias_type, sat, obs_list, value, units)

            except Exception as e:
                self.log.warning(f"Unable to process line: {line} due to error: {e}")

    def get_epoch(self, epoch_str: str) -> Epoch:
        """ Convert the provided time tag in the file to an Epoch object
        The time tags in these files follow the convention YYYY:DDD:SSSSS, where:
            * YYYY: 4-digit year
            * DDD: 3-digit day in year
            * SSSSS: 5-digit seconds in day
        """
        if len(epoch_str) < 14:
            raise FileError(f"Illegal epoch entry {epoch_str}. Does not follow the defined format yyyy:ddd:ssss.")

        tokens = epoch_str.split(":")
        if len(tokens) != 3:
            raise FileError(f"Illegal epoch entry {epoch_str}. Does not follow the defined format yyyy:ddd:ssss.")

        year = int(tokens[0])
        doy = int(tokens[1])
        seconds = int(tokens[2])
        return Epoch.from_year_and_doy(year, doy, seconds, scale=self._time_sys)
