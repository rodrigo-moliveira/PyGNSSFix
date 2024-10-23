""" Parser of RINEX Navigation files (version 3.XX)
The official documentation of these files may be found in https://igs.org/wg/rinex/#documents-formats
"""

from . import utils
from src.data_types.date import Epoch, EopDb
from src.data_types.gnss import get_satellite
from src.errors import FileError
from ...data_mng.gnss.navigation_data import NavigationPointGPS, NavigationPointGAL, NavigationData
from src import WORKSPACE_PATH
from src.common_log import IO_LOG, get_logger


class RinexNavReader:
    """
    Parser of Rinex Navigation files
    """

    def __init__(self, file, nav, gal_nav_type):
        """
        Reads the provided navigation file and stores its content in the `NavigationData` instance.

        Args:
            file(str): path to the input RINEX Navigation file to load
            nav(NavigationData): the `NavigationData` object to store the navigation information extracted from the file
            gal_nav_type(str): the GAL navigation message (I/NAV of F/NAV), as configured by the user
        """

        # instance variables
        self.nav = nav
        self.gal_nav_type = gal_nav_type

        f_handler = open(f"{WORKSPACE_PATH}/{file}", "r")
        self.log = get_logger(IO_LOG)
        self.log.info(f"Reading navigation file {WORKSPACE_PATH}/{file}...")

        # read header
        self._read_header(f_handler)

        # set time corrections in EopDb
        EopDb.set_time_correction(self.nav.header.time_correction)

        # read inputs
        self._read_data(f_handler)

        f_handler.close()

    def _read_header(self, file):
        """
        Method to read header data

        Header lines to look for (others are ignored):
            * RINEX VERSION / TYPE  -> rinex_version and satellite_system
            * IONOSPHERIC CORR      -> iono_corrections
            * LEAP SECONDS          -> leap_seconds
            * TIME SYSTEM CORR      -> read GPGA field (GPS-to-Galileo time offset)
            * END OF HEADER         -> end of header section
        """
        line = " "

        while line:
            line = file.readline()

            if "RINEX VERSION / TYPE" in line:
                rinex_version = utils.to_float(line[5:10])
                if rinex_version < 3 or rinex_version >= 4:
                    raise FileError("The provided rinex_parser file {} is of version {}. Only versions 3.XX "
                                    "are supported. Error!".format(file, self.nav.header.rinex_version))

                rinex_type = line[20]
                if rinex_type != 'N':
                    raise FileError("Rinex File {} should be a GNSS Navigation Data File. Instead, "
                                    "a {} was provided (code {})".format(file,
                                                                         utils.RINEX_FILE_TYPES.get
                                                                         (rinex_type, "Unknown Data File"), rinex_type))

            elif "LEAP SECONDS" in line:
                data = line[:utils.RINEX_OBS_END_OF_DATA_HEADER].split()
                self.nav.header.leap_seconds = int(data[0])

            elif "IONOSPHERIC CORR" in line:
                data = line[:utils.RINEX_OBS_END_OF_DATA_HEADER].split()
                if data[0] in {"GAL", "GPSA", "GPSB"}:
                    self.nav.header.iono_corrections[data[0]] = [utils.to_float(i) for i in data[1:]]

            elif "TIME SYSTEM CORR" in line:
                # only save GPGA (or GGTO)
                if "GPGA" in line:
                    a0 = utils.to_float(line[5:22])
                    a1 = utils.to_float(line[22:39])
                    T_ref = int(line[39:46])
                    week_ref = int(line[46:53])

                    self.nav.header.time_correction["GGTO"] = [a0, a1, T_ref, week_ref]

            elif "END OF HEADER" in line:
                break

    def _read_data(self, file):
        """
        Read GPS and GAL navigation data.

        Important Note: The time tags of the navigation messages (e.g., time of ephemeris, time of clock) are given in
            the respective satellite time systems!
        """
        line = " "

        while line:
            line = file.readline()

            if len(line) == 0:
                break

            if line[0] == "G":
                # new GPS navigation epoch inputs
                navMessage = NavigationPointGPS()

                # read 1st line
                satellite = get_satellite(line[0:3])
                year = int(line[4:8])
                month = int(line[9:11])
                day = int(line[12:14])
                hour = int(line[15:17])
                minute = int(line[18:20])
                second = int(line[21:23])

                toc = Epoch(year, month, day, hour, minute, second,
                            scale=str(satellite.sat_system))

                af0 = utils.to_float(line[23:42])
                af1 = utils.to_float(line[42:61])
                af2 = utils.to_float(line[61:80])
                setattr(navMessage, "satellite", satellite)
                setattr(navMessage, "toc", toc)
                setattr(navMessage, "af0", af0)
                setattr(navMessage, "af1", af1)
                setattr(navMessage, "af2", af2)

                # read 2nd line
                line = file.readline()
                iode = utils.to_float(line[4:23])
                crs = utils.to_float(line[23:42])
                deltaN = utils.to_float(line[42:61])
                M0 = utils.to_float(line[61:80])
                setattr(navMessage, "IODE", iode)
                setattr(navMessage, "crs", crs)
                setattr(navMessage, "deltaN", deltaN)
                setattr(navMessage, "M0", M0)

                # read 3rd line
                line = file.readline()
                cuc = utils.to_float(line[4:23])
                eccentricity = utils.to_float(line[23:42])
                cus = utils.to_float(line[42:61])
                sqrtA = utils.to_float(line[61:80])
                setattr(navMessage, "cuc", cuc)
                setattr(navMessage, "eccentricity", eccentricity)
                setattr(navMessage, "cus", cus)
                setattr(navMessage, "sqrtA", sqrtA)

                # read 4th line
                line = file.readline()
                sec_toe = utils.to_float(line[1:23])
                cic = utils.to_float(line[23:42])
                RAAN0 = utils.to_float(line[42:61])
                cis = utils.to_float(line[61:80])
                setattr(navMessage, "cic", cic)
                setattr(navMessage, "RAAN0", RAAN0)
                setattr(navMessage, "cis", cis)

                # read 5th line
                line = file.readline()
                inclination0 = utils.to_float(line[4:23])
                crc = utils.to_float(line[23:42])
                omega = utils.to_float(line[42:61])
                RAANDot = utils.to_float(line[61:80])
                setattr(navMessage, "i0", inclination0)
                setattr(navMessage, "crc", crc)
                setattr(navMessage, "omega", omega)
                setattr(navMessage, "RAANDot", RAANDot)

                # read 6th line
                line = file.readline()
                iDot = utils.to_float(line[4:23])
                codesL2 = utils.to_float(line[23:42])
                week_toe = utils.to_float(line[42:61])
                flagL2 = utils.to_float(line[61:80])
                setattr(navMessage, "toe", [week_toe, sec_toe])
                setattr(navMessage, "iDot", iDot)
                setattr(navMessage, "codesL2", codesL2)
                setattr(navMessage, "flagL2", flagL2)

                # read 7th line
                line = file.readline()
                SV_URA = utils.to_float(line[4:23])
                SV_health = utils.to_float(line[23:42])
                TGD = utils.to_float(line[42:61])
                IODC = utils.to_float(line[61:80])
                setattr(navMessage, "SV_URA", SV_URA)
                setattr(navMessage, "SV_health", SV_health)
                setattr(navMessage, "TGD", TGD)
                setattr(navMessage, "IODC", IODC)

                # read 8th line
                line = file.readline()
                TransmissionTime = utils.to_float(line[4:23])
                setattr(navMessage, "TransmissionTime", TransmissionTime)

                self.nav.set_data(toc, satellite, navMessage)

            if line[0] == "E":
                # new GAL navigation epoch inputs
                navMessage = NavigationPointGAL()

                # read 1st line
                satellite = get_satellite(line[0:3])
                year = int(line[4:8])
                month = int(line[9:11])
                day = int(line[12:14])
                hour = int(line[15:17])
                minute = int(line[18:20])
                second = int(line[21:23])

                toc = Epoch(year, month, day, hour, minute, second,
                            scale=str(satellite.sat_system))

                af0 = utils.to_float(line[23:42])
                af1 = utils.to_float(line[42:61])
                af2 = utils.to_float(line[61:80])
                setattr(navMessage, "satellite", satellite)
                setattr(navMessage, "toc", toc)
                setattr(navMessage, "af0", af0)
                setattr(navMessage, "af1", af1)
                setattr(navMessage, "af2", af2)

                # read 2nd line
                line = file.readline()
                IODnav = utils.to_float(line[4:23])
                crs = utils.to_float(line[23:42])
                deltaN = utils.to_float(line[42:61])
                M0 = utils.to_float(line[61:80])
                setattr(navMessage, "IODnav", IODnav)
                setattr(navMessage, "crs", crs)
                setattr(navMessage, "deltaN", deltaN)
                setattr(navMessage, "M0", M0)

                # read 3rd line
                line = file.readline()
                cuc = utils.to_float(line[4:23])
                eccentricity = utils.to_float(line[23:42])
                cus = utils.to_float(line[42:61])
                sqrtA = utils.to_float(line[61:80])
                setattr(navMessage, "cuc", cuc)
                setattr(navMessage, "eccentricity", eccentricity)
                setattr(navMessage, "cus", cus)
                setattr(navMessage, "sqrtA", sqrtA)

                # read 4th line
                line = file.readline()
                sec_toe = utils.to_float(line[1:23])
                cic = utils.to_float(line[23:42])
                RAAN0 = utils.to_float(line[42:61])
                cis = utils.to_float(line[61:80])
                setattr(navMessage, "cic", cic)
                setattr(navMessage, "RAAN0", RAAN0)
                setattr(navMessage, "cis", cis)

                # read 5th line
                line = file.readline()
                inclination0 = utils.to_float(line[4:23])
                crc = utils.to_float(line[23:42])
                omega = utils.to_float(line[42:61])
                RAANDot = utils.to_float(line[61:80])
                setattr(navMessage, "i0", inclination0)
                setattr(navMessage, "crc", crc)
                setattr(navMessage, "omega", omega)
                setattr(navMessage, "RAANDot", RAANDot)

                # read 6th line
                line = file.readline()
                iDot = utils.to_float(line[4:23])
                dataSrc = utils.to_float(line[23:42])
                week_toe = utils.to_float(line[42:61])
                setattr(navMessage, "toe", [week_toe, sec_toe])
                setattr(navMessage, "iDot", iDot)
                setattr(navMessage, "dataSrc", dataSrc)

                # read 7th line
                line = file.readline()
                SISA = utils.to_float(line[4:23])
                SV_health = utils.to_float(line[23:42])
                BGDe5a = utils.to_float(line[42:61])
                BGDe5b = utils.to_float(line[61:80])
                setattr(navMessage, "SISA", SISA)
                setattr(navMessage, "SV_health", SV_health)
                setattr(navMessage, "BGDE1E5a", BGDe5a)
                setattr(navMessage, "BGDE1E5b", BGDe5b)

                # read 8th line
                line = file.readline()
                TransmissionTime = utils.to_float(line[4:23])
                setattr(navMessage, "TransmissionTime", TransmissionTime)

                try:
                    nav_type = navMessage.find_message_type()
                    if nav_type == self.gal_nav_type:
                        self.nav.set_data(toc, satellite, navMessage)
                    # else: nav message is skipped
                except Exception as e:
                    self.log.warning(f"Error reading satellite {satellite} ephemeride for epoch {toc}: {e}")
