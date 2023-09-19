from . import utils
from src.data_types.date.date import Epoch
from src.data_types.gnss.satellite import get_satellite
from src.errors import FileError

from src.data_types.gnss.navigation_data import NavigationData, NavigationPoint
from ... import WORKSPACE_PATH

"""
Example of File Rinex Navigation File V3.03 (header + data):

     3.03           N: GNSS NAV DATA    G: GPS              RINEX VERSION / TYPE
Converto v3.5.1     IGN                 20190115 003605 UTC PGM / RUN BY / DATE
GPSA   7.4506E-09 -1.4901E-08 -5.9605E-08  1.1921E-07       IONOSPHERIC CORR
GPSB   8.8064E+04 -4.9152E+04 -1.9661E+05  3.2768E+05       IONOSPHERIC CORR
GPUT -1.8626451492E-09-2.664535259E-15 233472 2036          TIME SYSTEM CORR
    18    18  1929     7                                    LEAP SECONDS
                                                            END OF HEADER
G01 2019 01 13 16 00 00-1.447494141757D-04-6.366462912410D-12 0.000000000000D+00
     3.100000000000D+01-3.281250000000D+00 4.391254341941D-09 5.856830273633D-01
    -2.272427082062D-07 8.282583905384D-03 5.565583705902D-06 5.153660743713D+03
     5.760000000000D+04-2.607703208923D-08-5.974139377628D-01 1.303851604462D-07
     9.741765399631D-01 2.808750000000D+02 6.890211188724D-01-8.210699151434D-09
     1.610781381224D-10 1.000000000000D+00 2.036000000000D+03 0.000000000000D+00
     2.000000000000D+00 0.000000000000D+00 5.587935447693D-09 3.100000000000D+01
     5.427000000000D+04    

Important Note: The time tags of the navigation messages (e.g., time of ephemeris, time of clock) are given in the 
respective satellite time systems!

"""


class RinexNavReader:
    """
    Class RinexNavReader

    Attributes
        ----------
        nav : NavigationData
    """

    def __init__(self, file, nav):

        # instance variables
        self.nav = nav
        self._first_epoch = None
        self._first_epoch_set = False

        f_handler = open(f"{WORKSPACE_PATH}/{file}", "r")

        # read header
        self._read_header(f_handler)

        # read inputs
        self._read_data(f_handler)

        self.nav.header.first_epoch = self._first_epoch

        f_handler.close()

    def _read_header(self, file):
        """
        Method to read header data

        Tags to look for:
            * RINEX VERSION / TYPE  -> rinex_version and satellite_system
            * IONOSPHERIC CORR      -> iono_corrections
            * LEAP SECONDS          -> leap_seconds
            * END OF HEADER         -> end of header section
        """
        line = " "

        while line:
            line = file.readline()

            if "RINEX VERSION / TYPE" in line:
                self.nav.header.rinex_version = utils.to_float(line[5:10])
                if self.nav.header.rinex_version < 3:
                    raise FileError("The provided rinex file {} is of version {}. Only version 3.00 or "
                                    "higher is supported. Error!".format(file, self.nav.header.rinex_version))

                rinex_type = line[20]
                if rinex_type != 'N':
                    raise FileError("Rinex File {} should be a GNSS Navigation Data File. Instead, "
                                    "a {} was provided (code {})".format(file,
                                                                         utils.RINEX_FILE_TYPES.get
                                                                         (rinex_type, "Unknown Data File"), rinex_type))

                self.nav.header.satellite_system = utils.RINEX_SATELLITE_SYSTEM.get(line[40], "UNKNOWN")

            elif "LEAP SECONDS" in line:
                data = line[:utils.RINEX_OBS_END_OF_DATA_HEADER].split()
                self.nav.header.leap_seconds = int(data[0])

            elif "IONOSPHERIC CORR" in line:
                data = line[:utils.RINEX_OBS_END_OF_DATA_HEADER].split()
                if data[0] in {"GAL", "GPSA", "GPSB"}:
                    self.nav.header.iono_corrections[data[0]] = [utils.to_float(i) for i in data[1:]]

            elif "END OF HEADER" in line:
                break

    def _read_data(self, file):
        """
        Read GPS navigation data

        Data to fetch:
            * satellite, toc, af0, af1, af2                         [Line 1]
            * IODE, crs, deltaN, M0                                 [Line 2]
            * Cuc, e, Cus, sqrtA                                    [Line 3]
            * Toe, Cic, RAAN0, Cis,                                 [Line 4]
            * i0, crc, omega, RAANDot,                              [Line 5]
            * iDot, codesL2, toe (sensors week), flagL2,            [Line 6]
            * SV_URA, SV_health, TGD, IODC                          [Line 7]
            * TransmissionTime (seconds of week)                    [Line 8]

        Control variables:
            * SV_URA - User Range Accuracy, in meters (the message data is already in meters!!! not the ura index)
                    URA INDEX       URA (meters)
                    0           0.00 < URA ≤ 2.40
                    1           2.40 < URA ≤ 3.40
                    2           3.40 < URA ≤ 4.85
                    3           4.85 < URA ≤ 6.85
                    4           6.85 < URA ≤ 9.65
                    5           9.65 < URA ≤ 13.65
                    6           13.65 < URA ≤ 24.00
                    7           24.00 < URA ≤ 48.00
                    8           48.00 < URA ≤ 96.00
                    9           96.00 < URA ≤ 192.00
                    10          192.00 < URA ≤ 384.00
                    11          384.00 < URA ≤ 768.00
                    12          768.00 < URA ≤ 1536.00
                    13          1536.00 < URA ≤ 3072.00
                    14          3072.00 < URA ≤ 6144.00
                    15          6144.00 < URA (or no accuracy prediction is available - unauthorized users are
                                advised to use the SV at their own risk.)

            * SV_health : Satellite health status:
                    0 = all NAV data are OK,
                    != 0 some or all NAV data are bad (ignore this entry!!).
        """
        line = " "

        while line:
            line = file.readline()

            if len(line) == 0:
                break

            if line[0] == "G":
                # new GPS navigation epoch inputs
                navMessage = NavigationPoint()

                # read 1st line
                satellite = get_satellite(line[0:3])
                year = int(line[4:8])
                month = int(line[9:11])
                day = int(line[12:14])
                hour = int(line[15:17])
                minute = int(line[18:20])
                second = int(line[21:23])

                toc = Epoch(year, month, day, hour, minute, second,
                            scale=str(satellite.sat_system))  # time_system="GPS"

                # set first epoch for this file
                if not self._first_epoch_set:
                    self._first_epoch = toc
                    self._first_epoch_set = True

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
