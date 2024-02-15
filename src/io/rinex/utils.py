# Utilities for reading RINEX files

RINEX_FILE_TYPES = {"M": "Meteorological Data",
                    "N": "Navigation Data",
                    "O": "Observation Data"}

RINEX_SATELLITE_SYSTEM = {"G": "GPS",
                          "R": "GLO",
                          "E": "GAL",
                          "J": "QZSS",
                          "C": "BDS",
                          "I": "IRNSS",
                          "S": "SBAS payload",
                          "M": "Mixed"}

# RINEX OBS
RINEX_OBS_END_OF_DATA_HEADER = 60
RINEX_OBS_TYPES_TO_READ = {"C",  # pseudo range
                           "L",  # carrier phase
                           "S",  # signal-to-noise ratio
                           "D"}  # doppler


def to_float(nmb: str):
    """
    convert numbers from rinex files to python floats

    rinex floats can either be:
        4.391254341941D-09
        4.391254341941E-09
    that is, the scientific notation can use the ´E´ or ´D´ char
    """
    try:
        flt = float(nmb)
    except ValueError:
        nmb = nmb.replace("D", "E")
        flt = float(nmb)

    return flt
