# Utilities for reading RINEX files

RINEX_FILE_TYPES = {"M": "Meteorological Data",
                    "N": "Navigation Data",
                    "O": "Observation Data",
                    "C": "Clock Data"}

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

RINEX_OBS_TYPES_UNITS = {
    "PR": "m",
    "SPR": "m",
    "CP": "m",  # carrier phase is converted internally from cycles to meters
    "D": "Hz",  # Doppler is provided in Hz
    "S": "dB-Hz"  # default units of SNR
}

def to_float(nmb: str) -> float:
    """
    Convert numbers from RINEX files (in string format) to python floats

    RINEX numbers are provided in scientific notation with the ´E´ or ´D´ char, for example:
        * 4.391254341941D-09
        * 4.391254341941E-09
    """
    try:
        flt = float(nmb)
    except ValueError:
        nmb = nmb.replace("D", "E")
        flt = float(nmb)

    return flt
