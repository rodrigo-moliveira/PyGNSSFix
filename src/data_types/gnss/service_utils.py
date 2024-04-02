"""Utility maps and dicts to deal with the constellations and observations.
"""

from .data_type import get_data_type


AvailableConstellations = {"GPS", "GAL", "GLO", "BDS"}

# Rinex format
GPSServices = {'1C', '1S', '1L', '1X', '1P', '1W', '1Y', '1M',
               '2C', '2D', '2S', '2L', '2X', '2P', '2W', '2Y', '2M',
               '5I', '5Q', '5X'}
GALServices = {'1A', '1B', '1C', '1X', '1Z',
               '5I', '5Q', '5X',
               '7I', '7Q', '7X',
               '8I', '8Q', '8X',
               '6A', '6B', '6C', '6X', '6Z'}

Services = {
    "GPS": GPSServices,
    "GAL": GALServices
}

ConstellationToCodeMap = {"GPS": "G", "GAL": "E", "GLO": "R", "BDS": "C", "UNKNOWN": "U"}

CodeToConstellationMap = {"G": "GPS",
                          "E": "GAL",
                          "R": "GLO",
                          "C": "BDS"}


def get_code_from_service(service, constellation):
    datatype = None

    if constellation == "GPS":
        if service in GPSServices:
            datatype = get_data_type("C" + service[0], constellation)

    elif constellation == "GAL":
        if service in GALServices:
            datatype = get_data_type("C" + service[0], constellation)

    return datatype


def get_carrier_from_service(service, constellation):
    datatype = None

    if constellation == "GPS":
        if service in GPSServices:
            datatype = get_data_type("L" + service[0], constellation)

    elif constellation == "GAL":
        if service in GALServices:
            datatype = get_data_type("L" + service[0], constellation)

    return datatype
