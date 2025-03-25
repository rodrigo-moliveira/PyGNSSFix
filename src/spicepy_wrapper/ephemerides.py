""" This module provides functions to compute the required ephemerides and attitude of celestial bodies (Sun). """

import spiceypy
from src.data_types.date import Epoch

def compute_body_frame():
    # TODO: apply equation to compute body frame matrix
    pass

def compute_sun_pos(date):
    """
    Compute the position of the Sun at the given date.
    :param date: datetime object
    :return: position of the Sun in ECEF coordinates
    """
    # Compute the position of the Sun at the given date
    #sun_pos = spiceypy.spkpos('SUN', date, 'J2000', 'NONE', 'EARTH')
    sun_pos = spiceypy.spkpos('SUN', date, 'ITRF93', 'NONE', 'EARTH')
    # TODO convert to ITRF2020
    return sun_pos[0]

if __name__ == "__main__":
    from src.spicepy_wrapper import setup_cspice
    setup_cspice(".")
    # Example date
    utc = Epoch(2019, 10, 1, 12, 13, 39, 0, scale="UTC")
    print(utc.ephemeris_time)

    # Compute the position of the Sun
    sun_pos = compute_sun_pos(utc.ephemeris_time)

    # Print the result
    print(f"Position of the Sun on {utc}: {sun_pos}")
    print(sun_pos[0])
