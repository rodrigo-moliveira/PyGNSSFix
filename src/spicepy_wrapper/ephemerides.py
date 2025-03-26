""" This module provides functions to compute the required ephemerides of celestial bodies (namely, the Sun). """

import spiceypy
from src.data_types.date import Epoch


def compute_sun_pos(date, frame="ITRF93", obs="EARTH"):
    """
    Compute the position of the Sun at the given date in the ECEF Frame.

    Args:
        date(Epoch or float): the date for the computation, either as an Epoch instance or
            directly in float ephemeris time
        frame(str): the reference frame of the output position. Default is ITRF93.
        obs(str): the observer of the ephemeris. Default is EARTH.

    Returns:
        numpy.ndarray: the position of the Sun in the given frame.

    """
    if isinstance(date, Epoch):
        et = date.ephemeris_time
    else:
        et = date

    sun_pos = spiceypy.spkpos('SUN', et, frame, 'NONE', obs)
    # TODO convert to ITRF2020
    return sun_pos[0]
