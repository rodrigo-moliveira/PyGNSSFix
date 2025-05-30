""" This module provides functions to compute the required ephemerides of celestial bodies (namely, the Sun). """

import spiceypy
from src.data_types.date import Epoch
from src.io.config import config_dict
from src.models.frames.itrf_transformations import ITRF_Transformation, ITRFError


def compute_sun_pos(date, obs="EARTH"):
    """
    Compute the position of the Sun at the given date in the ECEF Frame.

    Applies the ITRF Transformation from the CSpice Frame (typically ITRF93) to the IGS Frame
    (typically ITRF2014 or ITRF2020 for the newer scenarios).

    Args:
        date(Epoch): the date for the computation
        obs(str): the observer of the ephemeris. Default is EARTH.

    Returns:
        numpy.ndarray: the position of the Sun in the given frame.

    """
    if not isinstance(date, Epoch):
        raise TypeError("The date must be an instance of the Epoch class.")
    decimal_year = date.decimal_year
    et = date.ephemeris_time

    # Create or load static ITRF Transformation class.
    cspice_frame = config_dict.get("inputs", "cspice_kernels", "cspice_ecef_frame")
    target_frame = config_dict.get("inputs", "IGS_ecef_frame")
    result = spiceypy.spkpos('SUN', et, cspice_frame, 'NONE', obs)
    sun_pos = result[0]

    try:
        tf = ITRF_Transformation.factory(cspice_frame, target_frame)
        sun_pos = tf.transform(sun_pos, decimal_year)
    except ITRFError as e:
        from src.common_log import MODEL_LOG, get_logger
        log = get_logger(MODEL_LOG)
        log.error(f"Unable to apply ITRF Transformation in the computation of CSpice Sun Ephemerides: {e}")

    return sun_pos
