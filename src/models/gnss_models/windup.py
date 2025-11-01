""" Module with computation of Phase Windup """

from datetime import timedelta

import numpy as np

from src.constants import PI
from src.io.config import config_dict
from src.models.frames import latlon2dcm_e_enu

_cache = {}


class WindupEntry:
    """ Class to store the cached Phase Windup"""
    def __init__(self, epoch, windup):
        self.epoch = epoch
        self.windup = windup


def _get_delta_prev_cached(epoch, sat):
    if sat not in _cache:
        return False, 0.0
    else:
        rate = config_dict.get("inputs", "rate")
        prev_epoch = epoch - timedelta(seconds=rate)
        cached_data = _cache[sat]
        if cached_data.epoch == prev_epoch:
            return True, cached_data.windup
        # else: there is a gap in the measurement arc, so we reset the windup (0)
    return False, 0.0


def compute_phase_windup(epoch, sat, dcm_b_e, los, lat, long):
    """ Function to compute the phase windup correction for a given satellite.
    Continuity is guaranteed due to the implementation of a cache that stores the previously computed windup value.

    For reference see Chapter 5.5 of:
        * ESA GNSS DATA PROCESSING, Volume I: Fundamentals and Algorithms, J. Sanz Subirana, J.M. Juan Zornoza and
            M. Hernández-Pajares

    Args:
        epoch(src.data_types.date.Epoch): epoch instance
        sat(src.data_types.gnss.satellite.Satellite): satellite instance
        dcm_b_e(numpy.ndarray): rotation matrix from GNSS body frame to ECEF frame
        los(numpy.ndarray): line of sight vector from satellite to receiver
        lat(float): user latitude, in rad
        long(float): user longitude, in rad

    Returns:
        float: computed phase windup, in cycles
    """
    dcm_e_enu = latlon2dcm_e_enu(lat, long)

    # receiver vectors
    a = dcm_e_enu.T @ [1, 0, 0]
    b = dcm_e_enu.T @ [0, 1, 0]

    # satellite vectors
    a_prime = dcm_b_e @ [1, 0, 0]
    b_prime = dcm_b_e @ [0, 1, 0]

    # compute dipoles
    d = a - los * np.dot(los, a) + np.cross(los, b)
    d_prime = a_prime - los * np.dot(los, a_prime) + np.cross(los, b_prime)

    # compute delta windup
    eta = np.dot(los, np.cross(d_prime, d))
    num = np.dot(d_prime, d)
    den = np.linalg.norm(d_prime) * np.linalg.norm(d)
    delta_phase = np.sign(eta) * np.arccos(num / den)

    # get N for continuity
    found, delta_prev = _get_delta_prev_cached(epoch, sat)
    if found:
        N_prev = round((delta_prev - delta_phase) / (2 * PI))
    else:
        N_prev = 0

    windup = delta_phase + 2 * PI * N_prev

    # update cache
    _cache[sat] = WindupEntry(epoch, windup)

    return windup
