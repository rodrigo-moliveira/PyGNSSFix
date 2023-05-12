import numpy as np
from math import cos, sin, sqrt, atan2
from src import constants
from src.utils.math_utils import rot2, rot3


def geodetic2cartesian(lat, long, height):
    """
    Convert from geodetic coordinates to cartesian coordinates. Both refer to ECEF frame
    geodetic coordinates are:
        * latitude [rad]
        * longitude [rad]
        * height [m]
    cartesian coordinates are:
        * x, y, z [m]

    Args:
        lat: float [rad]
        long: float [rad]
        height: float [m]
    Return:
        list : [x, y, z] [m]
    """

    a = constants.EARTH_SEMI_MAJOR_AXIS
    e2 = constants.EARTH_ECCENTRICITY_SQ

    # compute prime vertical radius of curvature at a given latitude for a given ellipsoid
    N = a / sqrt(1 - e2 * sin(lat) * sin(lat))

    x = (N + height) * cos(lat) * cos(long)
    y = (N + height) * cos(lat) * sin(long)
    z = ((1 - e2) * N + height) * sin(lat)

    return [x, y, z]


def latlon2dcm_e_n(lat, lon):
    """
    transformation matrix from the ECEF frame to the NED frame defined by lat and lon.
    Args:
        lat: latitude, rad
        lon: longitude, rad
    """
    # m = np.array([[-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)],
    #              [-sin(lon), cos(lon), 0],
    #              [-cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)]])

    return rot2(-np.pi / 2.0 - lat) @ rot3(lon)


def cartesian2geodetic(x, y, z):
    """
    Convert from cartesian coordinates to geodetic coordinates. Both refer to ECEF frame
    geodetic coordinates are:
        * latitude [rad]
        * longitude [rad]
        * height [m]
    cartesian coordinates are:
        * x, y, z [m]
    Algorithm B.1.2:  From Cartesian to Ellipsoidal Coordinates from **REF[1]**


    Args:
        x: float [m]
        y: float [m]
        z: float [m]
    Return:
        list : [lat, long, height] [rad,rad,m]
    """
    e2 = constants.EARTH_ECCENTRICITY_SQ
    a = constants.EARTH_SEMI_MAJOR_AXIS

    # computation of longitude
    long = atan2(y, x)

    # initial value of latitude
    p = sqrt(x * x + y * y)
    lat = 0 if p == 0 else atan2(z / p, 1 - e2)

    # iterative process to refine latitude
    MAX_ITERS = 10
    ABS_TOL = 1E-10
    i = 0
    height = 0
    while i < MAX_ITERS:
        lat_prev = lat

        N = a / sqrt(1 - e2 * sin(lat) * sin(lat))
        height = p / cos(lat) - N
        lat = 0 if p == 0 else atan2(z / p, 1 - N / (N + height) * e2)

        i += 1
        if abs(lat - lat_prev) < ABS_TOL:
            break

    return [lat, long, height]


def ecef2ned(arr_rover_ecef, arr_obs_ecef, origin_llh):
    """
    Convert from ECEF frame to NED frame (cartesian coordinates). The NED frame is centered (origin O) at the ground
    observer (for example, the GNSS receiver)
    NED stands for North East Down -> local components from the origin O to the rover

    Args:
        arr_rover_ecef (numpy.ndarray) : position/velocity of rover in ecef
        arr_obs_ecef (numpy.ndarray) : reference or observer position/velocity
        origin_llh (numpy.ndarray or list) : lat long height coordinates of the local topocentric frame origin
    Return:
        numpy.ndarray : [x_enu, y_enu, z_enu]  [m] or [m/s]
    """
    lat = origin_llh[0]
    lon = origin_llh[1]

    # get rotation matrix from ECEF to NED
    R = latlon2dcm_e_n(lat, lon)

    # apply transform
    arr_enu = R @ (arr_rover_ecef - arr_obs_ecef)

    return arr_enu
