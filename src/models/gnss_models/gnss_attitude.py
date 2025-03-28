""" Module for GNSS Attitude Models """

import numpy


def gnss_attitude(gnss_pos: numpy.ndarray, sun_pos: numpy.ndarray):
    """ Compute the usual yaw-steering attitude model for GNSS satellites.

    The IGS definition for the body-fixed (satellite) frame is:
        * coordinate system’s origin is at the satellite’s CoM
        * z-axis points toward the Earth center (bore sight direction of the antenna)
        * y-axis (rotation axis of the solar panels) corresponds to the cross product of the
            z-axis with the vector from the satellite to the Sun
        * x-axis completes the right-handed system (x cross y = z)

    This function computes the rotation matrix from body-fixed frame to ECEF frame.

    Args:
        gnss_pos (numpy.ndarray): the ECEF position of the GNSS satellite.
        sun_pos (numpy.ndarray): the ECEF position of the Sun.

    Returns:
        numpy.ndarray: the rotation matrix (DCM) from the body-fixed frame to the E
    """
    e_z = -gnss_pos / numpy.linalg.norm(gnss_pos)
    e_s = (sun_pos - gnss_pos) / numpy.linalg.norm(sun_pos - gnss_pos)
    e_y = numpy.cross(e_z, e_s) / numpy.linalg.norm(numpy.cross(e_z, e_s))
    e_x = numpy.cross(e_y, e_z)
    dcm = numpy.array([e_x, e_y, e_z]).T
    return dcm
