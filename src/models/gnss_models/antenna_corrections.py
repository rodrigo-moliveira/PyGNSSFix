""" Module with several utility functions to compute receiver and satellite antenna corrections:
    * Phase Center Offset (PCO) correction.
    * Phase Center Variation (PCV) correction.
    * Antenna Reference Point (ARP) correction.
    * (TO BE DONE): Phase windup correction.
"""


import numpy as np

from src.constants import PI
from src.data_types.gnss import DataType
from src.models.frames import latlon2dcm_e_enu
from src.common_log import get_logger, MODEL_LOG

_warning_cache = set()
_model_logger = get_logger(MODEL_LOG)

# TODO:
#  * NOTE: if CSpice is disabled, the satellite PCO / PCV corrections are not computed.
#  * NOTE: if ITRF data is not provided, the Cspice Sun ephemerides are not rotated.
# TODO: review caches for different satellites and reces.

def receiver_phase_center_correction(receiver_antenna, datatype, los, lat, long, azim, elev):
    """
    Compute the receiver antenna corrections for the given datatype and satellite position.

    The corrections are:
        * Antenna Reference Point (ARP) correction.
        * Phase Center Offset (PCO) correction.
        * Phase Center Variation (PCV) correction.

    Args:
        receiver_antenna (src.data_types.gnss.antenna.ReceiverAntenna): the receiver antenna.
        datatype (DataType): the data type to compute the corrections for.
        los (np.ndarray): the Line-Of-Sight vector from the receiver to the satellite.
        lat (float): the latitude of the receiver.
        long (float): the longitude of the receiver.
        azim (float): the azimuth of the satellite (seen from the receiver).
        elev (float): the elevation of the satellite (seen from the receiver).

    Returns:
        float: the correction in meters
    """
    correction = 0.0
    enu2ecef = latlon2dcm_e_enu(lat, long).T

    # Antenna Reference Point Correction
    arp_offset = receiver_antenna.arp_offset  # in DELTA H/E/N
    if receiver_antenna.arp_enabled:
        arp_offset_enu = np.array([arp_offset[1], arp_offset[2], arp_offset[0]])  # in ENU
        correction += los @ (enu2ecef @ arp_offset_enu)

    # PCO/PCV correction
    if not DataType.is_iono_free(datatype):
        # datatype is not iono free -> use the PCO/PCV for this frequency

        if datatype.freq in receiver_antenna.freq_data:

            # PCO correction
            if receiver_antenna.pco_enabled:
                pco = compute_rec_pco_for_freq(los, enu2ecef, receiver_antenna, datatype.freq)
                correction += pco

            # PCV correction
            if receiver_antenna.pcv_enabled:
                pcv = receiver_antenna.get_pcv(datatype.freq, azim, elev) / 1000.0  # PCV in meters for this freq
                correction += pcv
        else:
            if datatype.freq not in _warning_cache:
                _model_logger.warning(f"Receiver Phase Data (PCO/PCV) correction not found for freq {datatype.freq}.")
                _warning_cache.add(datatype.freq)

    else:
        # datatype is iono free -> compute the PCO/PCV for the iono free combination
        freq1, freq2 = DataType.get_iono_free_base_frequencies(datatype)
        # computations for IF
        f1 = freq1.freq_value
        f2 = freq2.freq_value
        gama1 = f1 * f1 / (f1 * f1 - f2 * f2)
        gama2 = f2 * f2 / (f1 * f1 - f2 * f2)

        if freq1 in receiver_antenna.freq_data and freq2 in receiver_antenna.freq_data:
            # PCO correction
            pco1 = compute_rec_pco_for_freq(los, enu2ecef, receiver_antenna, freq1)
            pco2 = compute_rec_pco_for_freq(los, enu2ecef, receiver_antenna, freq2)
            pco = gama1 * pco1 - gama2 * pco2
            correction += pco

            # PCV correction
            pcv1 = receiver_antenna.get_pcv(freq1, azim, elev) / 1000.0
            pcv2 = receiver_antenna.get_pcv(freq2, azim, elev) / 1000.0
            pcv = gama1 * pcv1 - gama2 * pcv2
            correction += pcv
        else:
            if freq1 not in _warning_cache:
                _model_logger.warning(f"Receiver Phase Data (PCO/PCV) correction not found for freq {freq1}.")
                _warning_cache.add(freq1)
            if freq2 not in _warning_cache:
                _model_logger.warning(f"Receiver Phase Data (PCO/PCV) correction not found for freq {freq2}.")
                _warning_cache.add(freq2)

    return correction


def satellite_phase_center_correction(sat_antenna, dcm_b_e, datatype, los, nadir, azimuth):
    correction = 0.0

    # PCO/PCV correction
    if not DataType.is_iono_free(datatype):
        # datatype is not iono free -> use the PCO/PCV for this frequency
        if datatype.freq in sat_antenna.freq_data:

            # PCO correction
            if sat_antenna.pco_enabled:
                pco = compute_sat_pco_for_freq(los, dcm_b_e, sat_antenna, datatype.freq)
                correction += pco

            # PCV correction
            if sat_antenna.pcv_enabled:
                pcv = sat_antenna.get_pcv(datatype.freq, azimuth, PI/2.0 - nadir) / 1000.0  # PCV in meters for this freq
                correction += pcv
        else:
            if datatype.freq not in _warning_cache:
                _model_logger.warning(f"Phase Data (PCO/PCV) correction not found for freq {datatype.freq} "
                                      f"for satellite {str(sat_antenna.satellite)}")
                _warning_cache.add(datatype.freq)

    else:
        # datatype is iono free -> compute the PCO/PCV for the iono free combination
        freq1, freq2 = DataType.get_iono_free_base_frequencies(datatype)
        # computations for IF
        f1 = freq1.freq_value
        f2 = freq2.freq_value
        gama1 = f1 * f1 / (f1 * f1 - f2 * f2)
        gama2 = f2 * f2 / (f1 * f1 - f2 * f2)

        if freq1 in sat_antenna.freq_data and freq2 in sat_antenna.freq_data:
            # PCO correction
            pco1 = compute_sat_pco_for_freq(los, dcm_b_e, sat_antenna, freq1)
            pco2 = compute_sat_pco_for_freq(los, dcm_b_e, sat_antenna, freq2)
            pco = gama1 * pco1 - gama2 * pco2
            correction += pco

            # PCV correction
            pcv1 = sat_antenna.get_pcv(freq1, azimuth, PI / 2.0 - nadir) / 1000.0
            pcv2 = sat_antenna.get_pcv(freq2, azimuth, PI / 2.0 - nadir) / 1000.0
            pcv = gama1 * pcv1 - gama2 * pcv2
            correction += pcv
        else:
            if freq1 not in _warning_cache:
                _model_logger.warning(f"Phase Data (PCO/PCV) correction not found for freq {freq1}.")
                _warning_cache.add(freq1)
            if freq2 not in _warning_cache:
                _model_logger.warning(f"Phase Data (PCO/PCV) correction not found for freq {freq2}.")
                _warning_cache.add(freq2)

    return correction

def compute_rec_pco_for_freq(los, enu2ecef, receiver_antenna, freq):
    """
    Compute the Phase Center Offset (PCO) correction for the given frequency.

    Args:
        los (np.ndarray): the Line-Of-Sight vector from the receiver to the satellite.
        enu2ecef (np.ndarray): the transformation matrix from ENU to ECEF.
        receiver_antenna (src.data_types.gnss.antenna.ReceiverAntenna): the receiver antenna.
        freq (DataType): the frequency to compute the PCO for.

    Returns:
        float: the PCO correction in meters.
    """
    rec_antenna_pco = receiver_antenna.freq_data[freq].pco  # NORTH / EAST / UP in millimeters
    pco_vec = np.array(
        [rec_antenna_pco[1], rec_antenna_pco[0], rec_antenna_pco[2]]) / 1000.0  # PCO in meters and ENU
    return los @ (enu2ecef @ pco_vec)

def compute_sat_pco_for_freq(los, dcm_b_e, sat_antenna, freq):
    """
    Compute the Phase Center Offset (PCO) correction for the given frequency.

    Args:
        los (np.ndarray): the Line-Of-Sight vector from the receiver to the satellite.
        enu2ecef (np.ndarray): the transformation matrix from ENU to ECEF.
        receiver_antenna (src.data_types.gnss.antenna.ReceiverAntenna): the receiver antenna.
        freq (DataType): the frequency to compute the PCO for.

    Returns:
        float: the PCO correction in meters.
    """
    antenna_pco = sat_antenna.freq_data[freq].pco  # X / Y / Z in millimeters
    pco_vec = np.array(antenna_pco) / 1000.0  # PCO in meters and body frame
    return -los @ (dcm_b_e @ pco_vec)