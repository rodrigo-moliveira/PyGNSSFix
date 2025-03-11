import numpy as np

from src.models.frames import latlon2dcm_e_enu
def compute_receiver_correction(receiver_antenna, datatype, los, lat, long, azim, elev):
    enu2ecef = latlon2dcm_e_enu(lat, long).T

    # Antenna Reference Point Correction
    arp_offset = receiver_antenna.arp_offset  # in DELTA H/E/N
    if receiver_antenna.arp_enabled:
        arp_offset_enu = np.array([arp_offset[1], arp_offset[2], arp_offset[0]])  # in ENU
    else:
        arp_offset_enu = np.zeros(3)

    # TODO: check if the antenna has a PCO for this frequency (it may even be iono free type...)
    # PCO correction
    if receiver_antenna.pco_enabled:
        rec_antenna_pco = receiver_antenna.freq_data[datatype.freq].pco  # NORTH / EAST / UP in millimeters
        pco_vec = np.array(
            [rec_antenna_pco[1], rec_antenna_pco[0], rec_antenna_pco[2]]) / 1000  # PCO vector in meters (ENU)
    else:
        pco_vec = np.zeros(3)

    # PCV correction
    if receiver_antenna.pcv_enabled:
        pcv = receiver_antenna.get_pcv(datatype.freq, azim, elev)/1000  # PCV in meters for this freq
    else:
        pcv = 0

    correction = los @ (enu2ecef @ arp_offset_enu) + los @ (enu2ecef @ pco_vec) + pcv
    return correction
