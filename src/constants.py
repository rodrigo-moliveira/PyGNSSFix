""" Definition of some useful constants. """

import numpy as np

"""
Reference:
    [1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
        Springer Cham, 2017 
"""

#######################
# RAD - DEG Conversions
#######################
PI = 3.141592653589793
DEG2RAD = PI / 180
RAD2DEG = 180 / PI

#######################################
# Orbital Mechanics and Earth constants
#######################################
MU_WGS84 = 3.986005E14  # [m^3/sec^2] (WGS84 value, for GPS)
MU_GTRF = 3.986004418e14  # (GTRF value, for GAL)
EARTH_ROTATION = 7.292115E-5  # [rad/sec]
SPEED_OF_LIGHT = 299792458  # [m/s]
EARTH_ANGULAR_RATE = np.array([0.0, 0.0, EARTH_ROTATION])
EARTH_FLATNESS = 1 / 298.257223563
EARTH_ECCENTRICITY_SQ = 2 * EARTH_FLATNESS - EARTH_FLATNESS * EARTH_FLATNESS
EARTH_SEMI_MAJOR_AXIS = 6378137.0  # [m]  Equatorial radius Re
EARTH_J2 = 1.08262668355315130e-3  # J2 harmonic
EARTH_MASS = 5.97237e24  # [kg]
EARTH_G0 = 9.80665  # [m/s^2] standard value cf. https://www.convertunits.com/from/g-unit/to/m/s%5E2

##################
# Time Definitions
##################
SECONDS_IN_DAY = 86400  # In a Julian day
MINUTES_IN_DAY = 1440  # In a Julian day
MINUTES_IN_HOUR = 60.0  # In a Julian day
HOURS_IN_DAY = 24.0
SECONDS_IN_GNSS_WEEK = 604800
DAYS_PER_WEEK = 7
AVERAGE_DAYS_IN_YEAR = 365.25
SECONDS_IN_HOUR = 3600


#######################
# GNSS Frequencies (Hz)
#######################
# GPS Frequencies
GPS_L1_FREQ = 1575.42e6
GPS_L2_FREQ = 1227.60e6
GPS_L5_FREQ = 1176.45e6

# GAL Frequencies
GAL_E1_FREQ = 1575.42e6
GAL_E5A_FREQ = 1176.45e6
GAL_E5B_FREQ = 1207.14e6
GAL_E5ALTBOC_FREQ = 1191.795e6
GAL_E6_FREQ = 1278.75e6

# Iono-Free Frequencies
# See Eqs. (20.45) and (20.46) from [1]
# f0 = 10.23 MHz
# GPS L1 = 154 * f0
# GPS L2 = 120 * f0
# GPS L5 = 115 * f0
# i_A / i_B = -f_A / f_B
# -> 154 / 120 = 77 / 60 => For L1-L2 i_A = 77, i_B = 60
# -> 154 / 115 is irreducible => For L1-L5 i_A = 154, i_B = 115
GPS_L1L2_FREQ = 77 * GPS_L1_FREQ - 60 * GPS_L2_FREQ  # L1-L2 IF frequency
GPS_L1L5_FREQ = 154 * GPS_L1_FREQ - 115 * GPS_L5_FREQ  # L1-L5 IF frequency

# GAL E1 = 154 * f0
# GAL E5A = 115 * f0
# GAL E5B = 118 * f0
# GAL E6 = 125 * f0
# GAL E5AltBOC = 116.5 * f0. NOTE: This is not integer
# -> 154 / 115 is irreducible => For E1-E5A i_A = 154, i_B = 115
# -> 154 / 118 = 77 / 59 => For E1-E5B i_A = 77, i_B = 59
# -> 154 / 125 is irreducible => For E1-E6 i_A = 154, i_B = 125
# -> 154 / 116.5 is irreducible => For E1-E5AltBOC i_A = 154, i_B = 116.5
GAL_E1E5a_FREQ = 154 * GAL_E1_FREQ - 115 * GAL_E5A_FREQ  # E1-E5a IF Frequency
GAL_E1E5b_FREQ = 77 * GAL_E1_FREQ - 59 * GAL_E5B_FREQ  # E1-E5b IF Frequency
GAL_E1E6_FREQ = 154 * GAL_E1_FREQ - 125 * GAL_E6_FREQ  # E1-E6 IF Frequency
GAL_E1E5AltBOC_FREQ = 154 * GAL_E1_FREQ - 116.5 * GAL_E5ALTBOC_FREQ  # E1-E5AltBOC IF Frequency

# Wide-Lane Frequencies
# See Eq. (20.23) of [1]
GPS_WL_L1L2_FREQ = GPS_L1_FREQ - GPS_L2_FREQ
GPS_WL_L1L5_FREQ = GPS_L1_FREQ - GPS_L5_FREQ
GAL_WL_E1E5a_FREQ = GAL_E1_FREQ - GAL_E5A_FREQ
GAL_WL_E1E5b_FREQ = GAL_E1_FREQ - GAL_E5B_FREQ
GAL_WL_E1E6_FREQ = GAL_E1_FREQ - GAL_E6_FREQ
GAL_WL_E1E5AltBOC_FREQ = GAL_E1_FREQ - GAL_E5ALTBOC_FREQ

# Narrow-Lane Frequencies
# See Eq. (20.26) of [1]
GPS_NL_L1L2_FREQ = GPS_L1_FREQ + GPS_L2_FREQ
GPS_NL_L1L5_FREQ = GPS_L1_FREQ + GPS_L5_FREQ
GAL_NL_E1E5a_FREQ = GAL_E1_FREQ + GAL_E5A_FREQ
GAL_NL_E1E5b_FREQ = GAL_E1_FREQ + GAL_E5B_FREQ
GAL_NL_E1E6_FREQ = GAL_E1_FREQ + GAL_E6_FREQ
GAL_NL_E1E5AltBOC_FREQ = GAL_E1_FREQ + GAL_E5ALTBOC_FREQ

# Geometry-free combination
# GPS_GF_L1L2 = SPEED_OF_LIGHT / (SPEED_OF_LIGHT/GPS_L2_FREQ - SPEED_OF_LIGHT/GPS_L1_FREQ)
GPS_GF_L1L2_FREQ = (GPS_L1_FREQ*GPS_L2_FREQ)/(GPS_L1_FREQ-GPS_L2_FREQ)
GPS_GF_L1L5_FREQ = (GPS_L1_FREQ*GPS_L5_FREQ)/(GPS_L1_FREQ-GPS_L5_FREQ)
GAL_GF_E1E5a_FREQ = (GAL_E1_FREQ*GAL_E5A_FREQ)/(GAL_E1_FREQ-GAL_E5A_FREQ)
GAL_GF_E1E5b_FREQ = (GAL_E1_FREQ*GAL_E5B_FREQ)/(GAL_E1_FREQ-GAL_E5B_FREQ)
GAL_GF_E1E6_FREQ = (GAL_E1_FREQ*GAL_E6_FREQ)/(GAL_E1_FREQ-GAL_E6_FREQ)
GAL_GF_E1E5AltBOC_FREQ = (GAL_E1_FREQ*GAL_E5ALTBOC_FREQ)/(GAL_E1_FREQ-GAL_E5ALTBOC_FREQ)

