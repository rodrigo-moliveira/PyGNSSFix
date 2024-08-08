"""Definition of output filenames
"""

# Output file names
AZEL = "satellite_azel.txt"
DOP_ECEF = "DOP_ECEF.txt"
DOP_LOCAL = "DOP_ENU.txt"
PRE_RESIDUALS = "prefit_residuals.txt"
POST_RESIDUALS = "postfit_residuals.txt"
VEL_PRE_RESIDUALS = "vel_prefit_residuals.txt"
VEL_POST_RESIDUALS = "vel_postfit_residuals.txt"
POSITION = "position.txt"
VELOCITY = "velocity.txt"
CLOCK = "clock_bias.txt"
TIME = "time.txt"
IONO = "iono.txt"
ISB = "isb.txt"
TROPO = "tropo.txt"
CLOCK_BIAS_RATE = "clock_bias_rate.txt"
OBSERVATIONS = "observations.txt"

OUTPUT_FILENAME_MAP = {
    "satellite_azel": AZEL,
    "dop_ecef": DOP_ECEF,
    "dop_local": DOP_LOCAL,
    "prefit_residuals": PRE_RESIDUALS,
    "postfit_residuals": POST_RESIDUALS,
    "vel_prefit_residuals": VEL_PRE_RESIDUALS,
    "vel_postfit_residuals": VEL_POST_RESIDUALS,
    "position": POSITION,
    "clock_bias": CLOCK,
    "time": TIME,
    "iono": IONO,
    "isb": ISB,
    "tropo_wet": TROPO,
    "velocity": VELOCITY,
    "clock_bias_rate": CLOCK_BIAS_RATE,
    "obs": OBSERVATIONS
}
