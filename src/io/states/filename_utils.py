""" Definition of output filenames
"""

# Output file names
AZEL = "satellite_azel.txt"
AMBIGUITY = "ambiguity.txt"
PHASE_BIAS = "receiver_phase_bias.txt"
DOP_ECEF = "DOP_ECEF.txt"
DOP_LOCAL = "DOP_ENU.txt"
POS_PRE_RESIDUALS = "pos_prefit_residuals.txt"
POS_POST_RESIDUALS = "pos_postfit_residuals.txt"
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
MW_OBSERVATIONS = "melbourne_wubbena_obs.txt"
GF_OBSERVATIONS = "geometry_free_obs.txt"

OUTPUT_FILENAME_MAP = {
    "satellite_azel": AZEL,
    "dop_ecef": DOP_ECEF,
    "dop_local": DOP_LOCAL,
    "pos_prefit_residuals": POS_PRE_RESIDUALS,
    "pos_postfit_residuals": POS_POST_RESIDUALS,
    "pr_prefit_residuals": POS_PRE_RESIDUALS,
    "pr_postfit_residuals": POS_POST_RESIDUALS,
    "cp_prefit_residuals": POS_PRE_RESIDUALS,
    "cp_postfit_residuals": POS_POST_RESIDUALS,
    "vel_prefit_residuals": VEL_PRE_RESIDUALS,
    "vel_postfit_residuals": VEL_POST_RESIDUALS,
    "pr_rate_prefit_residuals": VEL_PRE_RESIDUALS,
    "pr_rate_postfit_residuals": VEL_POST_RESIDUALS,
    "position": POSITION,
    "clock_bias": CLOCK,
    "time": TIME,
    "iono": IONO,
    "isb": ISB,
    "tropo_wet": TROPO,
    "velocity": VELOCITY,
    "clock_bias_rate": CLOCK_BIAS_RATE,
    "obs": OBSERVATIONS,
    "mw_obs": MW_OBSERVATIONS,
    "gf_obs": GF_OBSERVATIONS,
    "ambiguity": AMBIGUITY,
    "phase_bias": PHASE_BIAS
}
