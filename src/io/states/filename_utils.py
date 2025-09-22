""" Definition of output filenames
"""

# Output file names
AZEL = "satellite_azel.txt"
AMBIGUITY = "ambiguity.txt"
PHASE_BIAS = "receiver_phase_bias.txt"
DOP_ECEF = "DOP_ECEF.txt"
DOP_LOCAL = "DOP_ENU.txt"
PRE_RESIDUALS = "prefit_residuals.txt"
POST_RESIDUALS = "postfit_residuals.txt"
POSITION = "position.txt"
VELOCITY = "velocity.txt"
CLOCK = "clock_bias.txt"
TIME = "time.txt"
IONO = "iono.txt"
ISB = "isb.txt"
TROPO = "tropo.txt"
CLOCK_BIAS_RATE = "clock_bias_rate.txt"
OBSERVATIONS = "observations.txt"
RAW_OBSERVATIONS = "raw_observations.txt"
MW_OBSERVATIONS = "melbourne_wubbena_obs.txt"
GF_OBSERVATIONS = "geometry_free_obs.txt"
CYCLE_SLIPS = "cycle_slips.txt"

OUTPUT_FILENAME_MAP = {
    "satellite_azel": AZEL,
    "dop_ecef": DOP_ECEF,
    "dop_local": DOP_LOCAL,
    "prefit_residuals": PRE_RESIDUALS,
    "postfit_residuals": POST_RESIDUALS,
    "pr_prefit_residuals": PRE_RESIDUALS,
    "pr_postfit_residuals": POST_RESIDUALS,
    "cp_prefit_residuals": PRE_RESIDUALS,
    "cp_postfit_residuals": POST_RESIDUALS,
    "pr_rate_prefit_residuals": PRE_RESIDUALS,
    "pr_rate_postfit_residuals": POST_RESIDUALS,
    "position": POSITION,
    "clock_bias": CLOCK,
    "time": TIME,
    "iono": IONO,
    "isb": ISB,
    "tropo_wet": TROPO,
    "velocity": VELOCITY,
    "clock_bias_rate": CLOCK_BIAS_RATE,
    "obs": OBSERVATIONS,
    "raw_obs": RAW_OBSERVATIONS,
    "mw_obs": MW_OBSERVATIONS,
    "gf_obs": GF_OBSERVATIONS,
    "ambiguity": AMBIGUITY,
    "phase_bias": PHASE_BIAS,
    "cycle_slips": CYCLE_SLIPS
}
