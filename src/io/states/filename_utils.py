""" Definition of output filenames
"""

# Output file names
AZEL = "satellite_azel.txt"
AMBIGUITY = "ambiguity.txt"
PHASE_BIAS = "receiver_phase_bias.txt"
DOP_ECEF = "DOP_ECEF.txt"
DOP_LOCAL = "DOP_ENU.txt"
PR_PRE_RESIDUALS = "pr_prefit_residuals.txt"
PR_POST_RESIDUALS = "pr_postfit_residuals.txt"
PR_RATE_PRE_RESIDUALS = "pr_rate_prefit_residuals.txt"
PR_RATE_POST_RESIDUALS = "pr_rate_postfit_residuals.txt"
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
    "pr_prefit_residuals": PR_PRE_RESIDUALS,
    "pr_postfit_residuals": PR_POST_RESIDUALS,
    "pr_rate_prefit_residuals": PR_RATE_PRE_RESIDUALS,
    "pr_rate_postfit_residuals": PR_RATE_POST_RESIDUALS,
    "position": POSITION,
    "clock_bias": CLOCK,
    "time": TIME,
    "iono": IONO,
    "isb": ISB,
    "tropo_wet": TROPO,
    "velocity": VELOCITY,
    "clock_bias_rate": CLOCK_BIAS_RATE,
    "obs": OBSERVATIONS,
    "ambiguity": AMBIGUITY,
    "phase_bias": PHASE_BIAS
}
