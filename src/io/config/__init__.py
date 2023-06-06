"""import sys

this = sys.modules[__name__]

# we can explicitly make assignments on it
this.config = None


def initialize_config(config):
    if this.config is None:
        # also in local function scope. no scope specifier like global is needed
        this.config = config
    else:
        msg = "Database is already initialized to {0}."
        raise RuntimeError(msg.format(this.config))


def get_config():
    return this.config
"""
from .config import config_dict
