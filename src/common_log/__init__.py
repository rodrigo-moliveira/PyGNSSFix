""" The `common_log` package implements an interface for the
`Logging Python package <https://docs.python.org/3/library/logging.html>`__.

All log objects required by this project are set and initialized here. To get the logs, the user of this package
first has to call function `set_logs` below (to set the log level and path).
Then, the user modules and functions import the loggers and simply use them.

The names of the loggers represent different high-level procedures in the PyGNSSFix project:
    * MAIN_LOG: main logger of the program
    * IO_LOG: logger to be used in input/output operations
    * PREPROCESSOR_LOG: logger to be used in data preprocessing operations
    * MODEL_LOG: logger to be used in the different models (physical models, attitude functions, etc.)
    * GNSS_ALG_LOG: logger of the GNSS Solver algorithm
    * INS_ALG_LOG: logger of the INS algorithms (sensor emulator, INS propagation, navigation filter)
    * POST_PROC_LOG: logger of the Post-Processing / Performance Evaluation algorithm

Examples:
    This is an example of how to set up and use the logger objects:

    >>> set_logs("DEBUG", "sample_log.txt")
    >>> io_log = get_logger(IO_LOG)
    >>> io_log.info("Information Message")
    >>> io_log.error("Error Message")
"""
from .logger import get_logger

__all__ = ["get_logger", "set_logs", "IO_LOG", "PREPROCESSOR_LOG", "MAIN_LOG", "GNSS_ALG_LOG", "INS_ALG_LOG",
           "POST_PROC_LOG", "MODEL_LOG"]


MODEL_LOG = str("MODEL_LOG")
IO_LOG = str("IO_LOG")
PREPROCESSOR_LOG = str("PREPROCESSOR_LOG")
MAIN_LOG = str("MAIN_LOG")
GNSS_ALG_LOG = str("GNSS_ALG_LOG")
INS_ALG_LOG = str("INS_ALG_LOG")
POST_PROC_LOG = str("POST_PROC_LOG")


def set_logs(severity_level, log_path=""):
    """ Initializes the logger objects for the PyGNSSFix program

    Args:
        severity_level(str): severity level of the log. Available levels are 'CRITICAL', 'ERROR', 'WARNING', 'INFO',
            'DEBUG', 'NOTSET'
        log_path(str): path to save the log file
    """
    for log in (IO_LOG, PREPROCESSOR_LOG, MAIN_LOG, GNSS_ALG_LOG, INS_ALG_LOG, POST_PROC_LOG, MODEL_LOG):
        get_logger(log, severity_level, log_path)
