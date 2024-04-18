from .logger import get_logger

__all__ = ["get_logger", "set_logs", "IO_LOG", "PREPROCESSOR_LOG", "MAIN_LOG", "GNSS_ALG_LOG", "INS_ALG_LOG",
           "PERFORMANCE_LOG", "MODEL_LOG"]

""" All log objects required by this package are set and initialized here. To get the logs, the user first 
has to call function `set_logs` below (to set the log level and path). 
Then, in the user modules and functions, you import the logger and simply use it.

Example:
io_log = get_logger(IO_LOG)
io_log.info("Info")
io_log.error("Error")
"""

MODEL_LOG = str("MODEL_LOG")
IO_LOG = str("IO_LOG")
PREPROCESSOR_LOG = str("PREPROCESSOR_LOG")
MAIN_LOG = str("MAIN_LOG")
GNSS_ALG_LOG = str("GNSS_ALG_LOG")
INS_ALG_LOG = str("INS_ALG_LOG")
PERFORMANCE_LOG = str("PERFORMANCE_LOG")


# initialize common application loggers
def set_logs(file_level, log_path=""):
    for log in (IO_LOG, PREPROCESSOR_LOG, MAIN_LOG, GNSS_ALG_LOG, INS_ALG_LOG, PERFORMANCE_LOG, MODEL_LOG):
        get_logger(log, file_level, log_path)

