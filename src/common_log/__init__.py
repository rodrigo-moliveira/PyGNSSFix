from .logger import get_logger

__all__ = ["get_logger", "set_logs", "IO_LOG", "PREPROCESSOR_LOG", "MAIN_LOG", "GNSS_ALG_LOG", "INS_ALG_LOG",
           "PERFORMANCE_LOG"]

IO_LOG = str("IO_LOG")
PREPROCESSOR_LOG = str("PREPROCESSOR_LOG")
MAIN_LOG = str("MAIN_LOG")
GNSS_ALG_LOG = str("GNSS_ALG_LOG")
INS_ALG_LOG = str("INS_ALG_LOG")
PERFORMANCE_LOG = str("PERFORMANCE_LOG")


# initialize common application loggers
def set_logs(file_level, log_path):
    for log in (IO_LOG, PREPROCESSOR_LOG, MAIN_LOG, GNSS_ALG_LOG, INS_ALG_LOG, PERFORMANCE_LOG):
        get_logger(log, file_level, log_path)

