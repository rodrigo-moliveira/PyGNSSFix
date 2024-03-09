from .logger import get_logger

__all__ = ["get_logger", "set_logs"]


# initialize common application loggers
def set_logs(file_level, log_path):
    get_logger("IO", file_level, log_path)
    get_logger("PREPROCESSOR", file_level, log_path)
    get_logger("MAIN_LOG", file_level, log_path)
    get_logger("GNSS_SOLVER", file_level, log_path)
    get_logger("INS_ALG", file_level, log_path)
    get_logger("PERFORMANCE", file_level, log_path)
