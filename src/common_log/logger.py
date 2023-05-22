import logging

"""
from ``logging`` docs:
    * Logger: This is the class whose objects will be used in the application code directly to call the functions.
    * LogRecord: Loggers automatically create LogRecord objects that have all the information related to the event being
        logged, like the name of the logger, the function, the line number, the message, and more.
    * Handler: Handlers send the LogRecord to the required output destination, like the console or a file. Handler is a
        base for subclasses like StreamHandler, FileHandler, SMTPHandler, HTTPHandler, and more. These subclasses send
        the logging outputs to corresponding destinations, like sys.stdout or a disk file.
    * Formatter: This is where you specify the format of the output by specifying a string format that lists out the
        attributes that the output should contain.
"""

__logs__ = {}  # global dict for all logs
__file_level__ = logging.INFO
_file_path = ""


def clean_logs():
    __logs__.clear()


def set_log(log_str, logger):
    if log_str not in __logs__:
        __logs__[log_str] = logger


def get_logger(log_str, file_level=logging.INFO, file_path=""):
    """Returns a logger. This logger writes both to the consoled and to a file
    The severity level of the console is ERROR. For the file it is defined by the user

    Args:
        log_str (str): logger name string
        file_level(logging.LEVEL): severity level of log events to be printed in the file
        file_path (str): path to the log file
    Return:
        logger (logging)
    """
    if log_str in __logs__:
        return __logs__[log_str]
    else:
        return setup(log_str, file_level, file_path)


def setup(log_str, file_level, file_path):
    if log_str in __logs__:
        return __logs__[log_str]

    # Create custom loggers for provided modules (in ´´list_of_loggers´´)
    new_log = logging.getLogger(log_str)
    new_log.handlers = []  # remove all handlers
    new_log.setLevel(logging.DEBUG)  # default logging level DEBUG

    # add new_log to __logs__
    set_log(log_str, new_log)

    # Create handlers
    c_handler = logging.StreamHandler()
    f_handler = logging.FileHandler(file_path, mode='a')
    c_handler.setLevel(logging.WARN)  # console log level
    f_handler.setLevel(file_level)  # file log level

    # Create formatters and add it to handlers
    c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    f_format = logging.Formatter('[%(asctime)s] :: [%(name)s --- %(levelname)s] :: %(message)s')

    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)

    new_log.addHandler(c_handler)
    new_log.addHandler(f_handler)

    return new_log
