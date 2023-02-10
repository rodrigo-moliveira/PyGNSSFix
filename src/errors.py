"""Common errors declarations
"""


class Error(Exception):
    """Generic error"""
    pass


class DateError(Error):
    pass


class ConfigError(Error):
    pass


class EopError(Error):
    """Earth Orientation Parameters error (lack of data)"""
    pass
