"""Common errors declarations
"""
MODULE = "PyGNSSFix"


class PyGNSSFixError(Exception):
    """Generic error"""
    __module__ = MODULE


class UnknownScaleError(PyGNSSFixError):
    """Unknown timescale selected by user"""
    __module__ = MODULE


class DateError(PyGNSSFixError):
    """Error when converting between timescales (no connection defined)"""
    __module__ = MODULE


class ConfigError(PyGNSSFixError):
    """Error when parsing the json configuration file"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Error reading the configuration file -> " + message
        super().__init__(message)


class ConfigTypeError(TypeError):
    """Error when parsing the json configuration file"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Error reading the configuration file -> " + message
        super().__init__(message)


class UnknownConstellationError(PyGNSSFixError):
    """Unknown constellation defined by user"""
    __module__ = MODULE


class EopError(PyGNSSFixError):
    """Earth Orientation Parameters error (lack of data)"""
    __module__ = MODULE


class DataTypeError(PyGNSSFixError):
    """Unknown GNSS data type """
    __module__ = MODULE
