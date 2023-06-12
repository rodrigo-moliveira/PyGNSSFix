"""Common errors declarations
"""
MODULE = "PyGNSSFix"


class PyGNSSFixError(Exception):
    """Generic error"""
    __module__ = MODULE

    # def __init__(self, message):
    #    # Call the base class constructor with the parameters it needs
    #    message = "Generic Error -> " + message
    #    super().__init__(message)


class UnknownScaleError(PyGNSSFixError):
    """Unknown timescale selected by user"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "UnknownScaleError Error -> " + message
        super().__init__(message)


class DateError(PyGNSSFixError):
    """Error when converting between timescales (no connection defined)"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Date Error -> " + message
        super().__init__(message)


class ConfigError(PyGNSSFixError):
    """Error when parsing the json configuration file"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Configuration Error -> " + message
        super().__init__(message)


class ConfigTypeError(PyGNSSFixError):
    """Error when parsing the json configuration file"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Configuration Type Error -> " + message
        super().__init__(message)


class UnknownConstellationError(PyGNSSFixError):
    """Unknown constellation defined by user"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Unknown Constellation Error -> " + message
        super().__init__(message)


class EopError(PyGNSSFixError):
    """Earth Orientation Parameters error (lack of data)"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "EOP Error -> " + message
        super().__init__(message)


class DataTypeError(PyGNSSFixError):
    """Unknown GNSS data type """
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "GNSS Data Type Error -> " + message
        super().__init__(message)


class TimeSeriesError(PyGNSSFixError):
    """Error in time series class"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Time Series Error -> " + message
        super().__init__(message)


class FileError(PyGNSSFixError):
    """Error reading input file"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "File Error -> " + message
        super().__init__(message)


class UnknownConversionError(PyGNSSFixError):
    """Error in time series class"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Unknown Conversions Error -> " + message
        super().__init__(message)


class NonExistentObservable(PyGNSSFixError):
    """Observable not found"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Non Existent Observable Error -> " + message
        super().__init__(message)


class DuplicateObservable(PyGNSSFixError):
    """Observable not found"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Duplicate Observable Error -> " + message
        super().__init__(message)


class EmptyObservationData(PyGNSSFixError):
    """Observation Data is empty"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Empty Observation Data Error -> " + message
        super().__init__(message)


class ArraySizeError(PyGNSSFixError):
    """Array has wrong size/shape"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Array Size Error -> " + message
        super().__init__(message)


class PreprocessorError(PyGNSSFixError):
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Preprocessor Error -> " + message
        super().__init__(message)
