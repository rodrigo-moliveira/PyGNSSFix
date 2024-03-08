"""Common errors declarations
"""
MODULE = "GNSSNavPy"

# TODO: simplificar um pouco os erros. Criar erros mais genéricos...

class GNSSNavPyError(Exception):
    """Generic error"""
    __module__ = MODULE


class UnknownScaleError(GNSSNavPyError):
    """Unknown timescale selected by user"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "UnknownScaleError Error -> " + message
        super().__init__(message)


class DateError(GNSSNavPyError):
    """Error when converting between timescales (no connection defined)"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Date Error -> " + message
        super().__init__(message)


class ConfigError(GNSSNavPyError):
    """Error when parsing the json configuration file"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Configuration Error -> " + message
        super().__init__(message)


class ConfigTypeError(GNSSNavPyError):
    """Error when parsing the json configuration file"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Configuration Type Error -> " + message
        super().__init__(message)


class UnknownConstellationError(GNSSNavPyError):
    """Unknown constellation defined by user"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Unknown Constellation Error -> " + message
        super().__init__(message)


class EopError(GNSSNavPyError):
    """Earth Orientation Parameters error (lack of data)"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "EOP Error -> " + message
        super().__init__(message)


class DataTypeError(GNSSNavPyError):
    """Unknown GNSS data type """
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "GNSS Data Type Error -> " + message
        super().__init__(message)


class TimeSeriesError(GNSSNavPyError):
    """Error in time series class"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Time Series Error -> " + message
        super().__init__(message)


class FileError(GNSSNavPyError):
    """Error reading input file"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "File Error -> " + message
        super().__init__(message)


class UnknownConversionError(GNSSNavPyError):
    """Error in time series class"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Unknown Conversions Error -> " + message
        super().__init__(message)


class NonExistentObservable(GNSSNavPyError):
    """Observable not found"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Non Existent Observable Error -> " + message
        super().__init__(message)


class DuplicateObservable(GNSSNavPyError):
    """Observable not found"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Duplicate Observable Error -> " + message
        super().__init__(message)


class EmptyObservationData(GNSSNavPyError):
    """Observation Data is empty"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Empty Observation Data Error -> " + message
        super().__init__(message)


class ArraySizeError(GNSSNavPyError):
    """Array has wrong size/shape"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Array Size Error -> " + message
        super().__init__(message)


class PreprocessorError(GNSSNavPyError):
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Preprocessor Error -> " + message
        super().__init__(message)


class NavigationError(GNSSNavPyError):
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Navigation Error -> " + message
        super().__init__(message)


class PVTComputationFail(GNSSNavPyError):
    """Error in PVT Computation"""
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Error in PVT Computation -> " + message
        super().__init__(message)


class EnumError(GNSSNavPyError):
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Enum Error -> " + message
        super().__init__(message)


class UnknownModel(GNSSNavPyError):
    __module__ = MODULE

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Unknown Model Error -> " + message
        super().__init__(message)
