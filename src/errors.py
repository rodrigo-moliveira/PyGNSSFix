""" Common errors declarations of the library """
MODULE = "GNSSNavPy"


class MetaErrorClass(type):
    def __new__(mcs, name, bases, dct):
        dct['__module__'] = MODULE
        return super().__new__(mcs, name, bases, dct)


class GNSSNavPyError(Exception, metaclass=MetaErrorClass):
    """ Generic error. """


class TimeScaleError(GNSSNavPyError):
    """ Timescale configuration error. """

    def __init__(self, message):
        message = "TimeScaleError Error -> " + message
        super().__init__(message)


class EpochError(GNSSNavPyError):
    """ Error when dealing with Epoch objects. """

    def __init__(self, message):
        message = "Epoch Error -> " + message
        super().__init__(message)


class ConfigError(GNSSNavPyError):
    """ Error when parsing the json configuration file. """

    def __init__(self, message):
        message = "Configuration Error -> " + message
        super().__init__(message)


class EopError(GNSSNavPyError):
    """ Earth Orientation Parameters error (lack of data). """

    def __init__(self, message):
        message = "EOP Error -> " + message
        super().__init__(message)


class SignalError(GNSSNavPyError):
    """ Unknown GNSS data signal type. """

    def __init__(self, message):
        message = "GNSS Data Type Signal Error -> " + message
        super().__init__(message)


class TimeSeriesError(GNSSNavPyError):
    """ Error in time series class. """

    def __init__(self, message):
        message = "Time Series Error -> " + message
        super().__init__(message)


class FileError(GNSSNavPyError):
    """ Error in IO management. """

    def __init__(self, message):
        message = "File Error -> " + message
        super().__init__(message)


class ArraySizeError(GNSSNavPyError):
    """ Array has wrong size/shape. """

    def __init__(self, message):
        message = "Array Size Error -> " + message
        super().__init__(message)


class PreprocessorError(GNSSNavPyError):
    """ Error in Preprocessor Execution. """

    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        message = "Preprocessor Error -> " + message
        super().__init__(message)


class EphemerideError(GNSSNavPyError):
    """ Error in computation of ephemerides. """

    def __init__(self, message):
        message = "Ephemeride Error -> " + message
        super().__init__(message)


class SolverError(GNSSNavPyError):
    """ Error in the PVT Computation of the Solver. """

    def __init__(self, message):
        message = "Error in PVT Computation -> " + str(message)
        super().__init__(message)


class ReconstructionError(GNSSNavPyError):
    """ Error in the reconstruction of observables. """

    def __init__(self, message):
        message = "Error in ObservationReconstructor Computation -> " + message
        super().__init__(message)


class EnumError(GNSSNavPyError):
    """ Enumeration Error. """

    def __init__(self, message):
        message = "Enumeration Error -> " + message
        super().__init__(message)


class CSpiceWrapperError(GNSSNavPyError):
    """ Error in the CSpiceWrapper class. """

    def __init__(self, message):
        message = "CSpiceWrapper Error -> " + message
        super().__init__(message)


class ITRFError(GNSSNavPyError):
    """ Error in the `ITRF_Transformations` class. """

    def __init__(self, message):
        message = "ITRF Transformations Error -> " + message
        super().__init__(message)
