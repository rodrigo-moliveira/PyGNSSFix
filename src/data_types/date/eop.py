"""Retrieve and interpolate data for Earth Orientation and timescales conversions
"""

from pathlib import Path
from inspect import isclass

from src.errors import EopError, ConfigError, FileError

__all__ = ["register", "EopDb", "TaiUtc", "Finals", "Finals2000A", "GGTO"]


class TaiUtc:
    """File listing all leap seconds throughout history

    This file could be retrieved `here <http://maia.usno.navy.mil/ser7/tai-utc.dat>`
    """

    def __init__(self, path):

        self.path = Path(path)
        self.data = []

        try:
            f_handler = open(self.path, 'r')
        except OSError:
            raise FileError(f"Could not open/read file: {self.path}", )

        with f_handler:
            lines = f_handler.read().splitlines()
        for line in lines:
            if not line:
                continue

            line = line.split()
            mjd = int(float(line[4]) - 2400000.5)
            value = float(line[6])
            self.data.append((mjd, value))

    def __getitem__(self, date):
        for mjd, value in reversed(self.data):
            if mjd <= date:
                return value

    def get_last_next(self, date):
        """Provide the last and next leap-second events relative to a date

        Args:
            date (float): Date in MJD
        Return:
            tuple:
        """
        past, future = (None, None), (None, None)

        for mjd, value in reversed(self.data):
            if mjd <= date:
                past = (mjd, value)
                break
            future = (mjd, value)

        return past, future

    def set_leap_seconds_rinex(self):
        pass
        # TODO: add this


class GGTO:
    def __init__(self):
        self.data = list()
        self.global_ggto = 0

    def set_ggto(self, date, ggto):
        for index, (mjd, value) in enumerate(reversed(self.data)):
            i = len(self.data) - index
            if date >= mjd:
                self.data.insert(i, (date, ggto))
                return
        self.data.insert(0, (date, ggto))

    def set_global_ggto(self, global_ggto):
        self.global_ggto = global_ggto

    def __getitem__(self, date):
        for mjd, value in reversed(self.data):
            if mjd <= date:
                return value
        return 0.0

    def print(self):
        for mjd, value in self.data:
            print(mjd, value)


class Finals2000A:
    """History of Earth orientation correction for IAU2000 model

    Three files are available `here <https://datacenter.iers.org/eop.php>`__ for this model:

        - **finals2000A.all**, from 1976-01-02 to present, updated weekly
        - **finals2000A.data**, from 1992-01-01 to present, updated weekly
        - **finals2000A.daily**, last 90 days + 90 days of prediction, updated daily

    See the associated metadata for more information about the content of these files.
    """

    deltas = ("dx", "dy")

    def __init__(self, path):

        self.path = Path(path)
        d1, d2 = self.deltas

        try:
            f_handler = open(self.path, 'r')
        except OSError:
            raise FileError(f"Could not open/read file: {self.path}", )

        with f_handler:
            lines = f_handler.read().splitlines()

        self.data = {}
        for line in lines:
            line = line.rstrip()
            mjd = int(float(line[7:15]))

            try:
                self.data[mjd] = {
                    "mjd": mjd,
                    # 'flag': line[16],
                    "x": float(line[18:27]),
                    d1: None,
                    # 'Xerror': float(line[27:36]),
                    "y": float(line[37:46]),
                    d2: None,
                    # 'Yerror': float(line[46:55]),
                    "lod": None,
                    "ut1_utc": float(line[58:68]),
                }
            except ValueError:
                # Common values (X, Y, UT1-UTC) are not available anymore
                break
            else:
                try:
                    self.data[mjd][d1] = float(line[97:106])
                    self.data[mjd][d2] = float(line[116:125])
                except ValueError:
                    # dX and dY are not available for this date, so we take
                    # the last value available
                    self.data[mjd][d1] = self.data[mjd - 1][d1]
                    self.data[mjd][d2] = self.data[mjd - 1][d2]
                    pass
                try:
                    self.data[mjd]["lod"] = float(line[79:86])
                except ValueError:
                    # LOD is not available for this date, so we take the last value available
                    self.data[mjd]["lod"] = self.data[mjd - 1]["lod"]
                    pass

    def __getitem__(self, key):
        return self.data[key]

    def items(self):
        return self.data.items()


class Finals(Finals2000A):
    """History of Earth orientation correction for IAU1980 model

    Three files are available `here <https://datacenter.iers.org/eop.php>`__ for this model:

        - **finals.all**, from 1976-01-02 to present, updated weekly
        - **finals.data**, from 1992-01-01 to present, updated weekly
        - **finals.daily**, last 90 days + 90 days of prediction, updated daily

    See the associated metadata for more information about the content of these files.
    """

    deltas = ("dpsi", "deps")


class Eop:
    """Earth Orientation Parameters"""

    def __init__(self, **kwargs):
        self.x = kwargs.get("x", 0)
        self.y = kwargs.get("y", 0)
        self.dx = kwargs.get("dx", 0)
        self.dy = kwargs.get("dy", 0)
        self.deps = kwargs.get("deps", 0)
        self.dpsi = kwargs.get("dpsi", 0)
        self.lod = kwargs.get("lod", 0)
        self.ut1_utc = kwargs.get("ut1_utc", 0)
        self.tai_utc = kwargs.get("tai_utc", 0)
        self.ggto = kwargs.get("ggto", 0)

    def __repr__(self):
        return "{name}(x={x}, y={y}, dx={dx}, dy={dy}, deps={deps}, dpsi={dpsi}, lod={lod}, ut1_utc={ut1_utc}, " \
               "tai_utc={tai_utc}, ggto={ggto})".format(name=self.__class__.__name__, **self.__dict__)


class EopDb:
    """Class handling the different EOP databases available, in a simple abstraction layer.

    By defining a simple parameter in the config dict, this class will handle the instantiation
    of the database and queries in a transparent manner.

    see :ref:`dbname <eop-dbname>` and :ref:`missing policy <eop-missing-policy>` configurations.
    """

    _dbs = {}
    DEFAULT_DBNAME = "default"
    """Default name used for EOP database lookup."""

    PASS = "pass"
    WARN = "warning"
    ERROR = "error"

    MIS_DEFAULT = WARN
    """Default behaviour in case of missing value"""

    @classmethod
    def db(cls, dbname=DEFAULT_DBNAME):
        """Retrieve the database

        Args:
            dbname: Specify the name of the database to retrieve. If set to `None`, take the name
                from the configuration (see :ref:`configuration <eop-dbname>`)
        Return:
            object
        """

        if dbname not in cls._dbs.keys():
            raise EopError(f"Unknown database '{dbname}'")

        if isclass(cls._dbs[dbname]):
            # Instantiation
            try:
                cls._dbs[dbname] = cls._dbs[dbname]()
            except Exception as e:
                # Keep the exception in cache in order to not retry instantiation
                # every single time EopDb.db() is called, as instantiation
                # of database is generally a time-consuming operation.
                # If it failed once, it will most probably fail again
                cls._dbs[dbname] = e

        if isinstance(cls._dbs[dbname], Exception):
            raise cls._dbs[dbname]

        return cls._dbs[dbname]

    @classmethod
    def get(cls, mjd: float, dbname: str = DEFAULT_DBNAME) -> Eop:
        """Retrieve Earth Orientation Parameters and timescales differences
        for a given date

        Args:
            mjd: Date expressed as Modified Julian Date
            dbname: Name of the database to use
        Return:
            Eop: Interpolated data for this particular MJD
        """
        try:
            value = cls.db(dbname)[mjd]
        except KeyError as e:
            msg = f"Missing EOP data for mjd = '{e}'"
            if cls.policy() == cls.WARN:
                print(msg)  # TODO print warning log.warning(msg)
            elif cls.policy() == cls.ERROR:
                raise e

            value = Eop(
                x=0, y=0, dx=0, dy=0, deps=0, dpsi=0, lod=0, ut1_utc=0, tai_utc=0, ggto=0
            )

        return value

    @classmethod
    def policy(cls):
        pol = cls.MIS_DEFAULT  # TODO add policy to config config.get("eop", "missing_policy", fallback=cls.MIS_DEFAULT)
        if pol not in (cls.PASS, cls.WARN, cls.ERROR):
            raise ConfigError("Unknown config value for 'eop.missing_policy'")

        return pol

    @classmethod
    def register(cls, klass, name=DEFAULT_DBNAME):
        """Register an Eop Database

        The only requirement of this database is that it should have ``__getitem__``
        method accepting MJD as float.
        """
        if name in cls._dbs:
            msg = f"'{name}' is already registered for an Eop database. Skipping"
            print(msg)  # TODO: raise warning
        else:
            cls._dbs[name] = klass


def register(name=EopDb.DEFAULT_DBNAME):
    """Decorator for registering an Eop Database

    Example:

    .. code-block:: python

        @register
        class SqliteEnvDatabase:
            # sqlite implementation
            # this database will be known as 'default'

        @register('json')
        class JsonEnvDatabase:
            # JSON implementation

        EopDb.get(58090.2)                    # get Eop from SqliteEnvDatabase
        EopDb.get(58090.2, dbname='default')  # same as above
        EopDb.get(58090.2, dbname='json')     # get Eop from JsonEnvDatabase
    """

    # I had a little trouble setting this function up, due to the fact that
    # I wanted it to be usable both as a simple decorator (``@register``)
    # and a decorator with arguments (``@register('mydatabase')``).
    # The current implementation allows this dual-use, but it's a bit hacky.

    # In the simple decorator mode, when the @register decorator is called
    # the argument passed is the class to decorate. So it *is* the decorated
    # function

    # In the decorator-with-arguments mode, the @register decorator should provide
    # a callable that will be the decorated function. This callable takes
    # the class you want to decorate

    if isinstance(name, str):
        # decorator with argument
        def wrapper(klass):
            EopDb.register(klass, name)
            return klass

        return wrapper

    else:
        # simple decorator mode
        klass = name

        EopDb.register(klass)
        return klass


@register
class SimpleEopDatabase:
    """Simple implementation of database

    Uses ``tai-utc.dat``, ``finals.all`` and ``finals2000A.all`` files directly
    without caching nor interpolation.

    In order to use these files, you have to provide the directory containing them as a config
    variable. Optionally, you can provide the type of data you want to extract from finals files
    ('all', 'data' or 'daily').

    .. code-block:: python

        from beyond.config import config
        config.update({
            'eop': {
                'folder': "/path/to/eop/data/",
                'type': "all"
            }
        })
    """

    def __init__(self):
        from src import WORKSPACE_PATH
        from src.io.config import config_dict

        leap_file = config_dict.get("inputs", "leap_file")
        finals_file = config_dict.get("inputs", "finals_file")

        # Data reading
        f = Finals(WORKSPACE_PATH / f"{finals_file}")
        # f2 = Finals2000A(WORKSPACE_PATH / f"finals2000A")
        t = TaiUtc(WORKSPACE_PATH / f"{leap_file}")

        # Extracting data from finals files
        self._finals = {}
        for date, values in f.items():
            self._finals[date] = values
            # self._finals[date].update(f2[date])

        self._tai_utc = t.data.copy()
        self._ggto = GGTO()

    def __getitem__(self, mjd):
        data = self.finals(mjd)
        data["tai_utc"] = self.tai_utc(mjd)
        data["ggto"] = self._ggto[mjd]

        return Eop(**data)

    def set_ggto(self, mjd: float, ggto):
        self._ggto.set_ggto(mjd, ggto)

    def set_global_ggto(self, global_ggto):
        self._ggto.set_global_ggto(global_ggto)

    def finals(self, mjd: float):
        return self._finals[int(mjd)].copy()

    def tai_utc(self, mjd: float):
        for date, value in reversed(self._tai_utc):
            if date <= mjd:
                return value
        else:
            raise KeyError(mjd)
