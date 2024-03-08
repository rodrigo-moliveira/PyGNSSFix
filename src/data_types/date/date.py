"""Date module
"""

from numpy import sin, radians
from datetime import datetime, timedelta, date

from ...errors import EpochError, TimeScaleError
from .eop import EopDb
from ...utils.node import Node

__all__ = ["Epoch", "timedelta"]


"""TAI = UTC + LEAPS
TAI = GPS + 19
UTC + LEAPS = GPS + 19
UTC + 37 = GPS + 19
UTC + 18 = GPS
UTC = GPS -18"""


class Timescale(Node):
    """Definition of a timescale and its interactions with others"""

    def __repr__(self):  # pragma: no cover
        return f"<Scale '{self.name}'>"

    def __str__(self):
        return self.name

    def _scale_ut1_minus_utc(self, mjd, eop):
        """Definition of Universal Time relatively to Coordinated Universal Time"""
        return eop.ut1_utc

    def _scale_tai_minus_utc(self, mjd, eop):
        """Definition of International Atomic Time relatively to Coordinated Universal Time
        This is also known as the DAT (leap seconds)"""
        return eop.tai_utc

    def _scale_tt_minus_tai(self, mjd, eop):
        """Definition of Terrestrial Time relatively to International Atomic Time"""
        return 32.184

    def _scale_tai_minus_gpst(self, mjd, eop):
        """Definition of International Atomic Time relatively to GPS time"""
        return 19.0

    def _scale_gst_minus_gpst(self, mjd, eop):
        """Definition of GST relatively to GPST (GPST = GST + EOP.GGTO)"""
        # From theory GGTO = GST - GPST is GPS-to-Galileo Time Offset
        return eop.ggto

    def _scale_tdb_minus_tt(self, mjd, eop):
        """Definition of the Barycentric Dynamic Time scale relatively to Terrestrial Time"""
        # NOTE: This is the tdb_opt 3 from Vallado script (ast alm approach (2012) bradley email)
        jd = mjd + Epoch.JD_MJD
        jj = Epoch._julian_century(jd)
        m = radians(357.5277233 + 35999.05034 * jj)
        delta_lambda = radians(246.11 + 0.90251792 * (jd - Epoch.J2000))
        return 0.001657 * sin(m) + 0.000022 * sin(delta_lambda)

    def offset(self, mjd, new_scale, eop):
        """Compute the offset necessary in order to convert from one timescale to another

        Args:
            mjd (float):
            new_scale (str): Name of the desired scale
            eop (Eop): class with Eop data
        Return:
            float: offset to apply in seconds
        """

        delta = 0
        for one, two in self.steps(new_scale):
            one = one.name.lower()
            two = two.name.lower()
            # find the operation
            oper = f"_scale_{two}_minus_{one}"
            # find the reverse operation
            roper = f"_scale_{one}_minus_{two}"
            if hasattr(self, oper):
                delta += getattr(self, oper)(mjd, eop)
            elif hasattr(self, roper):
                delta -= getattr(self, roper)(mjd, eop)
            else:  # pragma: no cover
                raise EpochError(f"Unknown conversion {one} => {two}")

        return delta


UT1 = Timescale("UT1")      # Universal Time
GPST = Timescale("GPST")    # GPS Time
TDB = Timescale("TDB")      # Barycentric Dynamical Time
UTC = Timescale("UTC")      # Coordinated Universal Time
TAI = Timescale("TAI")      # International Atomic Time
TT = Timescale("TT")        # Terrestrial Time
GST = Timescale("GST")      # Galileo System Time

GPST + TAI + UTC + UT1
TDB + TT + TAI
GPST + GST

_cache = {"UT1": UT1, "GPST": GPST, "TDB": TDB, "UTC": UTC, "TAI": TAI, "TT": TT, "GPS": GPST,
          "GST": GST, "GAL": GST}


def get_scale(name):
    if name in _cache.keys():
        return _cache[name]
    else:
        raise TimeScaleError(name)


class Epoch:
    """Date object

    All computations and in-memory saving are made in
    `MJD <https://en.wikipedia.org/wiki/Julian_day>`__ and
    `UTC.
    In the current implementation, the Date object does not handle the
    leap second.

    The constructor can take:

        * the same arguments as the standard library's datetime object (year, month, day, hour,
          minute, second, microsecond)
        * MJD as :py:class:`float`
        * MJD as :py:class:`int` for days and :py:class:`float` for seconds
        * a :py:class:`Date` or :py:class:`datetime` object

    Keyword Arguments:
        scale (str) : One of the following scales : "UT1", "UTC", "GPS", "TDB", "TAI", "TT"

    Examples:

        .. code-block:: python

            Date(2016, 11, 17, 19, 16, 40)
            Date(2016, 11, 17, 19, 16, 40, scale="TAI")
            Date(57709.804455)  # MJD
            Date(57709, 69540.752649)
            Date(datetime(2016, 11, 17, 19, 16, 40))  # built-in datetime object
            Date.now()

    Date objects interact with :py:class:`timedelta` as datetime do.

    Attributes:
        eop (Eop): Value of the Earth Orientation Parameters for this particular date (see
            :ref:`eop`)
        scale: Scale in which this date is represented
    """

    __slots__ = ["_d", "_s", "_offset", "scale", "_cache", "eop"]

    MJD_T0 = datetime(1858, 11, 17)
    """Origin of MJD"""

    GPS_ORIGIN = datetime(1980, 1, 6, 0, 0, 0)
    """Origin of GPST (6th January 1980 00:00:00 midnight)"""

    JD_MJD = 2400000.5
    """Offset between JD and MJD"""

    J2000 = 2451545.0
    """Offset between JD and J2000"""

    REF_SCALE = "GPST"
    """Scale used internally to store the epoch"""

    DEFAULT_SCALE = "UTC"
    """Default scale"""

    DEFAULT_FORMAT = "%Y-%m-%d %H:%M:%S"

    def __init__(self, *args, scale=DEFAULT_SCALE, **kwargs):

        # TODO: constructor with seconds of week and week number???

        if type(scale) is str:
            scale = get_scale(scale.upper())
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, datetime):
                # Python datetime.datetime object
                d, s = self._convert_dt(arg)
            elif isinstance(arg, Epoch):
                # Date object
                d = arg.d
                s = arg.s
                scale = arg.scale
            elif isinstance(arg, (float, int)):
                # Modified Julian Day
                if isinstance(arg, int):
                    d = arg
                    s = 0.0
                else:
                    d = int(arg)
                    s = (arg - d) * 86400
            else:
                raise TypeError(f"Unknown type '{type(arg)}'")
        elif len(args) == 2 and (
                isinstance(args[0], int) and isinstance(args[1], (int, float))
        ):
            # Julian day and seconds in the day
            d, s = args
        elif len(args) in range(3, 8) and list(map(type, args)) == [int] * len(args):
            # Same constructor as datetime.datetime
            # (year, month, day, hour=0, minute=0, second=0, microsecond=0, tzinfo=None)
            dt = datetime(*args, **kwargs)
            d, s = self._convert_dt(dt)
        else:
            raise TypeError(
                f"Unknown type sequence {', '.join(str(type(x)) for x in args)}"
            )

        mjd = d + s / 86400.0

        # Retrieve EOP for the given date and store
        eop = EopDb.get(mjd)

        # Retrieve the offset from REF_SCALE for the current date
        offset = scale.offset(mjd, self.REF_SCALE, eop)

        d += int((s + offset) // 86400)
        s = (s + offset) % 86400.0

        # As Date acts like an immutable object, we can't set its attributes normally
        # like when we do ``self._d = _d``. Furthermore, those attribute represent the date with
        # respect to REF_SCALE
        super().__setattr__("_d", d)
        super().__setattr__("_s", s)
        super().__setattr__("_offset", offset)
        super().__setattr__("scale", scale)
        super().__setattr__("eop", eop)
        super().__setattr__("_cache", {})

    def __getstate__(self):  # pragma: no cover
        """Used for pickling"""
        return {
            "d": self._d,
            "s": self._s,
            "offset": self._offset,
            "scale": self.scale,
            "eop": self.eop,
        }

    def __setstate__(self, state):  # pragma: no cover
        """Used for unpickling"""
        super().__setattr__("_d", state["d"])
        super().__setattr__("_s", state["s"])
        super().__setattr__("_offset", state["offset"])
        super().__setattr__("scale", state["scale"])
        super().__setattr__("eop", state["eop"])
        super().__setattr__("_cache", {})

    def __setattr__(self, *args):  # pragma: no cover
        raise TypeError("Cannot modify attributes of immutable object")

    def __delattr__(self, *args):  # pragma: no cover
        raise TypeError("Cannot modify attributes of immutable object")

    def __add__(self, other):
        if isinstance(other, timedelta):
            days, sec = divmod(other.total_seconds() + self.s, 86400)
        else:
            raise TypeError(f"Unknown operation with {type(other)}")

        return self.__class__(self.d + int(days), sec, scale=self.scale)

    def __sub__(self, other):
        if isinstance(other, timedelta):
            other = timedelta(seconds=-other.total_seconds())
        elif isinstance(other, datetime):
            return self.datetime - other
        elif isinstance(other, Epoch):
            return self._datetime - other._datetime
        else:
            raise TypeError(f"Unknown operation with {type(other)}")

        return self.__add__(other)

    def __gt__(self, other):
        return self._mjd > other._mjd

    def __ge__(self, other):
        return self._mjd >= other._mjd

    def __lt__(self, other):
        return self._mjd < other._mjd

    def __le__(self, other):
        return self._mjd <= other._mjd

    def __eq__(self, other):
        return self._mjd == other._mjd

    def __repr__(self):  # pragma: no cover
        return f"<{self.__class__.__name__} '{self}'>"

    def debug_print(self, scale="ref_scale"):
        if scale == "ref_scale":
            return f"({self._d},{self._s}) {self.REF_SCALE}"
        else:
            return f"{self._convert_to_scale()} {self.scale}"

    def __str__(self):  # pragma: no cover
        if "str" not in self._cache.keys():
            self._cache["str"] = f"{self.datetime.isoformat()} {self.scale}"
        return self._cache["str"]

    def __format__(self, fmt):  # pragma: no cover
        if fmt:
            return self.datetime.__format__(fmt)
        else:
            return str(self)

    def __hash__(self):
        return hash((self._d, self._s))

    # TODO: add just a function to convert between GPST and GST

    @classmethod
    def _convert_dt(cls, dt):
        if dt.tzinfo is None:
            delta = dt - cls.MJD_T0
        else:
            tz = dt.utcoffset()
            delta = dt.replace(tzinfo=None) - cls.MJD_T0 - tz

        return delta.days, delta.seconds + delta.microseconds * 1e-6

    def _convert_to_scale(self):
        """Convert the inner value (defined with respect to REF_SCALE) into the given scale
        of the object
        """
        d = self._d
        s = (self._s - self._offset) % 86400.0
        d -= int((s + self._offset) // 86400)
        return d, s

    @property
    def d(self):
        return self._convert_to_scale()[0]

    @property
    def s(self):
        return self._convert_to_scale()[1]

    @property
    def datetime(self):
        """Conversion of the Date object into a ``datetime.datetime``

        The resulting object is a timezone-naive instance with the same scale
        as the originating Date object.
        """
        if "dt_scale" not in self._cache.keys():
            self._cache["dt_scale"] = self._datetime - timedelta(seconds=self._offset)
        return self._cache["dt_scale"]

    @property
    def _datetime(self):
        """Conversion of the Date object into a :py:class:`datetime.datetime`.

        The resulting object is a timezone-naive instance in the REF_SCALE timescale
        """
        if "dt" not in self._cache.keys():
            self._cache["dt"] = self.MJD_T0 + timedelta(days=self._d, seconds=self._s)
        return self._cache["dt"]

    @classmethod
    def strptime(cls, data, format=DEFAULT_FORMAT, scale=DEFAULT_SCALE):  # pragma: no cover
        """Convert a string representation of a date to a Date object"""
        return cls(datetime.strptime(data, format), scale=scale)

    @classmethod
    def now(cls, scale=DEFAULT_SCALE):
        """
        Args:
            scale (str)
        Return:
            Date: Current time in the chosen scale
        """
        return cls(datetime.utcnow()).change_scale(scale)

    def strftime(self, fmt):  # pragma: no cover
        """Format the date following the given format"""
        return self.datetime.strftime(fmt)

    def change_scale(self, new_scale):
        """
        Args:
            new_scale (str)
        Return:
            Date: new Date object representing the same instant, with a different timescale
        """
        offset = self.scale.offset(self._mjd, new_scale, self.eop)
        result = self.datetime + timedelta(seconds=offset)

        return self.__class__(result, scale=new_scale)

    @classmethod
    def _julian_century(cls, jd):
        return (jd - cls.J2000) / 36525.0

    @property
    def julian_century(self):
        """Compute the julian_century of the Date object relatively to its
        scale

        Return:
            float
        """
        return self._julian_century(self.jd)

    @property
    def jd(self):
        """Compute the Julian Date, which is the number of days from the
        January 1, 4712 B.C., 12:00.

        Return:
            float
        """
        return self.mjd + self.JD_MJD

    @property
    def _mjd(self):
        """
        Return:
            float: Date in terms of MJD in the REF_SCALE timescale (internal timescale)
        """
        return self._d + self._s / 86400.0

    @property
    def mjd(self):
        """Date in terms of MJD (in scale defined by the Epoch)

        Return:
            float
        """
        return self.d + self.s / 86400.0

    @property
    def gnss_time(self):
        """
        Computes GPS Time GPST (Week number and Seconds of Week) for this Epoch
        The origin of GPST (week = 0, seconds of week = 0) is 6 January 1980 at 00:00 midnight

        In RINEX files the GPS week number is provided as a continuous number, without roll-over

        NOTE: This is also valid for GST (Galileo System Time). According to RINEX documentation,
        The GAL week number is a continuous number, aligned to (and hence identical to) the
        continuous GPS week number used in the RINEX navigation message files
        """
        if "gnss_time" not in self._cache.keys():
            # convert this epoch to GPS Time
            gps_epoch = self.change_scale(GPST)
            dt = (gps_epoch - Epoch.GPS_ORIGIN).total_seconds()
            week = (dt/3600/24) // 7
            sow = dt - week*3600*24*7
            self._cache["gnss_time"] = (int(week), sow)
        return self._cache["gnss_time"]

    @property
    def doy(self):
        """Computes day of year"""
        if "doy" not in self._cache.keys():
            self._cache["doy"] = self.datetime.timetuple().tm_yday
        return self._cache["doy"]


# This part is here to allow matplotlib to display Date objects directly
# in the plot, without any other conversion by the developer
# If matplotlib is importable, then a converter class is registered
# for converting all Date objects on the fly
try:
    import matplotlib.dates as mdates
    import matplotlib.units as munits
except ImportError:  # pragma: no cover
    pass
else:  # pragma: no cover

    class DateConverter(mdates.DateConverter):
        @staticmethod
        def _conv(v):
            if isinstance(v, (datetime, date)):
                v = mdates.date2num(v)
            else:
                v = mdates.date2num(v.datetime)
            return v

        @staticmethod
        def convert(values, unit, axis):

            try:
                iter(values)
            except TypeError:
                values = [values]

            values = [DateConverter._conv(v) for v in values]

            return values


    munits.registry.setdefault(Epoch, DateConverter())
