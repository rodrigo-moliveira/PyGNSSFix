import math

from src.data_types.date import Epoch
from src.constants import SECONDS_IN_DAY, DS2R, PI, DAS2R

__all__ = ["get_gmst", "iauGMST06", "iauGMST82", "iauGMST00"]


def anp(theta):
    """ Normalize angle into the range [0, 2π). """
    return theta % (2*PI)


def get_gmst(epoch, model):
    """
    Compute Greenwich Mean Sidereal Time (GMST).

    Depending on the selected input parameter, three models are available:
        * IAU1982 (`model=IAU82`)
        * IAU2000 (`model=IAU00`)
        * IAU2006 (`model=IAU06`)

    Args:
        epoch(Epoch): input epoch instance
        model(str): model to compute GMST

    Raises:
        ValueError: if the input model is invalid, a `ValueError` exception is raised.
    """
    epoch_ut1 = epoch.change_scale("UT1")
    epoch_tt = epoch.change_scale("TT")
    if model == 'IAU82':
        gmst = iauGMST82(epoch_ut1.jd, 0)
    elif model == 'IAU00':
        gmst = iauGMST00(epoch_ut1.jd, 0, epoch_tt.jd, 0)
    elif model == 'IAU06':
        gmst = iauGMST06(epoch_ut1.jd, 0, epoch_tt.jd, 0)
    else:
        raise ValueError(f"Invalid model {model}. Available models are IAU82, IAU00 or IAU06.")
    return gmst


def iauGMST82(dj1, dj2):
    """
    Compute Greenwich Mean Sidereal Time (IAU 1982 model).

    This is the canonical IAU 1982 expression for GMST, as implemented in
    the IAU SOFA software collection.

    Args:
        dj1 (float): First part of UT1 Julian Date.
        dj2 (float): Second part of UT1 Julian Date. The sum dj1+dj2 is the UT1 Julian Date.

    Returns:
        float: Greenwich Mean Sidereal Time (GMST) in radians, normalized to [0, 2π).

    Notes:
        - The UT1 Julian Date may be split between `dj1` and `dj2` in various ways
          (JD, MJD, J2000, date+time method, etc.), as long as their sum is correct.
        - The algorithm follows the IAU 1982 GMST-UT1 model.
        - Accuracy is preserved by using the fractional part of the day
          in seconds rather than whole days.
    """
    # Coefficients of IAU 1982 GMST-UT1 model
    A = 24110.54841 - SECONDS_IN_DAY / 2.0
    B = 8640184.812866
    C = 0.093104
    D = -6.2e-6

    # Split to preserve precision
    if dj1 < dj2:
        d1 = dj1
        d2 = dj2
    else:
        d1 = dj2
        d2 = dj1

    # Julian centuries since J2000.0
    t = (d1 + (d2 - Epoch.J2000)) / Epoch.DJC

    # Fractional part of UT1 day, in seconds
    f = SECONDS_IN_DAY * (math.fmod(d1, 1.0) + math.fmod(d2, 1.0))

    # GMST at this UT1
    gmst = anp(DS2R * ((A + (B + (C + D * t) * t) * t) + f))

    return gmst


def iauERA00(dj1, dj2):
    """
    Compute Earth Rotation Angle (IAU 2000 model).

    Args:
        dj1 (float): First part of UT1 Julian Date.
        dj2 (float): Second part of UT1 Julian Date. The sum dj1+dj2 is the UT1 Julian Date.

    Returns:
        float: Earth Rotation Angle (radians), normalized to [0, 2π).
    """
    # Order to preserve precision
    if dj1 < dj2:
        d1, d2 = dj1, dj2
    else:
        d1, d2 = dj2, dj1

    # Days since J2000.0
    t = d1 + (d2 - Epoch.J2000)

    # Fractional part of day
    f = math.fmod(d1, 1.0) + math.fmod(d2, 1.0)

    # Earth rotation angle
    theta = anp(2*PI * (f + 0.7790572732640 + 0.00273781191135448 * t))

    return theta


def iauGMST00(uta, utb, tta, ttb):
    """
    Compute Greenwich Mean Sidereal Time (IAU 2000 model).

    Args:
        uta (float): First part of UT1 Julian Date.
        utb (float): Second part of UT1 Julian Date. The sum uta+utb is UT1.
        tta (float): First part of TT Julian Date.
        ttb (float): Second part of TT Julian Date. The sum tta+ttb is TT.

    Returns:
        float: Greenwich Mean Sidereal Time (GMST) in radians, normalized to [0, 2π).

    Notes:
        - UT1 is required for Earth rotation (ERA).
        - TT is required for precession.
        - Algorithm is from Capitaine et al. (2003), consistent with IERS Conventions 2003.
    """
    # TT Julian centuries since J2000.0
    t = ((tta - Epoch.J2000) + ttb) / Epoch.DJC

    # GMST, IAU 2000 model
    gmst = anp(
        iauERA00(uta, utb) +
        (0.014506 +
         (4612.15739966 +
          (1.39667721 +
           (-0.00009344 +
            0.00001882 * t) * t) * t) * t) * DAS2R)

    return gmst


def iauGMST06(uta, utb, tta, ttb):
    """
    Compute Greenwich Mean Sidereal Time (IAU 2006 model).

    Args:
        uta (float): First part of UT1 Julian Date.
        utb (float): Second part of UT1 Julian Date. The sum uta+utb is UT1.
        tta (float): First part of TT Julian Date.
        ttb (float): Second part of TT Julian Date. The sum tta+ttb is TT.

    Returns:
        float: Greenwich Mean Sidereal Time (GMST) in radians, normalized to [0, 2π).

    Notes:
        - UT1 is required for Earth rotation (ERA).
        - TT is required for precession.
        - This GMST is consistent only with the IAU 2006 precession model.
    """
    # TT Julian centuries since J2000.0
    t = ((tta - Epoch.J2000) + ttb) / Epoch.DJC

    # Greenwich Mean Sidereal Time, IAU 2006
    gmst = anp(
        iauERA00(uta, utb) +
        (0.014506 +
         (4612.156534 +
          (1.3915817 +
           (-0.00000044 +
            (-0.000029956 +
             (-0.0000000368) * t) * t) * t) * t) * t) * DAS2R)

    return gmst
