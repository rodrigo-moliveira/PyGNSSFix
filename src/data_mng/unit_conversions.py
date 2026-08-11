import numpy as np
import pandas as pd

from src.constants import DEG2RAD, EARTH_G0


def convert_unit(data, src_unit, out_unit):
    """
    Unit conversion. Notice not to change values in data
    Args:
        data: convert data units from src_unit to dst_unit. Data should be a scalar,
            a numpy array of size(n,) or (n,m). n is data length, m is data dimension.
        src_unit: a list of unit of the input data.
        out_unit: a list of unit we want to convert the data to (output).
    Returns:
        x: data after unit conversion (from src_unit to out_unit).
    """

    # if isinstance(src_unit, str):
    #     src_unit = [src_unit]
    # if isinstance(out_unit, str):
    #     out_unit = [out_unit]

    scale = unit_conversion_scale(src_unit, out_unit)

    if isinstance(data, pd.DataFrame):
        x = data.to_numpy().copy()
    else:
        x = data.copy()  # avoid changing values in data

    x = convert_unit_ndarray_scalar(x, scale)
    return x


def unit_conversion_scale(src_unit, out_unit):
    """
    Calculate unit conversion scale.
    """
    len_in = len(src_unit)
    len_out = len(out_unit)
    if len_out != len_in:
        raise ValueError(f"Units {src_unit} and {out_unit} should have consistent dimensions")

    scale = np.ones((len_out,))
    for i in range(len_out):

        # deg to rad
        if src_unit[i] == 'deg' and out_unit[i] == 'rad':
            scale[i] = DEG2RAD
        elif src_unit[i] == 'deg/s' and out_unit[i] == 'rad/s':
            scale[i] = DEG2RAD
        elif src_unit[i] == 'deg/hr' and out_unit[i] == 'rad/s':
            scale[i] = DEG2RAD / 3600.0
        elif src_unit[i] == 'deg/sqrt(hr)' and out_unit[i] == 'rad/sqrt(s)':
            scale[i] = DEG2RAD / 60.0

        # rad to deg
        elif src_unit[i] == 'rad' and out_unit[i] == 'deg':
            scale[i] = 1.0 / DEG2RAD
        elif src_unit[i] == 'rad/s' and out_unit[i] == 'deg/s':
            scale[i] = 1.0 / DEG2RAD
        elif src_unit[i] == 'rad/s' and out_unit[i] == 'deg/hr':
            scale[i] = 3600.0 / DEG2RAD

        # g to m/s^2
        elif src_unit[i] == 'g' and out_unit[i] == 'm/s^2':
            scale[i] = EARTH_G0
        elif src_unit[i] == 'mg' and out_unit[i] == 'm/s^2':
            scale[i] = 1E-3 * EARTH_G0
        elif src_unit[i] == 'm/s/sqrt(hr)' and out_unit[i] == 'm/s/sqrt(s)':
            scale[i] = 1 / 60

        elif src_unit[i] == 'ppm' and out_unit[i] == '':
            scale[i] = 1E-6

        # Unknown convertion!
        else:
            if src_unit[i] != out_unit[i]:
                from src.errors import UnknownConversion
                raise UnknownConversion(f"Cannot convert unit from {src_unit[i]} in {src_unit} to {out_unit[i]}.")
    return scale


def convert_unit_ndarray_scalar(x, scale):
    """
    Unit conversion of numpy array or a scalar.
    Args:
        x: convert x units from src_unit to out_unit. x should be a scalar,
            a numpy array of size(m,) or (n,m). n is x length (in time), m is x dimension.
        scale: 1D numpy array of unit convertion scale. x = x * scale
    Returns:
        x: x after unit conversion.
    """
    scale_m = scale.shape[0]

    if isinstance(x, np.ndarray) and (x.ndim == 2 or x.ndim == 1):
        if x.ndim == 2:
            # 2D array (time, vector)
            x_m = x.shape[1]

            for i in range(x_m):
                if scale_m == 1:
                    if scale[0] != 1.0:
                        x[:, i] = x[:, i] * scale[0]
                else:
                    if scale[i] != 1.0:
                        x[:, i] = x[:, i] * scale[i]

        elif x.ndim == 1:
            # 1D array (vector)
            if len(x) == scale_m:
                x = x * scale
            else:
                x = x * scale[0]
    elif isinstance(x, (int, float)):
        x = x * scale[0]
    else:
        raise ValueError('Input x should be a scalar, 1D or 2D array, ndim = %s'% x.ndim)
    return x