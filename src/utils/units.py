import numpy as np
import src.constants as constant
from src.errors import UnknownConversionError

# TODO: add this functionality to TimeSeries
# TODO: consider using a units library instead
def convert_unit(data, src_unit, out_unit):
    """
    Unit conversion. Notice not to change values in data
    Args:
        data: convert data units from src_unit to out_unit. Data should be a scalar,
            a numpy array of size(n,1) or (n,m). n is data length, m is data dimension.
        src_unit: a list of unit of the input data.
        out_unit: a list of unit we want to convert the data to (output).
    Returns:
        x: data after unit conversion (from src_unit to out_unit).
    """

    if isinstance(src_unit, str):
        src_unit = [src_unit]
    if isinstance(out_unit, str):
        out_unit = [out_unit]

    scale_to_apply = unit_conversion_scale(src_unit, out_unit)

    # unit conversion
    x = data.copy()  # avoid changing values in data

    x = convert_unit_ndarray_scalar(x, scale_to_apply)
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
            scale[i] = constant.DEG2RAD
        elif src_unit[i] == 'deg/s' and out_unit[i] == 'rad/s':
            scale[i] = constant.DEG2RAD
        elif src_unit[i] == 'deg/hr' and out_unit[i] == 'rad/s':
            scale[i] = constant.DEG2RAD / 3600.0
        elif src_unit[i] == 'deg/sqrt(hr)' and out_unit[i] == 'rad/sqrt(s)':
            scale[i] = constant.DEG2RAD / 60.0

        # rad to deg
        elif src_unit[i] == 'rad' and out_unit[i] == 'deg':
            scale[i] = 1.0 / constant.DEG2RAD
        elif src_unit[i] == 'rad/s' and out_unit[i] == 'deg/s':
            scale[i] = 1.0 / constant.DEG2RAD
        elif src_unit[i] == 'rad/s' and out_unit[i] == 'deg/hr':
            scale[i] = 3600.0 / constant.DEG2RAD

        # g to m/s^2
        elif src_unit[i] == 'g' and out_unit[i] == 'm/s^2':
            scale[i] = constant.EARTH_G0
        elif src_unit[i] == 'mg' and out_unit[i] == 'm/s^2':
            scale[i] = 1E-3 * constant.EARTH_G0
        elif src_unit[i] == 'm/s/sqrt(hr)' and out_unit[i] == 'm/s/sqrt(s)':
            scale[i] = 1 / 60

        elif src_unit[i] == 'ppm' and out_unit[i] == '':
            scale[i] = 1E-6

        # m to km
        elif src_unit[i] == 'm' and out_unit[i] == 'km':
            scale[i] = 1 / 1000

        # km to m
        elif src_unit[i] == 'km' and out_unit[i] == 'm':
            scale[i] = 1000

        # Unknown conversion!
        else:
            if src_unit[i] != out_unit[i]:
                raise UnknownConversionError(f"Cannot convert unit from {src_unit[i]} in {src_unit} to {out_unit[i]}.")
    return scale


def convert_unit_ndarray_scalar(x, scale):
    """
    Unit conversion of numpy array or a scalar.
    Args:
        x: convert x units from src_unit to out_unit. x should be a scalar,
            a numpy array of size(m,) or (n,m). n is x length (in time), m is x dimension.
        scale: 1D numpy array of unit conversion scale. x = x * scale
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
        raise ValueError('Input x should be a scalar, 1D or 2D array, ndim = %s' % x.ndim)
    return x
