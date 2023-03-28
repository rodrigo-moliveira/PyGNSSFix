import numpy as np

from src.utils.units import convert_unit


class CSVData:
    """
    Container to store data from a CSV file. The
    """
    def __init__(self, name, description, units=None, output_units=None, legend=None, ignore_first_row=False,
                 title=""):
        """
        Set up data properties (input/output data). The provided data can be a TimeSeries, numpy array, scalar
        * In case of numpy array, it is of shape (m,n)
            m is the time dimension
            n is the length of the vector for each epoch (may be scalar, tri-dimensional, etc.)
        * In case of TimeSeries see :class:`src.data_mng.containers.timeseries.TimeSeries`
        * In case of scalar, it is just a constant scalar (no timeseries associated)

        Args:
            name (str): string name of the data
            description (str): string description of the data
            units (list): a list of strings to specify units of data that is going to be used in internal
                computations. The length of units must be consisted with the provided data in `self.data`
            output_units (list): a list of strings to specify units of data when being plotted or saved to
                file. If this is set to None, output_units will be the same as units, and no unit conversion is needed.
            legend (list): list of strings to specify legend of data
                To be stored in the first line of output csv files
                The length of units is the same as columns of each set of data in `self.data`.
            ignore_first_row(bool): whether to ignore the first row (first epoch) of the data matrix.
                This is used to ignore the data used for initialization purposes
                (for example, when computing gyro readouts)
        """
        self.name = name
        self.description = description
        self.title = title
        # units of self.data
        if units is None:
            self.units = ['']
        else:
            self.units = list(units)
        # output units should have same length as units
        if output_units is None:
            self.output_units = self.units
        else:
            self.output_units = list(output_units)
            len_in = len(self.units)
            len_out = len(self.output_units)
            if len_in != len_out:
                raise ValueError(f"units {units} and output units {output_units} should have the same dimension")

        if len(legend) != len(units):
            raise ValueError(f"legend {legend} and units {units} should have the same dimension")
        self.legend = legend

        self.data = None
        self.ignore_first_row = ignore_first_row

    def add_data(self, data, units=None):
        """
        Add data to SimulatedData.
        Args:
            data: a scalar, a numpy array or TimeSeries of the above two.
            units: Units of the input data. If you know clearly no units conversion is needed, set
                units to None. If you do not know what units are used in the algorithm,
                you'd better provide the units of the data. Units conversion will be done
                automatically here.
                If data is a scalar, units should be a list of one string to define its unit.
                If data is a numpy of size(m,n), units should be a list of n strings
                to define the units.
        """
        if units is not None:
            units = list(units)
            if len(units) == len(self.units):
                if units != self.units:
                    # data units are different from units in the manager, need to be converted
                    print(f"converting {self.name} from {units} to {self.units}")
                    data = convert_unit(data, units, self.units)
            else:
                raise ValueError(f'Units {units} and {self.units} are of different lengths. Error')

        self.data = data

    def __str__(self):
        return self.name

    def save_to_file(self, directory, time_obj=None):
        print(f"saving {self.name} to file")
        # TODO: log message...
        f = open(directory+f"\\{self.name}.csv", "w")

        time_data = None
        time_header = None

        if time_obj is not None:
            time_data = time_obj.data
            time_header = f"{time_obj.legend[0]}[{time_obj.units[0]}]"

        # Write Header
        f.write(f"{self._header_to_file(time_str=time_header)}\n")

        # Write Data
        _data = self._data_to_file(time_data)
        for line in _data:
            f.write(f"{line}\n")
        f.close()

    def _header_to_file(self, time_str=None):
        # check dimension _dim of vector
        try:
            _len, _dim = self.data.shape
        except:
            # scalar dimension
            _dim = 1
            _len = len(self.data)

        header = ""

        if time_str is not None:
            header += str(time_str)+","

        for i in range(_dim):
            header += f"{self.legend[i]}[{self.output_units[i]}]" + ","

        # remove trailing ','
        header = header[0:-1]

        return header

    def _data_to_file(self, time_arr=None):
        # NOTE: I'm currently assuming that it's always np.array

        data_converted = convert_unit(self.data, self.units, self.output_units)

        # update self.data and self.units
        # self.data = data_converted
        # self.units = self.output_units

        # check dimension _dim of vector
        try:
            _len, _dim = data_converted.shape
        except:
            # scalar dimension
            _dim = 1
            _len = len(data_converted)

        _data = []

        # whether to ignore the first row (epoch) of data
        _start = 0 if self.ignore_first_row is False else 1

        for i in range(_start, _len):
            line = ""

            # check if there are NAN entries in this line (ignore line if True)
            if np.isnan(np.sum(data_converted[i])):
                continue

            if time_arr is not None:
                line += f"{time_arr[i,0]}" + ","

            for j in range(_dim):
                line += str(data_converted[i, j]) + ","

            # remove trailing ','
            line = line[0:-1]
            _data.append(line)

        return _data

    def is_empty(self):
        return self.data is None
