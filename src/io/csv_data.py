import csv
import numpy as np

from src.utils.units import convert_unit


def read_csv_file(filepath, **options):
    print(f"reading file {filepath}...")
    cols = options.get("columns", None)
    skip_first = 1 if options.get("header", True) else 0

    data = np.loadtxt(filepath,
                      skiprows=skip_first,
                      delimiter=',',
                      usecols=cols)
    return data


class CSVData:
    """
    Container to store data from a CSV file. The
    """

    def __init__(self, name, description, units=None, legend=None, title=""):
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
            legend (list): list of strings to specify legend of data
                To be stored in the first line of output csv files
                The length of units is the same as columns of each set of data in `self.data`.
        """
        self.name = name
        self.description = description
        self.title = title
        # units of self.data
        if units is None:
            self.units = ['']
        else:
            self.units = list(units)

        if len(legend) != len(units):
            raise ValueError(f"legend {legend} and units {units} should have the same dimension")
        self.legend = legend

        self.data = None

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

    def is_empty(self):
        return self.data is None
