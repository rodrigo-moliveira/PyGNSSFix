"""Module with CSVData to store data from a CSV File"""
import pandas as pd
import os

# from src.io.states.csv_utils import read_csv_file, parse_satellite_data_file


class CSVData:
    """
    CSVData class.

    A class to represent and manage data from a CSV file.

    This class provides methods to load, access, and manipulate data from CSV files, which can include
    position, velocity, timeseries measurements, or other types of numerical data.

    """

    def __init__(self, name, description, time_cols, data_cols, units=None, legend=None, title=""):
        """
        Initializes the CSVData object.

        Arguments:
            name(str): Short name of the dataset.
            description(str): Description of the dataset.
            units(list[str] or None): optional. List of strings to specify the units of data. The dimension of the list
                must be consistent with the dimensions of the data in `self.data`. Default is None.
            legend(list[str] or None): optional. List of strings to specify the legend of the data for plotting purposes
                This is the name of each axis in the plot. The dimension of the list must be consistent with the
                dimensions of the data in `self.data`. Default is None.
            title(str): optional. Name of the title of the plot. Default is an empty string.
            time_cols(list[int] or None): optional. List with the column indexes in the csv file containing the time
                information
            data_cols(list[int] or None): optional. List with the column indexes in the csv file containing the data
        """
        self.name = name
        self.description = description
        self.title = title
        self.data = None
        self.data_cols = data_cols
        self.time_cols = time_cols
        self._n_epochs = 0
        self._dim = len(data_cols)  # dimension of the data for each array (1D, 2D, 3D, etc.)
        self.units = units
        self.legend = legend

        # if self.units is not None and len(units) != self._dim:
        #     raise ValueError(f"Inconsistency between length of arguments 'units' {units} and 'data_cols' "
        #                      f"{self.data_cols}. They should have the same dimension")
        # if self.legend is not None and self._dim is not None and len(legend) != self._dim:
        #     raise ValueError(f"Inconsistency between length of arguments 'legend' {legend} and 'data_cols' "
        #                      f"{self.data_cols}. They should have the same dimension")

    def read_data(self, file):
        """
        Read data from the provided file. Currently, no unit conversion is performed.

        Arguments:
            file(pathlib.Path or str): path to the input CSV file
        Raises:
            ValueError : an exception is raised if an inconsistency between the shape of the input CSV data and the
                configured dimension of the data is verified
        """
        if not os.path.exists(file):
            raise ValueError(f"File path {file} does not exist")
        print(f"reading file {file}...")  # TODO: add log message

        cols = list(self.time_cols) + list(self.data_cols)
        data = pd.read_csv(file, usecols=cols)

        n_epochs, n_dim = data.shape  # (n,m) shape, where n is the number of epochs and m is the row dimension

        # unpack shape
        n_dim -= len(self.time_cols)

        if self._dim != n_dim:
            raise ValueError(f"Inconsistency between the defined dimension of data (dim={self._dim}) and the shape "
                             f"of the input data {data.shape}. Please review the input CSV file")

        # Rename columns by index
        if self.legend is not None and len(self.data_cols) == len(self.legend):
            for i, new_name in zip(self.data_cols, self.legend):
                data.columns.values[i] = new_name

        self.data = data

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self.data.to_string())

    def is_empty(self):
        return self.data is None

    def to_data_array(self):
        # Select the clock_bias and sigma columns
        data_matrix = self.data.iloc[:, list(self.data_cols)].values
        return data_matrix

    def to_time_array(self):
        # Select the clock_bias and sigma columns
        time_matrix = self.data.iloc[:, list(self.time_cols)].values
        return time_matrix

    # def format_sat_data(self):
    #    cols = list(self.time_cols) + list(self.data_cols)
    #    self.data = parse_satellite_data_file(self.data, cols=cols)
