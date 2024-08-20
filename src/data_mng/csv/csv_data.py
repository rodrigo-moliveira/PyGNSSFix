"""Module with CSVData to store data from a CSV File"""
import pandas as pd
import os


class CSVData:
    """
    CSVData class.

    A class to represent, store and manage data from a CSV file.

    This class provides methods to load, access, and manipulate data from CSV files, which can include
    position, velocity, timeseries measurements, or other types of numerical data.

    """

    def __init__(self, name, description, time_cols, data_cols, units=None, legend=None, title=""):
        """
        Initializes the CSVData object.

        Arguments:
            name(str): Short name of the dataset.
            description(str): Description of the dataset.
            time_cols(list[int] or None): List with the column indexes in the csv file containing the time
                information
            data_cols(list[int] or None): List with the column indexes in the csv file containing the data
            units(list[str] or None): optional. List of strings to specify the units of data. The dimension of the list
                must be consistent with the dimension of the `data_cols` list. Default is None.
            legend(list[str] or None): optional. List of strings to specify the legend of the data for plotting purposes
                This is the name of each axis in the plot. The dimension of the list must be consistent with the
                dimension of the `data_cols` list. Default is None.
            title(str): optional. Name of the title of the plot. Default is an empty string.

        """
        self.name = name
        self.description = description
        self.title = title
        self.data = None
        self.data_cols = data_cols
        self.time_cols = time_cols
        self._dim = len(data_cols)  # dimension of the data array
        self.units = units
        self.legend = legend

        if self.units and len(self.units) != self._dim:
            raise ValueError(f"Inconsistency between length of arguments 'units' {units} and 'data_cols' "
                             f"{self.data_cols}. They should have the same dimension")
        if self.legend and len(self.legend) != self._dim:
            raise ValueError(f"Inconsistency between length of arguments 'legend' {legend} and 'data_cols' "
                             f"{self.data_cols}. They should have the same dimension")

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

        cols = list(self.time_cols) + list(self.data_cols)
        data = pd.read_csv(file, usecols=cols)

        n_epochs, n_dim = data.shape  # (n,m) shape, where n is the number of epochs and m is the row dimension

        # unpack shape
        n_dim -= len(self.time_cols)

        if self._dim != n_dim:
            raise ValueError(f"Inconsistency between the defined dimension of data (dim={self._dim}) and the shape "
                             f"of the input data {data.shape}. Please review the input CSV file")

        # Rename columns by index
        if self.legend:
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
        """ Return a matrix with the data columns """
        data_matrix = self.data.iloc[:, list(self.data_cols)].values
        return data_matrix

    def to_time_array(self):
        """ Return a matrix with the time columns """
        time_matrix = self.data.iloc[:, list(self.time_cols)].values
        return time_matrix
