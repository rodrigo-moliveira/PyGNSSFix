"""Module with CSVData to store data from a CSV File"""
from src.io.states.csv_utils import read_csv_file, parse_satellite_data_file


class CSVData:
    """
    CSVData class.

    A class to represent and manage data from a CSV file.

    This class provides methods to load, access, and manipulate data from CSV files, which can include
    position, velocity, timeseries measurements, or other types of numerical data.

    """

    def __init__(self, name, description, units=None, legend=None, title="", cols=None, dtype=None):
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
            cols(list[int] or None): optional. List with the column indexes to read the associated csv file
            dtype(str or None): string specifying the data type to read the csv cells. (See `numpy.loadtxt`
                for more information)
        """
        self.name = name
        self.description = description
        self.title = title
        self.data = None
        self.cols = cols
        self._n_epochs = 0
        self._dim = None if cols is None else len(cols)  # dimension of the data for each array (1D, 2D, 3D, etc.)
        self.units = units
        self.legend = legend
        self.dtype = dtype

        if self.units is not None and self._dim is not None and len(units) != self._dim:
            raise ValueError(f"Inconsistency between length of arguments 'units' {units} and 'cols' {self.cols}. "
                             f"They should have the same dimension")
        if self.legend is not None and self._dim is not None and len(legend) != self._dim:
            raise ValueError(f"Inconsistency between length of arguments 'legend' {legend} and 'cols' {self.cols}. "
                             f"They should have the same dimension")

    def read_data(self, file):
        """
        Read data from the provided file. Currently, no unit conversion is performed.

        Arguments:
            file(pathlib.Path or str): path to the input CSV file
        Raises:
            ValueError : an exception is raised if an inconsistency between the shape of the input CSV data and the
                configured dimension of the data is verified
        """
        print(f"reading file {file}...")  # TODO: add log message
        data = read_csv_file(file, header=True, cols=self.cols, dtype=self.dtype)
        shape = data.shape  # (n,m) shape, where n is the number of epochs and m is the data dimension

        # unpack shape
        self._n_epochs = shape[0]
        n_dim = 1 if len(shape) == 1 else shape[1]

        if self._dim is None:
            self._dim = n_dim  # no data dimension has been previously specified -> simply update it
        else:
            if self._dim != n_dim:
                raise ValueError(f"Inconsistency between the defined dimension of data (dim={self._dim}) and the shape "
                                 f"of the input data {data.shape}. Please review the input CSV file")
        self.data = data

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self.data)

    def is_empty(self):
        return self.data is None

    def format_sat_data(self):
        # TODO: parei aqui
        self.data = parse_satellite_data_file(self.data, cols)
