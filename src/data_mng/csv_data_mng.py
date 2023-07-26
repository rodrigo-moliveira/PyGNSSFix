from src.constants import OUTPUT_FILENAME_MAP
from src.data_mng.container import Container
from src.data_mng.csv_data import CSVData
from src.io.csv_reader import CSVReader


# TODO: quando chegar ao INS, criar um class mae com as features generices e depois separar em duas

class GnssRunStorageManager(Container):
    __slots__ = ["position", "velocity", "prefit_residuals", "clock_bias",
                 "postfit_residuals", "dop", "satellite_azel", "_available"]

    def __init__(self):
        super().__init__()

        # position in ECEF
        self.position = CSVData(name="position",
                                description="position in the ECEF frame",
                                units=['m', 'm', 'm'],
                                legend=['pos_x', 'pos_y', 'pos_z'],
                                title="Position (ECEF)")

        # clock bias
        self.clock_bias = CSVData(name="clock_bias",
                                  description="Receiver Clock Bias",
                                  units=['s'],
                                  legend=['time'],
                                  title="Receiver Clock Bias")

        # velocity in ECEF
        self.velocity = CSVData(name="velocity",
                                description="velocity in the ECEF frame",
                                units=['m/s', 'm/s', 'm/s'],
                                legend=['vel_x', 'vel_y', 'vel_z'],
                                title="Velocity (ECEF)")

        # Dilution of Precision (DOP)
        self.dop = CSVData(name="DOP",
                           description="Dilution of Precision in ECEF frame",
                           units=['', '', '', ''],
                           legend=[' ', ' ', ' ', ' '],
                           title="Dilution of Precision (ECEF)")

        # prefit residuals
        self.prefit_residuals = CSVData(name="prefit_residuals",
                                        description="Prefit Residuals",
                                        units=['m'],
                                        legend=[' '],
                                        title="Prefit Residuals")

        # postfit residuals
        self.postfit_residuals = CSVData(name="postfit_residuals",
                                         description="Postfit Residuals",
                                         units=['m'],
                                         legend=[' '],
                                         title="Postfit Residuals")

        # available data for the current simulation
        self._available = []

    def __str__(self):
        return f'{type(self).__name__}( DataManager for GNSS Run Performance Evaluation )'

    def add_data(self, data_name, data, units=None):
        """
        Add data to available.
        Args:
            data_name: data name, str
            data: a scalar, a numpy array or a dict of the above two. If data is a dict, each
                value in it should be of same type (scalar or numpy array), same size and same
                units.
            units: Units of the data. If you know clearly no units convertion is needed, set
                units to None. If you do not know what units are used in the class InsDataMgr,
                you'd better provide the units of the data. Units convertion will be done
                automatically.
                If data is a scalar, units should be a list of one string to define its unit.
                If data is a numpy of size(m,n), units should be a list of n strings
                to define the units.
        """
        if data_name in self.__slots__:
            sim = getattr(self, data_name, None)
            if sim is not None:
                sim.add_data(data, units)

                # add to 'available' list
                if data_name not in self._available:
                    self._available.append(data_name)
        else:
            raise ValueError(f"Unsupported data: {data_name}, not in {self.__slots__}")

    def get_data(self, data_names):
        """
        Get data section of data_names.
        Args:
            data_names: a list of data names
        Returns:
            data: a list of data corresponding to data_names.
            If there is any unavailable data in data_names, return None
        """
        # single data
        if isinstance(data_names, str):
            if data_names in self._available:
                return getattr(self, data_names).data
            else:
                raise ValueError(f'{data_names} is not available.')

        # vector data
        data = []
        for i in data_names:
            if i in self._available:
                data.append(getattr(self, i).data)
            else:
                raise ValueError(f'{i} is not available.')
        return data

    """   
        def save_data(self, directory):
        for sim_data in self._available:
            if sim_data not in self._do_not_save:

                # fetch sim_data
                sim = getattr(self, sim_data, None)

                if sim is not None:
                    # print("saving", sim_data)
                    sim.save_to_file(directory, self.time)
    """

    @property
    def available(self):
        return self._available

    def read_data(self, output_folder):
        position = CSVReader.read_csv_file(output_folder/OUTPUT_FILENAME_MAP["position"],
                                           ignore_header=True)
        dop = CSVReader.read_csv_file(output_folder / OUTPUT_FILENAME_MAP["dop"],
                                      ignore_header=True)
        clock_bias = CSVReader.read_csv_file(output_folder / OUTPUT_FILENAME_MAP["clock_bias"],
                                             ignore_header=True)
        satellite_azel = CSVReader.read_csv_file(output_folder / OUTPUT_FILENAME_MAP["satellite_azel"],
                                                 ignore_header=True)

        self.add_data("position", position)
        self.add_data("dop", dop)
        self.add_data("clock_bias", clock_bias)
        self.add_data("satellite_azel", satellite_azel)
