from src.data_mng.containers.container import Container
from src.data_mng.csv_data import CSVData


class InsDataManager(Container):
    __slots__ = ["navigation_data", "observation_data",
                 "_available", "_do_not_save", "error_dict", "_plottable"]

    def __init__(self):
        super().__init__()

        # Reference time
        self.time = CSVData(name="time", description="sample time", units=["s"], legend=["time"])

        # Input Navigation Data
        self.navigation_data = CSVData(name="navigation_data", description="Input Navigation Data")
        # Input Observation Data
        self.observation_data = CSVData(name="observation_data", description="Input Observation Data")

        # available data for the current simulation
        self._available = []
        self._plottable = {"projection"}

        # SimulationData which is not intended to be saved
        self._do_not_save = []

        self.error_dict = {}

    def __str__(self):
        return f'{type(self).__name__}( DataManager for GNSS algorithms )'

    def __repr__(self):
        return str(self)

    def add_data(self, data_name, data, units=None):
        """
        Add data to available.
        Args:
            data_name: data name, str
            data: a scalar, a numpy array or a dict of the above two. If data is a dict, each
                value in it should be of same type (scalar or numpy array), same size and same
                units.
            units: Units of the data. If you know clearly no units conversion is needed, set
                units to None. If you do not know what units are used in the class InsDataMgr,
                you'd better provide the units of the data. Units conversion will be done
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

    def save_data(self, directory):
        for sim_data in self._available:
            if sim_data not in self._do_not_save:

                # fetch sim_data
                sim = getattr(self, sim_data, None)

                if sim is not None:
                    # print("saving", sim_data)
                    sim.save_to_file(directory, self.time)

    """def performance_evaluation(self):
        # pos, vel, att
        summary = ""

        for data, ref_data in zip(["pos", "vel", "att"], ["ref_pos", "ref_vel", "ref_att"]):
            _data = getattr(self, data, None)
            _ref_data = getattr(self, ref_data, None)
            if not _data.is_empty() and not _ref_data.is_empty():
                summary += self._evaluate(_data, _ref_data)

        return summary
    """

    @property
    def available(self):
        return self._available

