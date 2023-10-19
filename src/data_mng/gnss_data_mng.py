from src.constants import OUTPUT_FILENAME_MAP
from src.data_mng.container import Container
from src.io.states.export_states import get_file_header, export_to_file


class GnssDataManager(Container):
    __slots__ = ["nav_data", "obs_data", "nav_solution", "_to_save",
                 "_available", "_do_not_save"]

    def __init__(self):
        super().__init__()

        # Input Navigation Data
        self.nav_data = None
        # Input Observation Data
        self.obs_data = None

        # available data for the current simulation
        self._available = []
        self._to_save = ["nav_solution"]

    def __str__(self):
        return f'{type(self).__name__}( DataManager for GNSS algorithms )'

    def __repr__(self):
        return str(self)

    def add_data(self, data_name, data):
        """
        Add data to available.
        Args:
            data_name: data name, str
            data: a scalar, a numpy array or a dict of the above two. If data is a dict, each
                value in it should be of same type (scalar or numpy array), same size and same
                units.
        """
        if data_name in self.__slots__:
            setattr(self, data_name, data)
            # add to 'available' list
            if data_name not in self._available:
                self._available.append(data_name)
        else:
            raise ValueError(f"Unsupported data: {data_name}, not in {self.__slots__[0:-2]}")

    def get_data(self, data_name):
        """
        Get data section of data_names.
        Args:
            data_name: a list of data names
        Returns:
            data: a list of data corresponding to data_names.
            If there is any unavailable data in data_names, return None
        """
        # single data
        if isinstance(data_name, str):
            if data_name in self._available:
                return getattr(self, data_name)
            else:
                raise ValueError(f'{data_name} is not available.')

    def save_data(self, directory, log):
        log.info(f"storing data to {directory}...")
        file_list = {}

        for sim_data in self._available:
            if sim_data in self._to_save:
                # fetch sim_data
                sim = getattr(self, sim_data, None)
                if sim is not None:

                    # iterate over estimated states
                    for state in sim:

                        week, sow = state.epoch.gnss_time
                        time_str = f"{week},{sow}"

                        # save estimated data for this epoch
                        exportable_lst = state.get_exportable_lst()
                        for ext in exportable_lst:

                            # add this estimable to the file list (only do this once)
                            if ext not in file_list:
                                filename = f"{directory}\\{OUTPUT_FILENAME_MAP[ext]}"
                                file_list[ext] = open(filename, "w")
                                file_list[ext].write(f"{get_file_header(ext)}\n")
                                log.info(f"creating output file {filename}")

                            # save this epoch data
                            data = export_to_file(state, ext)
                            if isinstance(data, str):
                                file_list[ext].write(f"{time_str},{data}\n")
                            elif isinstance(data, list):
                                for entry in data:
                                    file_list[ext].write(f"{time_str},{entry}\n")

    @property
    def available(self):
        return self._available
