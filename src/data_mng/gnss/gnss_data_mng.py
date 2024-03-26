from src.io.states import OUTPUT_FILENAME_MAP, get_file_header, export_to_file
from src.io.config import config_dict, EnumPositioningMode
from src.io.rinex_parser import RinexNavReader, RinexObsReader
from src.data_mng import Container
from src.common_log import IO_LOG, get_logger
from .navigation_data import NavigationData
from .observation_data import ObservationData


__all__ = ["GnssDataManager"]


class GnssDataManager(Container):
    __slots__ = ["nav_data", "obs_data", "isb_data", "nav_solution",
                 "smooth_obs_data", "iono_free_obs_data"]

    def __init__(self):
        super().__init__()

        self.nav_data = NavigationData()  # Input Navigation Data
        self.obs_data = ObservationData()  # Input Observation Data
        self.smooth_obs_data = ObservationData()  # Smooth Observation Data
        self.iono_free_obs_data = ObservationData()  # Iono Free Observation Data
        self.isb_data = None  # ISB (Inter-system bias) data
        self.nav_solution = None  # Navigation solution

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
            if data_name in self.__slots__:
                return getattr(self, data_name)
            else:
                raise ValueError(f'{data_name} is not available.')

    def get_clean_obs_data(self):
        # return the observation data for the PVT processing, depending on user configuration
        model = config_dict.get("model", "mode")
        if config_dict.get("preprocessor", "compute_smooth"):
            return self.smooth_obs_data
        elif model == EnumPositioningMode.SPS_IF:
            return self.iono_free_obs_data
        else:
            return self.obs_data

    def get_raw_obs_data(self):
        return self.obs_data

    def read_inputs(self, trace_dir):

        # TODO: add possibility of multiple files

        # TODO: check if precise products (clk and sp3) or navigation products (nav)

        log = get_logger(IO_LOG)

        # read navigation data
        nav_file = config_dict.get("inputs", "nav_files")[0]
        obs_file = config_dict.get("inputs", "obs_files")[0]
        services = config_dict.get_services()
        first_epoch = config_dict.get("inputs", "arc", "first_epoch")
        last_epoch = config_dict.get("inputs", "arc", "last_epoch")
        snr_check = config_dict.get("inputs", "snr_control")

        log.info("Launching RinexNavReader")
        RinexNavReader(nav_file, self.get_data("nav_data"))

        log.info("Launching RinexObsReader")
        RinexObsReader(self.get_data("obs_data"), obs_file, services, first_epoch, last_epoch, snr_check)

        self._trace_files(trace_dir)

    def _trace_files(self, trace_dir):
        # trace data files
        with open(f"{trace_dir}\\RawObservationData.txt", "w") as file:
            file.write(str(self.get_data("obs_data")))
        with open(f"{trace_dir}\\RawNavigationData.txt", "w") as file:
            file.write(str(self.get_data("nav_data")))

    def save_data(self, directory):
        log = get_logger(IO_LOG)
        log.info(f"storing data to {directory}...")
        file_list = {}

        for sim_data in ["nav_solution"]:
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
                            file_list[ext].write(f"{get_file_header(ext, state)}\n")
                            log.info(f"creating output file {filename}")

                        # save this epoch data
                        data = export_to_file(state, ext)
                        if isinstance(data, str):
                            file_list[ext].write(f"{time_str},{data}\n")
                        elif isinstance(data, list):
                            for entry in data:
                                file_list[ext].write(f"{time_str},{entry}\n")
