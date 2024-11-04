""" Data Manager for GNSS Algorithms
"""
import os

from src.io.states import OUTPUT_FILENAME_MAP, get_file_header, export_to_file
from src.io.config import config_dict, EnumAlgorithmPNT
from src.io.rinex_parser import RinexNavReader, RinexObsReader
from src.data_mng import Container
from src.common_log import IO_LOG, get_logger
from .navigation_data import NavigationData
from .observation_data import ObservationData
from .sat_clock_data import SatelliteClocks
from .sat_orbit_data import SatelliteOrbits

__all__ = ["GnssDataManager"]


class GnssDataManager(Container):
    """
    Data Manager class that stores all necessary data for a GNSS run. Derives from the :py:class:`Container`
    class.

    Attributes:
        nav_data(NavigationData): navigation data object containing RINEX NAV ephemerides
        obs_data(ObservationData): raw observation data from RINEX OBS
        sat_clocks(SatelliteClocks): manager of satellite clocks (precise or navigation clocks)
        sat_orbits(SatelliteOrbits): manager of satellite orbits (precise or navigation orbits)
        smooth_obs_data(ObservationData): processed smooth observation data
        iono_free_obs_data(ObservationData): processed iono-free observation data
        nav_solution(list): navigation solution, list of :py:class:`src.data_mng.gnss.state_space.GnssStateSpace`
            objects

    """
    __slots__ = [
        "nav_data",             # Input
        "obs_data",             # Input
        "sat_clocks",           # Input
        "sat_orbits",           # Input
        "smooth_obs_data",      # Internal data
        "iono_free_obs_data",   # Internal data
        "nav_solution"          # Output
    ]

    def __init__(self):
        """ Default constructor with no arguments """
        super().__init__()

        self.nav_data = NavigationData()  # Input Navigation Data
        self.obs_data = ObservationData()  # Input Observation Data
        self.smooth_obs_data = ObservationData()  # Smooth Observation Data
        self.sat_clocks = SatelliteClocks()  # Satellite clocks manager (precise or navigation clocks)
        self.sat_orbits = SatelliteOrbits()  # Satellite orbits manager (precise or navigation orbits)
        self.iono_free_obs_data = ObservationData()  # Iono Free Observation Data
        self.nav_solution = None  # Navigation solution

    def __str__(self):
        return f'{type(self).__name__}( DataManager for GNSS algorithms )'

    def __repr__(self):
        return str(self)

    def add_data(self, data_name, data):
        """
        Adds the provided data to this DataManager instance.

        Args:
            data_name(str): the string name of the data to fill, matching one of the class attributes
            data(Any): data object to be stored

        Raises:
            AttributeError: an `AttributeError` is raised if the `data_name` string does not match one of the
                class attributes
        """
        if data_name in self.__slots__:
            setattr(self, data_name, data)
        else:
            raise AttributeError(f"Unsupported data: {data_name}, not in {self.__slots__}")

    def get_data(self, data_name):
        """
        Args:
            data_name(str): the data to be fetched
        Returns:
            Any: returns the queried data for the provided `data_name` attribute
        """
        # single data
        if isinstance(data_name, str):
            if data_name in self.__slots__:
                return getattr(self, data_name)
            else:
                raise AttributeError(f'{data_name} is not available.')

    def get_clean_obs_data(self):
        """
        Returns:
            ObservationData: returns the observation data for the PVT processing, depending on user configuration
        """
        model = config_dict.get("model", "mode")
        if config_dict.get("preprocessor", "compute_smooth"):
            return self.smooth_obs_data
        elif model == EnumAlgorithmPNT.SPS_IF:
            return self.iono_free_obs_data
        else:
            return self.obs_data

    def get_raw_obs_data(self):
        """
        Returns:
            ObservationData: returns the raw observation data for the PVT processing
        """
        return self.obs_data

    def read_inputs(self, trace_dir):
        """
        Function to read the input data. The data files are read from the user configurations
        according to the chosen algorithm.
        There are mandatory inputs for each model (SPS, PPP).

        Args:
            trace_dir(str): path to write trace files
        Raises:
            IOError: an exception is raised if the trace directory is not created successfully
        """
        log = get_logger(IO_LOG)

        # read navigation data
        nav_files = config_dict.get("inputs", "nav_files")
        obs_files = config_dict.get("inputs", "obs_files")
        clock_files = config_dict.get("inputs", "clk_files")
        sp3_files = config_dict.get("inputs", "sp3_files")
        use_precise_products = config_dict.get("inputs", "use_precise_products")
        gal_nav_type = config_dict.get("model", "GAL", "nav_type")

        log.info("Launching RinexObsReader")
        for file in obs_files:
            RinexObsReader(file, self.get_data("obs_data"))

        log.info(f"Using precise orbit/clock products: {use_precise_products}")

        log.info(f'Galileo messages selected by user are {gal_nav_type}')
        log.info('Launching RinexNavReader.')
        for file in nav_files:
            RinexNavReader(file, self.get_data("nav_data"), gal_nav_type)

        log.info("Launching SatelliteClocks constructor")
        self.sat_clocks.init(self.get_data("nav_data"), clock_files, use_precise_products)

        log.info("Launching SatelliteOrbits constructor")
        self.sat_orbits.init(self.get_data("nav_data"), sp3_files, use_precise_products)

        if config_dict.get("inputs", "trace_files"):
            self._trace_files(trace_dir, use_precise_products)

    def _trace_files(self, trace_dir, use_precise_products):

        inputs_dir = f"{trace_dir}\\inputs"
        try:
            os.makedirs(inputs_dir)
        except:
            raise IOError(f"Cannot create dir: {inputs_dir}")
        # trace data files
        with open(f"{inputs_dir}\\RawObservationData.txt", "w") as file:
            file.write(str(self.get_data("obs_data")))
        with open(f"{inputs_dir}\\RawNavigationData.txt", "w") as file:
            file.write(str(self.get_data("nav_data")))
        if use_precise_products:
            with open(f"{inputs_dir}\\PreciseClocks.txt", "w") as file:
                file.write(str(self.get_data("sat_clocks")))
            with open(f"{inputs_dir}\\PreciseOrbits.txt", "w") as file:
                file.write(str(self.get_data("sat_orbits")))

    def save_data(self, directory):
        """ Saves the navigation solution (contained in the `nav_solution` attribute) to the output file """
        log = get_logger(IO_LOG)
        log.info(f"storing data to {directory}...")
        file_list = {}

        # Save state variables (navigation solution)
        sim = getattr(self, "nav_solution", None)
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
                    data = export_to_file(ext, state)
                    if isinstance(data, str):
                        file_list[ext].write(f"{time_str},{data}\n")
                    elif isinstance(data, list):
                        for entry in data:
                            file_list[ext].write(f"{time_str},{entry}\n")

        # save observation data
        obs_data = self.get_clean_obs_data()
        ext = 'obs'
        filename = f"{directory}\\{OUTPUT_FILENAME_MAP[ext]}"
        file_list[ext] = open(filename, "w")
        file_list[ext].write(f"{obs_data.to_csv_file()}")
        log.info(f"creating output file {filename}")

        # close all files
        for ext in file_list.keys():
            file_list[ext].close()
