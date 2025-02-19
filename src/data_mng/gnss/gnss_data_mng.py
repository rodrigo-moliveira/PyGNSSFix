""" Data Manager for GNSS Algorithms
"""
import os

from src.io.states import OUTPUT_FILENAME_MAP, get_file_header, export_to_file
from src.io.config import config_dict, EnumObservationModel, EnumAlgorithmPNT, EnumSatelliteBias
from src.io.rinex_parser import RinexNavReader, RinexObsReader, AntexReader
from src.data_mng import Container
from src.common_log import IO_LOG, get_logger
from .navigation_data import NavigationData
from .observation_data import ObservationData
from .phase_center_mng import PhaseCenterManager
from .sat_clock_data import SatelliteClocks
from .sat_orbit_data import SatelliteOrbits
from .bias_manager import BiasManager
from .global_iono_map import GlobalIonoMap

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
        iono_gim(GlobalIonoMap): manager of global ionospheric maps (VTEC)
        phase_center(PhaseCenterManager): manager of phase center data
        sat_bias(BiasManager): manager of satellite code and phase biases (precise or navigation bias)
        smooth_obs_data(ObservationData): processed smooth observation data
        iono_free_obs_data(ObservationData): processed iono-free observation data
        nav_solution(list): navigation solution, list of :py:class:`src.data_mng.gnss.state_space.GnssStateSpace`
            objects

    """
    __slots__ = [
        "nav_data",            # Input
        "obs_data",            # Input
        "sat_clocks",          # Input
        "sat_orbits",          # Input
        "iono_gim",            # Input
        "phase_center",        # Input
        "sat_bias",            # Input
        "smooth_obs_data",     # Internal data
        "iono_free_obs_data",  # Internal data
        "nav_solution"         # Output
    ]

    def __init__(self):
        """ Default constructor with no arguments """
        super().__init__()

        self.nav_data = NavigationData()  # Input Navigation Data
        self.obs_data = ObservationData()  # Input Observation Data
        self.smooth_obs_data = ObservationData()  # Smooth Observation Data
        self.sat_clocks = SatelliteClocks()  # Satellite clocks manager (precise or navigation clocks)
        self.iono_gim = GlobalIonoMap()  # Global Ionospheric (VTEC) Map Manager
        self.phase_center = PhaseCenterManager()  # Phase Center Data Manager
        self.sat_bias = BiasManager()  # Satellite Code and Phase Bias Manager
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
        obs_model = config_dict.get('obs_model')
        if config_dict.get("preprocessor", "compute_smooth"):
            return self.smooth_obs_data
        elif obs_model == EnumObservationModel.COMBINED:
            return self.iono_free_obs_data
        else:
            return self.obs_data

    def get_raw_obs_data(self):
        """
        Returns:
            ObservationData: returns the raw observation data for the PVT processing
        """
        return self.obs_data

    def read_inputs(self, gnss_alg: EnumAlgorithmPNT, trace_dir):
        """
        Function to read the input data. The data files are read from the user configurations
        according to the chosen algorithm.
        There are mandatory inputs for each GNSS Algorithm (SPS, PR-PPP).

        Args:
            gnss_alg(EnumAlgorithmPNT): GNSS Algorithm Enumeration
            trace_dir(str): path to write trace files
        Raises:
            IOError: an exception is raised if the trace directory is not created successfully
        """
        log = get_logger(IO_LOG)

        # read observation data
        obs_files = config_dict.get("inputs", "obs_files")
        GnssDataManager.check_input_list("inputs.obs_files", obs_files, log)
        log.info("Launching RinexObsReader.")

        for file in obs_files:
            RinexObsReader(file, self.get_data("obs_data"), self.get_data("phase_center"))

        # Load specific inputs for each PNT Algorithm
        if gnss_alg == EnumAlgorithmPNT.SPS:
            log.info("In SPS Mode, GNSS orbits and clocks are provided from broadcast ephemerides (RINEX NAV).")

            nav_files = config_dict.get("inputs", "nav_files")
            gal_nav_type = config_dict.get("model", "GAL", "nav_type")
            GnssDataManager.check_input_list("inputs.nav_files", nav_files, log)
            log.info(f'Galileo messages selected by user are {gal_nav_type}.')
            log.info('Launching RinexNavReader.')

            for file in nav_files:
                RinexNavReader(file, self.get_data("nav_data"), gal_nav_type)

            log.info("Launching SatelliteClocks constructor (clocks from broadcast ephemerides).")
            self.sat_clocks.init(self.get_data("nav_data"), None, False)

            log.info("Launching SatelliteOrbits constructor (orbits from broadcast ephemerides).")
            self.sat_orbits.init(self.get_data("nav_data"), None, False)

            log.info("Launching Satellite Code Bias Manager with BGD/TGD data from broadcast ephemerides.")
            self.sat_bias.init(self.get_data("nav_data"), None, EnumSatelliteBias.BROADCAST)

        elif gnss_alg == EnumAlgorithmPNT.PR_PPP:
            log.info("In PR-PPP Mode, GNSS orbits and clocks are provided from precise products (SP3 and CLK files).")

            clock_files = config_dict.get("inputs", "clk_files")
            sp3_files = config_dict.get("inputs", "sp3_files")
            dcb_files = config_dict.get("inputs", "dcb_files")
            osb_files = config_dict.get("inputs", "osb_files")
            ionex_files = config_dict.get("inputs", "ionex_files")
            antex_files = config_dict.get("inputs", "antex_files")
            nav_files = config_dict.get("inputs", "nav_files")
            gal_nav_type = config_dict.get("model", "GAL", "nav_type")

            GnssDataManager.check_input_list("inputs.clk_files", clock_files, log)
            GnssDataManager.check_input_list("inputs.sp3_files", sp3_files, log)

            # not mandatory file checks
            GnssDataManager.check_input_list("inputs.ionex_files", ionex_files, log, warning=True)
            GnssDataManager.check_input_list("inputs.antex_files", antex_files, log, warning=True)
            GnssDataManager.check_input_list("inputs.nav_files", nav_files, log, warning=True)

            log.info("Launching GlobalIonoMap constructor (ionex products).")
            self.iono_gim.init(ionex_files, trace_dir)

            log.info("Reading Antenna Exchange Files.")
            for file in antex_files:
                AntexReader(file, self.get_data("phase_center"))

            for file in nav_files:
                RinexNavReader(file, self.get_data("nav_data"), gal_nav_type)

            log.info("Launching SatelliteClocks constructor (clocks from precise products).")
            self.sat_clocks.init(self.get_data("nav_data"), clock_files, True)

            log.info("Launching SatelliteOrbits constructor (orbits from precise products).")
            self.sat_orbits.init(self.get_data("nav_data"), sp3_files, True)

            bias_type_str = config_dict.get("inputs", "bias_type")
            log.info(f"Launching Satellite Code and Phase Bias Manager with bias {bias_type_str}.")
            if bias_type_str.upper() == "DCB":
                GnssDataManager.check_input_list("inputs.dcb_files", dcb_files, log)
                self.sat_bias.init(self.get_data("nav_data"), dcb_files, EnumSatelliteBias.DCB)
            elif bias_type_str.upper() == "OSB":
                GnssDataManager.check_input_list("inputs.dcb_files", osb_files, log)
                self.sat_bias.init(self.get_data("nav_data"), osb_files, EnumSatelliteBias.OSB)
            else:
                raise IOError(f"Unknown bias type in inputs.bias_type ({bias_type_str})")

        else:
            raise IOError(f"Unknown Model {gnss_alg}")

        if config_dict.get("inputs", "trace_files"):
            self._trace_files(trace_dir, gnss_alg)

    def _trace_files(self, trace_dir, gnss_alg: EnumAlgorithmPNT):
        inputs_dir = f"{trace_dir}\\inputs"
        try:
            os.makedirs(inputs_dir)
        except:
            raise IOError(f"Cannot create dir: {inputs_dir}")
        # trace data files
        with open(f"{inputs_dir}\\RawObservationData.txt", "w") as file:
            file.write(str(self.get_data("obs_data")))
        if gnss_alg == EnumAlgorithmPNT.SPS:
            with open(f"{inputs_dir}\\RawNavigationData.txt", "w") as file:
                file.write(str(self.get_data("nav_data")))
        elif gnss_alg == EnumAlgorithmPNT.PR_PPP:
            with open(f"{inputs_dir}\\PreciseClocks.txt", "w") as file:
                file.write(str(self.get_data("sat_clocks")))
            with open(f"{inputs_dir}\\PreciseOrbits.txt", "w") as file:
                file.write(str(self.get_data("sat_orbits")))
            with open(f"{inputs_dir}\\PreciseBiasProducts.txt", "w") as file:
                file.write(str(self.get_data("sat_bias")))
            with open(f"{inputs_dir}\\GlobalIonoMap.txt", "w") as file:
                file.write(str(self.get_data("iono_gim")))
            with open(f"{inputs_dir}\\PhaseCenterVariations.txt", "w") as file:
                file.write(str(self.get_data("phase_center")))

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

    @staticmethod
    def check_input_list(folder_name: str, input_list: list, log, warning=False):
        if len(input_list) == 0:
            if warning:
                log.warning(f"Optional input files {folder_name} not provided. "
                            f"Please check if they are necessary for some of the selected models.")
            else:
                raise IOError(f"Mandatory input files {folder_name} not provided. Please check the configurations.")
        if len(input_list) == 1 and not input_list[0]:
            if warning:
                log.warning(f"Optional input files {folder_name} not provided. "
                            f"Please check if they are necessary for some of the selected models.")
            else:
                raise IOError(f"Mandatory input files {folder_name} not provided. Please check the configurations.")
