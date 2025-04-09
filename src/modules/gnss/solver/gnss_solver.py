""" Module with the GNSS PNT Solver Algorithm """
import os

from src.data_mng.gnss.state_space import GnssStateSpace
from src.data_types.gnss import DataType
from src.modules.gnss.solver.lsq_engine import LSQ_Engine_Position, LSQ_Engine_Velocity
from src.modules.gnss.solver.obs_reconstructor import PseudorangeReconstructor, RangeRateReconstructor
from src.common_log import get_logger, GNSS_ALG_LOG
from src.errors import SolverError, ConfigError
from src.io.config import config_dict
from src.io.config.enums import *
from src.data_mng.gnss.geometry import SystemGeometry
from src import constants
from src.models.gnss_models import TropoManager, IonoManager


class GnssSolver:
    """
    GnssSolver class.
    Implements GNSS-based estimation of the receiver Position, Velocity Time (PVT) states.

    The module is able to process currently GPS and GAL constellations in both single-frequency and dual-frequency
    configurations.

    The estimation process is achieved in two parts:
        1) First, the receiver position and clock bias are estimated, using pseudorange observations
        2) Then, the receiver velocity and clock drift are estimated, using doppler (converted to pseudorange rate
            observations)

    The state vector (estimated parameters) is dependent on the chosen user configuration and simulation setup. The
    available parameters for estimation are (for more information see :py:class:`GnssStateSpace`):
        * Receiver position [3x1]: this state is mandatory and always estimated
        * Receiver clock bias [1]: this state is mandatory and always estimated
        * Inter system bias (ISB) [1] : this state is enabled for dual-constellation scenarios
        * Troposphere [1] (Zenith Wet Delay ZWD coefficient): this state is enabled/disabled in the user configurations
        * Ionosphere [1xN_i]: this state is enabled for dual-frequency scenarios, where N_i is the number of satellites
            for constellation i
        * Receiver velocity [3x1]: the velocity estimation is enabled/disabled by the user in the configurations and
            requires Doppler measurements
        * Clock drift [1xN_C]: the clock drift state is a by-product of the velocity estimation process, where N_C is
            the number of available constellations.

    Attributes:
        obs_data_for_pos(src.data_mng.gnss.ObservationData): observation data for the position estimation
            process (may contain raw, smooth or iono-free pseudorange observations).
        obs_data_for_vel(src.data_mng.gnss.ObservationData): observation data for the velocity estimation
            process (contains Doppler measurements).
        nav_data(src.data_mng.gnss.NavigationData): navigation data object containing RINEX NAV ephemerides.
        sat_orbits(src.data_mng.gnss.sat_orbit_data.SatelliteOrbits): `SatelliteOrbits` object with orbit data.
        sat_clocks(src.data_mng.gnss.sat_clock_data.SatelliteClocks): `SatelliteClocks` object with clock data.
        sat_bias(src.data_mng.gnss.bias_manager.BiasManager): `BiasManager` object that manages code and phase
            satellite biases.
        write_trace(bool): control flag to enable trace files for output
        trace_dir(str): path to store the trace files
        log(logging.Logger): logger instance
        solution(list[GnssStateSpace]): a list with the solved PNT states for each epoch (output)
    """

    def __init__(self, data_manager, trace_dir):
        """
        Constructor of the GnssSolver class.

        Args:
            data_manager(src.data_mng.gnss.gnss_data_mng.GnssDataManager): data manager with all necessary data for
                the GNSS run
            trace_dir(str): path to store the trace files
        """
        self.obs_data_for_pos = data_manager.get_clean_obs_data()
        self.obs_data_for_vel = data_manager.get_raw_obs_data()
        self.nav_data = data_manager.get_data("nav_data")
        self.sat_orbits = data_manager.get_data("sat_orbits")
        self.sat_clocks = data_manager.get_data("sat_clocks")
        self.sat_bias = data_manager.get_data("sat_bias")
        self.phase_center = data_manager.get_data("phase_center")
        self.write_trace = config_dict.get("solver", "trace_files")
        if self.write_trace:
            self.trace_dir = f"{trace_dir}\\solver"
            try:
                os.makedirs(self.trace_dir)
            except:
                raise IOError(f"Cannot create dir: {self.trace_dir}")
        else:
            self.trace_dir = None

        self.log = get_logger(GNSS_ALG_LOG)
        self.log.info("Starting module GNSS Positioning Solver...")

        # user configurations
        self._set_solver_metadata(config_dict, data_manager)

        # solution list
        self.solution = list()

    def _set_solver_metadata(self, config, data_manager):
        """
        Set up the metadata dict with the scenario setup (user configurations)

        Args:
            config(src.io.config.config.Config): user configurations
            data_manager(src.data_mng.gnss.gnss_data_mng.GnssDataManager): data manager with all necessary data for
                the GNSS run
        Raises:
            ConfigError: an exception is raised if there are problems/inconsistencies with the user configuration

        """
        # Fetching user options
        MAX_ITER = config.get("solver", "n_iterations")  # maximum number of iterations
        STOP_CRITERIA = config.get("solver", "stop_criteria")  # RMS threshold for stop criteria
        SOLVER_ALG = EnumSolver.init_model(config.get("solver", "algorithm"))
        TX_TIME_ALG = EnumTransmissionTime(config.get("solver", "transmission_time_alg"))
        REL_CORRECTION = EnumOnOff(config.get("solver", "relativistic_corrections"))  # 0 disable, 1 enable
        CONSTELLATIONS = config.get("model", "constellations")
        ELEVATION_FILTER = config.get("solver", "elevation_filter")
        VELOCITY_EST = config.get("model", "estimate_velocity")

        APRIORI_CONSTRAIN = EnumOnOff(config.get("solver", "a_priori_constrain"))
        FEEDFORWARD = EnumOnOff(config.get("solver", "feedforward_solution"))
        INITIAL_POS = config.get("solver", "initial_pos_cov")
        INITIAL_VEL = config.get("solver", "initial_vel_cov")
        INITIAL_CLOCK_BIAS = config.get("solver", "initial_clock_cov")
        INITIAL_ISB = config.get("solver", "initial_isb_cov")
        INITIAL_IONO = config.get("solver", "initial_iono_cov")
        INITIAL_TROPO = config.get("solver", "initial_tropo_cov")
        INITIAL_CLOCK_RATE = config.get("solver", "initial_clock_rate_cov")
        INITIAL_STATE = {
            "pos": INITIAL_POS,
            "vel": INITIAL_VEL,
            "clock": INITIAL_CLOCK_BIAS,
            "isb": INITIAL_ISB,
            "iono": INITIAL_IONO,
            "tropo": INITIAL_TROPO,
            "clock_rate": INITIAL_CLOCK_RATE
        }

        TROPO = TropoManager()
        obs_model = config_dict.get("obs_model")

        IONO = {}
        MODEL = {}
        CODES = {}
        DOPPLER = {}

        # Set up the user models for each active constellation
        for const in CONSTELLATIONS:
            IONO[const] = IonoManager(const, data_manager.get_data("iono_gim"))

            code_types = self.obs_data_for_pos.get_code_types(const)
            doppler_types = self.obs_data_for_vel.get_doppler_types(const)

            # check if the model is single frequency or dual frequency
            if len(code_types) == 2 and (obs_model == EnumObservationModel.UNCOMBINED):
                # Dual-Frequency Uncombined Model
                self.log.info(f"Selected model for {const} is Dual-Frequency Uncombined with observations {code_types}")
                MODEL[const] = EnumFrequencyModel.DUAL_FREQ
                CODES[const] = [code_types[0], code_types[1]]  # setting main and second code type
            elif len(code_types) == 1 and obs_model == EnumObservationModel.COMBINED:
                # Iono-Free Model
                self.log.info(f"Selected model for {const} is Iono-Free with code observations {code_types}")
                iono_code = code_types[0]
                if not DataType.is_iono_free_code(iono_code) and not DataType.is_iono_free_smooth_code(iono_code):
                    raise ConfigError(f"Iono-Free Model is selected for constellation {const} but no iono-free code "
                                      f"observations are available. Available code observations: {code_types}")
                MODEL[const] = EnumFrequencyModel.SINGLE_FREQ  # Iono-free is treated as single frequency in the LS
                CODES[const] = [iono_code]
                IONO[const].iono_model = EnumIonoModel.DISABLED  # disable Iono model in Iono-Free scenarios
            elif len(code_types) == 1:
                # Single Frequency Model
                self.log.info(f"Selected model for {const} is Single Frequency with code types {code_types[0]}")
                MODEL[const] = EnumFrequencyModel.SINGLE_FREQ
                CODES[const] = [code_types[0]]
            else:
                raise ConfigError(f"Unable to initialize GNSS Solver Model for constellation {const}. Check configs.")

            # Velocity estimation checks
            if VELOCITY_EST:
                self.log.info(f"Velocity estimation is enabled. Available Doppler observations for constellation "
                              f"{const} are {doppler_types}")
                if len(doppler_types) == 0:
                    self.log.error(f"No Doppler observations for constellation {const} are available")
                DOPPLER[const] = doppler_types

        # fill info dict
        self._metadata = {
            "CONSTELLATIONS": CONSTELLATIONS,
            "MAX_ITER": MAX_ITER,
            "STOP_CRITERIA": STOP_CRITERIA,
            "SOLVER": SOLVER_ALG,
            "APRIORI_CONSTRAIN": APRIORI_CONSTRAIN,
            "FEEDFORWARD": FEEDFORWARD,
            "TROPO": TROPO,
            "IONO": IONO,
            "MODEL": MODEL,
            "CODES": CODES,
            "DOPPLER": DOPPLER,
            "REL_CORRECTION": REL_CORRECTION,
            "TX_TIME_ALG": TX_TIME_ALG,
            "INITIAL_STATES": INITIAL_STATE,
            "ELEVATION_FILTER": ELEVATION_FILTER,
            "VELOCITY_EST": VELOCITY_EST
        }

    def solve(self) -> None:
        """ Launch GNSS Solver algorithm (main function) """
        # available epochs
        epochs = self.obs_data_for_pos.get_epochs()

        # iterate over all available epochs
        for epoch in epochs:
            self.log.info(f"processing epoch {epoch}...")
            # fetch observation data for this epoch
            obs_for_epoch = self.obs_data_for_pos.get_epoch_data(epoch)

            # initialize solve-for variables (receiver position and bias) for the present epoch
            state = self._init_state(epoch, obs_for_epoch.get_satellites())

            # call lower level of position estimation
            success = self._estimate_position(epoch, obs_for_epoch, state)

            if success:
                # add solution to Output timeseries
                self.log.info(f"Successfully solved positioning for epoch {str(epoch)} with "
                              f"RMS = {state.get_additional_info('rms')} [m]")

                if self._metadata["VELOCITY_EST"]:
                    success = self._estimate_velocity(state, epoch)

                    if success:
                        self.log.info(f"Successfully solved velocity estimation for epoch {str(epoch)} with "
                                      f"RMS = {state.get_additional_info('vel_rms')} [m/s]")
                    else:
                        self.log.warning(f"No velocity solution computed for epoch {str(epoch)}.")

                # store data for this epoch
                self.solution.append(state)

            else:
                self.log.warning(f"No solution computed for epoch {str(epoch)}.")

        self.log.info("Successfully ending module GNSS Positioning Solver...")

    def _init_state(self, epoch, sat_list):
        """ Initialize state vector for this epoch.

        Args:
            epoch(src.data_types.date.Epoch): epoch to initialize the state vector
            sat_list(list[src.data_types.gnss.Satellite]) : list of available satellites
        Returns:
            GnssStateSpace : initialized state vector
        """
        # initialize GNSS state from user input configurations
        feedforward = self._metadata["FEEDFORWARD"]
        if len(self.solution) == 0 or feedforward == EnumOnOff.DISABLED:
            state = GnssStateSpace(self._metadata, epoch, sat_list)
            state.initial_state = state.clone()
        else:
            # initialize from previous state
            prev_state = self.solution[-1]
            state = prev_state.clone()  # deep copy
            state.epoch = epoch
            state.initial_state = state.clone()
        return state

    def _estimate_position(self, epoch, obs_data, state):
        """
        Sub function to launch the position estimation procedure

        Args:
            epoch(src.data_types.date.Epoch): epoch to process
            obs_data (src.data_mng.gnss.observation_data.EpochData) : instance of `EpochData` (GNSS observable database
                for a single epoch)
            state(GnssStateSpace) : state vector to process

        Returns:
            bool : True if the process succeeds and False if it fails. In that case debug information is written in the
                log file
        """
        # begin iterative process for this epoch
        iteration = 0
        rms = rms_prev = 1
        prefit_residuals = postfit_residuals = dop_matrix = None
        apply_elevation_filter = False if self._metadata["ELEVATION_FILTER"] == -1 else True
        sats_for_epoch = obs_data.get_satellites()

        # build system geometry for this epoch
        system_geometry = SystemGeometry(obs_data, self.sat_clocks, self.sat_orbits, self.phase_center)

        self.log.debug(f"Available Satellites: {sats_for_epoch}")

        # Iterated Least-Squares algorithm
        while iteration < self._metadata["MAX_ITER"]:
            # compute geometry-related data for each satellite link
            system_geometry.compute(epoch, state, self._metadata)

            # solve the Least Squares
            try:
                prefit_residuals, postfit_residuals, dop_matrix, rms = \
                    self._iterate_pos(str(iteration), system_geometry, obs_data, state, epoch)
            except SolverError as e:
                self.log.warning(f"Least Squares failed for {str(epoch)} on iteration {iteration}."
                                 f"Reason: {e}")
                return False

            # check stop condition
            if self._stop(rms_prev, rms, self._metadata["STOP_CRITERIA"]):
                self.log.debug(f"Least Squares was successful. Reached convergence at iteration {iteration}")
                break

            # increase iteration counter
            rms_prev = rms
            iteration += 1

        # end of iterative procedure
        if iteration == self._metadata["MAX_ITER"]:
            self.log.warning(f"PVT failed to converge for epoch {str(epoch)}, with RMS={rms}. "
                             f"No solution will be computed for this epoch.")
            return False

        # perform additional iteration after elevation filter
        if apply_elevation_filter:
            self.log.info(f"Applying elevation filter with threshold {self._metadata['ELEVATION_FILTER']} [deg] "
                          f"and performing additional iteration")
            self._elevation_filter(system_geometry, self._metadata["ELEVATION_FILTER"])

            # solve the Least Squares
            try:
                prefit_residuals, postfit_residuals, dop_matrix, rms = \
                    self._iterate_pos("final", system_geometry, obs_data, state, epoch)
            except SolverError as e:
                self.log.warning(f"Least Squares failed for {str(epoch)} on additional iteration "
                                 f"(after elevation filter). Reason: {e}")
                return False

        # save other iteration data to state variable
        state.add_additional_info("geometry", system_geometry)
        state.add_additional_info("pr_prefit_residuals", prefit_residuals)
        state.add_additional_info("pr_postfit_residuals", postfit_residuals)
        state.add_additional_info("rms", rms)
        state.add_additional_info("dop_matrix", dop_matrix)

        return True

    def _iterate_pos(self, iteration, system_geometry, obs_data, state, epoch):
        """ Low-level function to solve a single iteration of the position estimation iterative process """
        self._check_model_availability(system_geometry, epoch)

        satellite_list = system_geometry.get_satellites()
        state.build_index_map(satellite_list)

        if self.trace_dir is not None:
            trace_file = f"{self.trace_dir}\\PseudorangeReconstructionIter_{iteration}.txt"
        else:
            trace_file = None
        # TODO: move reconstructor creation to LSQ_Engine. Add Carrier Phase Reconstructor
        reconstructor = PseudorangeReconstructor(system_geometry, self._metadata, state, self.sat_bias, trace_file)

        # build LSQ Engine matrices for all satellites
        lsq_engine = LSQ_Engine_Position(satellite_list, self._metadata, epoch, obs_data, reconstructor, state)

        # solve LS problem
        return lsq_engine.solve_ls(state)

    def _check_model_availability(self, system_geometry, epoch) -> None:
        """
        Checks if it is possible to perform the dual frequency / single frequency PVT estimation, or, in contrast,
        we have not enough data to perform the computations.

        Raises:
            SolverError : raises an exception there is not enough information to process this epoch
        """
        MIN_SAT = 4  # default: estimate position + clock (dimension 4)
        if len(self._metadata["CONSTELLATIONS"]) != 1:
            MIN_SAT += 1  # estimate ISB
        if self._metadata["TROPO"].estimate_tropo():
            MIN_SAT += 1  # estimate tropo
        sat_list = system_geometry.get_satellites()

        # check number of available satellites
        if len(sat_list) < MIN_SAT:
            raise SolverError(
                f"Not enough satellites to compute model PVT positioning at {str(epoch)}. "
                f"Available satellites with data for provided codes: "
                f"{sat_list}. Minimum number of satellites is {MIN_SAT}. Aborting the solution for this epoch...")

    def _elevation_filter(self, system_geometry, cutoff: float):
        """ Apply elevation filter to remove satellites below the provided cutoff elevation """
        # get elevation threshold in radians
        sats_to_remove = []
        for sat, sat_info in system_geometry.items():

            if sat_info.el * constants.RAD2DEG < cutoff:
                sats_to_remove.append(sat)
                self.log.debug(f"Removing satellite {sat} due to elevation filter. Minimum threshold is {cutoff} "
                               f"[deg], computed elevation is {sat_info.el * constants.RAD2DEG} [deg]")

        # remove flagged satellites
        for sat in sats_to_remove:
            system_geometry.remove(sat)

        if sats_to_remove:
            self.log.debug(f"Removing satellites {sats_to_remove} due to elevation filter. ")

    def _iterate_vel(self, system_geometry, obs_data, state, epoch):
        """ Low-level function to solve a single iteration of the velocity estimation process """
        satellite_list = system_geometry.get_satellites()

        if self.trace_dir is not None:
            trace_file = f"{self.trace_dir}\\RangeRateReconstructor.txt"
        else:
            trace_file = None
        reconstructor = RangeRateReconstructor(system_geometry,
                                               self._metadata,
                                               state, trace_file)

        # build LSQ Engine matrices for all satellites
        lsq_engine = LSQ_Engine_Velocity(satellite_list, self._metadata, epoch, obs_data, reconstructor)

        # solve LS problem
        return lsq_engine.solve_ls(state)

    def _estimate_velocity(self, state, epoch):
        """
        Sub function to launch the velocity estimation procedure

        Args:
            epoch(src.data_types.date.Epoch): epoch to process
            state(GnssStateSpace) : state vector to process

        Returns:
            bool : True if the process succeeds and False if it fails. In that case debug information is written in the
                log file
        """
        # get observations and system geometry for this epoch
        obs_for_epoch = self.obs_data_for_vel.get_epoch_data(epoch)
        system_geometry = state.get_additional_info("geometry")

        self.log.info(f"Applying velocity estimation {str(epoch)}")
        self.log.debug(f"Available Satellites: {system_geometry.get_satellites()}")

        # Least-Squares algorithm
        try:
            prefit_residuals, postfit_residuals, _, rms = \
                self._iterate_vel(system_geometry, obs_for_epoch, state, epoch)
        except SolverError as e:
            self.log.warning(f"Least Squares failed for {str(epoch)}. Reason: {e}")
            return False

        # save other iteration data to state variable
        state.add_additional_info("pr_rate_prefit_residuals", prefit_residuals)
        state.add_additional_info("pr_rate_postfit_residuals", postfit_residuals)
        state.add_additional_info("vel_rms", rms)

        return True

    @staticmethod
    def _stop(rms_old, rms_new, stop_criteria):
        """ Apply stop condition for the iteration procedure """
        return abs((rms_old - rms_new) / rms_old) <= stop_criteria
