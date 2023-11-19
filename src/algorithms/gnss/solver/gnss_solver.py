import numpy as np

from src.algorithms.gnss.estimators.state_space import GnssStateSpace
from src.algorithms.gnss.solver.lsq_engine import LSQ_Engine
from src.algorithms.gnss.solver.obs_reconstructor import ObservationReconstruction
from src.common_log import get_logger
from src.errors import PVTComputationFail, ConfigError
from src.io.config import config_dict
from src.io.config.enums import *
from src.models.gnss_obs.geometry import SystemGeometry
from src import constants

np.set_printoptions(linewidth=np.inf)


class GnssSolver:
    """
        GPSSolver. Implements GPS Position Velocity Time (PVT) algorithms, to compute receiver position and clock bias
        (not yet velocity).

        The following algorithms are included:

            -> Single constellation Single Point Positioning (SPP). This is the most basic GNSS algorithm. It refers to
                an epoch-wise iterated Least Squares (LS) parameter adjustment, where the gnss_obs equations are
                linearized and solved with respect to the receiver position and clock bias.
                The linearized equation to be solved in a LS sense is:
                    dy = G @ dx (see Eq. 6.9 of **REF[1]**), where dy is the gnss_obs residual, G is the design
                    (geometry) matrix, and dx is the solve-for state vector x = [dx, dy, dz, dt, dl]
                    After x is computed for iteration i, the receiver position is updated as:
                        X(i) = X(i-1) + dx
                    The process is iterated in i until convergence. Then, the process restarts for the next epoch
                    evaluation. The state vector can also comprehend other parameters (dl), for example ionosphere
                    delays.

                There exist a couple of variations to this solution, depending on the user configurations.
                The implemented variations are:
                    - Single Frequency (SF):
                        * raw pseudorange observables (e.g.: C1 pseudorange)
                        * pseudorange smoothed with Hatch filter (e.g.: SPR1)
                        * iono free combination, using two frequencies to obtain, for example, C12 pseudorange
                        * smooth iono free observables (e.g.: SPR12)

                    - Dual Frequency (DF)
                        * in DF processing, the system of equations comprises observations from two frequencies
                            (e.g.: C1 and C2)

                The GPS gnss_obs model comprises:
                    * true range
                    * receiver clock bias
                    * satellite clock bias (corrected for relativistic corrections and TGDs)
                    * ionosphere - a priori Klobuchar Ionospheric Model or estimated in DF mode
                    * troposphere - a priori Saastamoinen model


        Iterated Least-Squares PVT - High Level Algorithm
        -------

        for each gnss_obs epoch (time tag in rinex obs):
            set initial pos = 0, dt = 0 (or with previous results...)

            get closest navigation message

            for each iteration:
                set r = 0, dt = 0 for first iteration or with results of previous iteration

                for each satellite that there exist measurements:
                    compute transit time
                    compute signal transmission epoch
                    with these, compute final satellite position at TX, and rotate to ECEF frame at t = RX
                    get true range
                    estimate azimuth, elevation
                    compute dt_sat and implement TGD
                    compute iono / tropo delays

                compute predicted observations using the quantities computed for this iteration:
                compute prefit_residual = gnss_obs - predicted_observation
                set up system LS system and solve for the state variables dx, dy, dz, dt, ...

                update (x,y,z) += (dx,dy,dz) and finish this iteration
                finish iterative procedure if process converged

    """

    def __init__(self, obs_data, nav_data):
        """
        Args:
            obs_data : gnss_obs data
            nav_data : navigation data
        """
        self.obs_data = obs_data
        self.nav_data = nav_data

        self.log = get_logger("GNSS_SOLVER")
        self.log.info("Starting module GNSS Positioning Solver...")

        # user configurations
        self._set_solver_metadata(config_dict)

        # solution dict
        self.solution = []

    def _set_solver_metadata(self, config):

        # TODO : add this to the config object rather than being here
        # Fetching user options
        MAX_ITER = config.get("solver", "n_iterations")  # maximum number of iterations
        STOP_CRITERIA = config.get("solver", "stop_criteria")  # RMS threshold for stop criteria
        SOLVER = EnumSolver.init_model(config.get("solver", "algorithm"))
        TX_TIME_ALG = config.get("solver", "transmission_time_alg")
        REL_CORRECTION = EnumOnOff(config.get("solver", "relativistic_corrections"))  # 0 disable, 1 enable
        INITIAL_POS = config.get("solver", "initial_pos_std")
        INITIAL_CLOCK_BIAS = config.get("solver", "initial_clock_std")
        CONSTELLATIONS = config.get("model", "constellations")
        _model = config.get_model()

        TROPO = {}
        IONO = {}
        MODEL = {}
        CODES = {}

        for const in CONSTELLATIONS:
            TROPO[const] = EnumTropo.init_model(config.get("model", const, "troposphere"))
            IONO[const] = EnumIono.init_model(config.get("model", const, "ionosphere"))

            # check if the model is single frequency or dual frequency
            code_types = self.obs_data.get_code_types(const)
            if len(code_types) == 2 and (_model != EnumPositioningMode.SPS_IF):
                # Dual-Frequency Model
                self.log.info(f"Selected model for {const} is Dual-Frequency Uncombined with observations {code_types}")
                MODEL[const] = EnumModel.DUAL_FREQ
                CODES[const] = [code_types[0], code_types[1]]  # setting main and second code type
            elif len(code_types) == 1 and _model == EnumPositioningMode.SPS_IF:
                # Iono-Free Model
                self.log.info(f"Selected model for {const} is Iono-Free with gnss_obs {code_types}")
                MODEL[const] = EnumModel.SINGLE_FREQ  # we process as single frequency
                CODES[const] = [code_types[0]]
                IONO[const] = EnumIono.DISABLED
            elif len(code_types) == 1:
                # Single Frequency Model
                self.log.info(f"Selected model for {const} is Single Frequency with gnss_obs {code_types[0]}")
                MODEL[const] = EnumModel.SINGLE_FREQ
                CODES[const] = [code_types[0]]
            else:
                raise ConfigError(f"Unable to initialize GNSS Solver Model for constellation {const}. Check configs.")

        # fill info dict
        self._metadata = {
            "CONSTELLATIONS": CONSTELLATIONS,
            "MAX_ITER": MAX_ITER,
            "STOP_CRITERIA": STOP_CRITERIA,
            "SOLVER": SOLVER,
            "TROPO": TROPO,
            "IONO": IONO,
            "MODEL": MODEL,
            "CODES": CODES,
            "REL_CORRECTION": REL_CORRECTION,
            "TX_TIME_ALG": EnumTransmissionTime(TX_TIME_ALG),
            "INITIAL_POS": INITIAL_POS,
            "INITIAL_CLOCK_BIAS": INITIAL_CLOCK_BIAS
        }

    def solve(self):
        # available epochs
        epochs = self.obs_data.get_epochs()

        # iterate over all available epochs
        for epoch in epochs:
            # fetch gnss_obs data for this epoch
            obs_for_epoch = self.obs_data.get_epoch_data(epoch)
            sats_for_epoch = self.obs_data.get_satellites()

            # initialize solve-for variables (receiver position and bias) for the present epoch
            state = self._init_state(epoch, sats_for_epoch)

            # call lower level of solve
            success = self._solve_for_epoch(epoch, obs_for_epoch, state)

            if success:
                # add solution to Output timeseries
                self.log.info(f"Successfully solved positioning for epoch {str(epoch)} with "
                              f"RMS = {state.get_additional_info('rms')} [m]")

                # store data for this epoch
                self.solution.append(state)

            else:
                self.log.warning(f"No solution computed for epoch {str(epoch)}.")

        self.log.info("Successfully ending module GNSS Positioning Solver...")

    def _init_state(self, epoch, sat_list):
        # initialize GNSS state
        if len(self.solution) == 0:
            position = np.array(self._metadata["INITIAL_POS"][0:3], dtype=np.float64)
            clock = self._metadata["INITIAL_CLOCK_BIAS"][0]
            state = GnssStateSpace(self._metadata,
                                   position=position,
                                   clock_bias=clock,
                                   epoch=epoch,
                                   sat_list=sat_list)
        else:
            # initialize from previous state
            prev_state = self.solution[-1]
            state = prev_state.clone()
            state.epoch = epoch
        return state

    @staticmethod
    def _stop(rms_old, rms_new, stop_criteria):
        return abs((rms_old - rms_new) / rms_old) <= stop_criteria

    def _solve_for_epoch(self, epoch, obs_data, state):
        # begin iterative process for this epoch
        iteration = 0
        success = False
        rms = rms_prev = 1
        prefit_residuals = postfit_residuals = dop_matrix = None

        # build system geometry for this epoch
        system_geometry = SystemGeometry(self.nav_data, obs_data)

        # check data availability for this epoch
        if not self._check_model_availability(system_geometry, obs_data, epoch):
            return False  # not enough data to process this epoch

        self.log.info(f"Processing epoch {str(epoch)}")
        self.log.debug(f"Available Satellites: {system_geometry.get_satellites()}")

        # Iterated Least-Squares algorithm
        while iteration < self._metadata["MAX_ITER"]:
            # compute geometry-related data for each satellite link
            system_geometry.compute(epoch, state, self._metadata)

            # solve the Least Squares
            try:
                postfit_residuals, prefit_residuals, dop_matrix, rms = \
                    self._compute(system_geometry, obs_data, state, epoch)
            except PVTComputationFail as e:
                self.log.warning(f"Least Squares failed for {str(epoch)} on iteration {iteration}."
                                 f"Reason: {e}")
                break

            # check stop condition
            if self._stop(rms_prev, rms, self._metadata["STOP_CRITERIA"]):
                self.log.debug(f"Least Squares was successful. Reached convergence at iteration {iteration}")
                success = True
                break

            # increase iteration counter
            rms_prev = rms
            iteration += 1

        # end of iterative procedure
        if iteration == self._metadata["MAX_ITER"]:
            self.log.warning(f"PVT failed to converge for epoch {str(epoch)}, with RMS={rms}. "
                             f"No solution will be computed for this epoch.")
            return False

        # save other iteration data to state variable
        state.add_additional_info("geometry", system_geometry)
        state.add_additional_info("prefit_residuals", prefit_residuals)
        state.add_additional_info("postfit_residuals", postfit_residuals)
        state.add_additional_info("rms", rms)
        state.add_additional_info("dop_matrix", dop_matrix)

        return success

    def _compute(self, system_geometry, obs_data, state, epoch):
        satellite_list = system_geometry.get_satellites()

        lsq_engine = LSQ_Engine(satellite_list, self._metadata)

        reconstructor = ObservationReconstruction(system_geometry,
                                                  self._metadata,
                                                  state,
                                                  self.nav_data.header)

        # build LSQ Engine matrices for all satellites
        lsq_engine.build_lsq(epoch, obs_data, reconstructor, self.nav_data)

        # solve LS problem
        return lsq_engine.solve_ls(state)

    def _check_model_availability(self, system_geometry, obs_data, epoch):
        """
        Checks if it is possible to perform the dual frequency / single frequency PVT estimation, or, in contrast,
        we have not enough data to perform the computations.

        If the user selected Dual Frequency model this functions:
            * checks if there are at least 4 satellites with dual frequency data
            * if that is the case, then eliminate satellites with single frequency data only (will be discarded)
            * if that is not the case, then fall back to the single frequency case

        If the user selected the Single frequency case (or we fell back to this model):
            * checks if, for the main code datatype, we have at least 4 satellites available
            * if that is not the case, a PVT solution is impossible to be computed -> Warn the user..

        """
        single_const = len(self._metadata["CONSTELLATIONS"]) == 1  # true if only one constellation
        MIN_SAT = 4 if single_const else 5
        sat_list = []

        for const in self._metadata["CONSTELLATIONS"]:
            model = self._metadata["MODEL"][const]
            codes = self._metadata["CODES"][const]
            sat_list_ = obs_data.get_sats_for_datatypes(codes)

            # Dual Frequency model
            if model == EnumModel.DUAL_FREQ:
                # removing satellites with data for only 1 frequency
                removed = []
                for sat in obs_data.get_satellites():
                    if sat not in sat_list_:
                        system_geometry.remove(sat)
                        removed.append(sat)
                if removed:
                    self.log.info(f"Removing satellites {removed} at epoch {str(epoch)} due to "
                                  f"inconsistencies in code data for both frequencies")

            sat_list += system_geometry.get_satellites()

        # check number of available satellites
        if len(sat_list) < MIN_SAT:
            self.log.warning(
                f"Not enough satellites to compute model PVT positioning at {str(epoch)}. "
                f"Available satellites with data for provided codes: "
                f"{sat_list}. Minimum number of satellites is {MIN_SAT}. Aborting the solution for this epoch...")
            return False

        return True

    def _elevation_filter(self, system_geometry):
        # TODO re-add this on last iteration

        # get elevation threshold in radians
        sats_to_remove = []
        threshold = self._info["ELEVATION_FILTER"]
        for sat, sat_info in system_geometry.items():

            if sat_info.el * constants.RAD2DEG < threshold:
                sats_to_remove.append(sat)
                # self.log.debug(f"Removing satellite {sat} in iteration {iteration} due to elevation filter. "
                #               f"Minimum threshold is {threshold} "
                #               f"[deg], computed elevation is {sat_info.el * Constant.RAD2DEG} [deg]")

        # remove flagged satellites
        for sat in sats_to_remove:
            system_geometry.remove(sat)

        if sats_to_remove:
            self.log.debug(f"Removing satellites {sats_to_remove} in due to "
                           f"elevation filter. ")
