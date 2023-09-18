from numpy.linalg import norm
import numpy as np

from src.algorithms.gnss.estimators.state_space import GnssStateSpace
from src.algorithms.gnss.estimators.weighted_ls import WeightedLeastSquares
from src.algorithms.gnss.solver.obs_reconstructor import ObservationReconstruction
from src.common_log import get_logger
from src.errors import PVTComputationFail
from src.io.config import config_dict
from src.io.config.enums import *
from src.models.observation.geometry import SystemGeometry
from src import constants

np.set_printoptions(linewidth=np.inf)


def get_weight(system_geometry, sat):
    # TODO: need to add here the user defined sigmas as a multiplication factor
    # "obs_std", and can add the possibility of this mask as well.
    sigma_elevation = np.e ** (-system_geometry.get("el", sat))
    w = (1 / sigma_elevation) ** 2

    return w


class GnssSolver:
    """
        GPSSolver. Implements GPS Position Velocity Time (PVT) algorithms, to compute receiver position and clock bias
        (not yet velocity).

        The following algorithms are included:

            -> Single constellation Single Point Positioning (SPP). This is the most basic GNSS algorithm. It refers to
                an epoch-wise iterated Least Squares (LS) parameter adjustment, where the observation equations are
                linearized and solved with respect to the receiver position and clock bias.
                The linearized equation to be solved in a LS sense is:
                    dy = G @ dx (see Eq. 6.9 of **REF[1]**), where dy is the observation residual, G is the design
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

                The GPS observation model comprises:
                    * true range
                    * receiver clock bias
                    * satellite clock bias (corrected for relativistic corrections and TGDs)
                    * ionosphere - a priori Klobuchar Ionospheric Model or estimated in DF mode
                    * troposphere - a priori Saastamoinen model


        Iterated Least-Squares PVT - High Level Algorithm
        -------

        for each observation epoch (time tag in rinex obs):
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
                compute prefit_residual = observation - predicted_observation
                set up system LS system and solve for the state variables dx, dy, dz, dt, ...

                update (x,y,z) += (dx,dy,dz) and finish this iteration
                finish iterative procedure if process converged

    """

    def __init__(self, obs_data, nav_data):
        """
        Args:
            obs_data : observation data
            nav_data : navigation data
        """
        self.obs_data = obs_data
        self.nav_data = nav_data

        self.log = get_logger("GNSS_SOLVER")
        self.log.info("Starting module GNSS Positioning Solver...")

        # user configurations
        self._set_solver_info(config_dict)

        # solution dict
        self.solution = []

    def _set_solver_info(self, config):

        # TODO : add this to the config object rather than being here
        # Fetching user options
        MAX_ITER = config.get("solver", "n_iterations")  # maximum number of iterations
        STOP_CRITERIA = config.get("solver", "stop_criteria")  # RMS threshold for stop criteria
        SOLVER = config.get("solver", "algorithm")  # 0 - LS, 1 - WLS
        TX_TIME_ALG = config.get("solver", "transmission_time_alg")
        REL_CORRECTION = EnumOnOff(config.get("solver", "relativistic_corrections"))  # 0 disable, 1 enable
        INITIAL_POS = config.get("solver", "initial_pos_std")
        INITIAL_CLOCK_BIAS = config.get("solver", "initial_clock_std")

        TROPO = {}
        IONO = {}
        MODEL = {}
        CODES = {}

        for const in self.obs_data.get_constellations():
            TROPO[const] = EnumTropo(config.get("model", const, "troposphere"))  # 0 - no model, 1 - Saastamoinen
            IONO[const] = EnumIono(config.get("model", const, "ionosphere"))

            # check if the model is single frequency or dual frequency
            code_types = self.obs_data.get_code_types(const)
            if len(code_types) > 1 and not config_dict.is_iono_free():
                # Dual-Frequency Model
                self.log.info(f"Selected model for {const} is Dual-Frequency Uncombined with observations {code_types}")
                MODEL[const] = EnumModel.DUAL_FREQ
                CODES[const] = [code_types[0], code_types[1]]  # setting main and second code type
            elif config_dict.is_iono_free():
                # Iono-Free Model
                self.log.info(f"Selected model for {const} is Iono-Free with observation {code_types}")
                MODEL[const] = EnumModel.SINGLE_FREQ  # we process as single frequency
                CODES[const] = [code_types[0], None]
            else:
                # Single Frequency Model
                self.log.info(f"Selected model for {const} is Single Frequency with observation {code_types[0]}")
                MODEL[const] = EnumModel.SINGLE_FREQ
                CODES[const] = [code_types[0], None]

        # fill info dict
        self._info = {
            "MAX_ITER": MAX_ITER,
            "STOP_CRITERIA": STOP_CRITERIA,
            "SOLVER": EnumSolver(SOLVER),
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
            # fetch observation data for this epoch
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
            initial_pos = self._info["INITIAL_POS"][0:3]
            initial_clock = self._info["INITIAL_CLOCK_BIAS"][0]
            state = GnssStateSpace(position=np.array(initial_pos, dtype=np.float64),
                                   clock_bias=initial_clock,
                                   epoch=epoch,
                                   sat_list=sat_list)
        else:
            prev_state = self.solution[-1]
            state = GnssStateSpace(position=prev_state.position.copy(),
                                   clock_bias=prev_state.clock_bias,
                                   epoch=epoch,
                                   sat_list=sat_list,
                                   iono=prev_state.iono)
        return state

    @staticmethod
    def _stop(rms_old, rms_new, stop_criteria):
        return abs((rms_old - rms_new) / rms_old) <= stop_criteria

    def _solve_for_epoch(self, epoch, obs_data, state):
        # begin iterative process for this epoch
        iteration = 0
        success = False
        rms = rms_prev = 1
        geometry_matrix = cov = prefit_residuals = postfit_residuals = None

        # TODO: tmp -> this needs to be revisited
        gps_codes = self._info['CODES']['GPS']
        # gal_codes

        # build system geometry for this epoch
        system_geometry = SystemGeometry(self.nav_data, obs_data)

        # check data availability for this epoch
        data_available = self._check_model_availability(system_geometry, obs_data, epoch)

        if not data_available:
            # the log message is issued in self._check_model_availability
            return False  # not enough data to process this epoch

        self.log.info(f"Processing epoch {str(epoch)}")
        self.log.debug(f"Available Satellites: {system_geometry.get_satellites()}")

        # Iterated Least-Squares algorithm
        # Note: in the eventuality of implementing other filters, e.g.: Kalman Filter, then we need to branch
        # the solver from here...
        while iteration < self._info["MAX_ITER"]:
            # compute geometry-related data for each satellite link
            system_geometry.compute(epoch, state.position, state.clock_bias, self._info["TX_TIME_ALG"],
                                    gps_codes[0], self._info["REL_CORRECTION"])

            # solve the Least Squares
            try:
                postfit_residuals, geometry_matrix, prefit_residuals, cov = \
                    self._compute(system_geometry, obs_data, state, epoch)
            except PVTComputationFail as e:
                self.log.warning(f"Least Squares failed for {str(epoch)} on iteration {iteration}."
                                 f"Reason: {e}")
                break

            # update RMS for this iteration
            rms = np.linalg.norm(postfit_residuals)

            # check stop condition
            if self._stop(rms_prev, rms, self._info["STOP_CRITERIA"]):
                self.log.debug(f"Least Squares was successful. Reached convergence at iteration {iteration}")
                success = True
                break

            # increase iteration counter
            rms_prev = rms
            iteration += 1

        # end of iterative procedure
        if iteration == self._info["MAX_ITER"]:
            self.log.warning(f"PVT failed to converge for epoch {str(epoch)}, with RMS={rms}. "
                             f"No solution will be computed for this epoch.")
            return False

        # save other iteration data to state variable
        state.add_additional_info("geometry", system_geometry)
        state.add_additional_info("geometry_matrix", geometry_matrix)
        state.add_additional_info("prefit_residuals", prefit_residuals)
        state.add_additional_info("postfit_residuals", postfit_residuals)
        state.add_additional_info("rms", rms)
        state.add_additional_info("cov", cov)

        return success

    def _initialize_ls(self, n_sats):
        # initializations
        gps_codes = self._info['CODES']['GPS']

        # least squares arrays
        # TODO: Currently this is GPS specific
        if self._info["MODEL"]["GPS"] == EnumModel.SINGLE_FREQ:
            y = np.zeros(n_sats)  # observation vector <=> prefit residuals
            design = np.ones((n_sats, 4))  # geometry matrix
            weight = np.eye(n_sats)  # diagonal weight matrix
        else:
            y = np.zeros(n_sats * 2)  # observation vector <=> prefit residuals
            geometry = np.ones((n_sats, 4))  # geometry matrix (state + clock)
            weight = np.eye(n_sats * 2)  # diagonal weight matrix
            ionoMatrix1 = np.eye(n_sats)
            ionoMatrix2 = np.eye(n_sats) * (gps_codes[0].freq.freq_value / gps_codes[1].freq.freq_value) ** 2

            design = np.block([
                [geometry, ionoMatrix1],
                [geometry, ionoMatrix2]]
            )
        return y, design, weight

    def _solve_ls(self, y, G, W, state, n_sats):
        # solve LS problem for this iteration
        try:
            solver = WeightedLeastSquares(y, G, W=W)
            solver.solve()

            cov = np.linalg.inv(G.T @ W @ G)  # covariance matrix of the LS estimator

        except (AttributeError, np.linalg.LinAlgError) as e:
            # possible error in the numpy.linalg.inv() function -> solution not possible
            raise PVTComputationFail(e)

        dX = solver.get_solution()

        # update state vector with incremental dX
        state.position += dX[0:3]
        state.clock_bias = dX[3] / constants.SPEED_OF_LIGHT  # receiver clock in seconds

        # if iono is estimated
        if self._info["MODEL"]["GPS"] == EnumModel.DUAL_FREQ:
            state.iono = [dX[i + 4] for i in range(n_sats)]

        # if isb is estimated
        # ...

        # get post-fit residuals
        post_fit = y - G[:, 0:3] @ dX[0:3]

        return post_fit, cov

    def _compute(self, system_geometry, obs_data, state, epoch):
        satellite_list = system_geometry.get_satellites()
        n_sats = len(satellite_list)  # number of available satellites

        y_vec, design_mat, weight_mat = self._initialize_ls(n_sats)

        datatypes = self._info["CODES"]["GPS"]
        reconstructor = ObservationReconstruction(system_geometry,
                                                  self._info["TROPO"],
                                                  self._info["IONO"],
                                                  self.nav_data.header,
                                                  self._info["REL_CORRECTION"])

        iFreq = 0
        for datatype in datatypes:

            iSat = 0
            for sat in satellite_list:

                # get observable
                obs = obs_data.get_observable(sat, datatype)  # main observable

                # fetch valid navigation message (closest to the current epoch)
                nav_message = self.nav_data.get_sat_data_for_epoch(sat, epoch)

                # compute predicted observation
                predicted_obs = reconstructor.compute(nav_message, sat, epoch, datatype)

                # prefit residuals (measured observation - predicted observation)
                prefit_residuals = obs - predicted_obs

                # get LOS vector w.r.t. ECEF frame (column in geometry matrix)
                line_sight = system_geometry.get_unit_line_of_sight(sat)

                # filling the corresponding entry of data vectors for the Least Squares quantities
                y_vec[iFreq * n_sats + iSat] = prefit_residuals.value
                design_mat[iFreq * n_sats + iSat][0:3] = line_sight

                # Weight matrix -> sigma = 1 / e^{-elevation}
                weight_mat[iFreq * n_sats + iSat][iFreq * n_sats + iSat] = get_weight(system_geometry, sat)

                iSat += 1
            iFreq += 1

        # solve LS problem
        post_fit, cov = self._solve_ls(y_vec, design_mat, weight_mat, state, n_sats)

        return post_fit, design_mat, y_vec, cov

    def _check_model_availability(self, system_geometry, obs_data, epoch):
        """
        TODO: this docs need to be revised
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
        # TODO: currently only GPS is processed. This will need to be revisited
        # TODO: esta é a ordem certa para dual frequency? se calhar o check tem de ser depois de eliminar os dados
        model = self._info["MODEL"]["GPS"]
        codes = self._info["CODES"]["GPS"]
        MIN_SAT = 4

        sat_list = obs_data.get_sats_for_datatypes(codes)
        if len(sat_list) < MIN_SAT:
            self.log.warning(
                f"Not enough satellites to compute model {model} PVT positioning at {str(epoch)}. "
                f"Available satellites with data for codes {codes}: "
                f"{sat_list}. Minimum number of satellites is {MIN_SAT}. Aborting the solution for this epoch...")
            return False

        # Dual Frequency model
        if model == EnumModel.DUAL_FREQ:
            # removing satellites with data for only 1 frequency
            removed = []
            for sat in obs_data.get_satellites():
                if sat not in sat_list:
                    system_geometry.remove(sat)
                    removed.append(sat)
            if removed:
                self.log.info(f"Removing satellites {removed} at epoch {str(epoch)} due to "
                              f"inconsistencies in code data for both frequencies")
        return True

    def _elevation_filter(self, system_geometry):
        # TODO: re-add this on last iteration

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
