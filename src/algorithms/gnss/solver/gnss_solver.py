from numpy.linalg import norm
import numpy as np

from src.algorithms.gnss.estimators.state_space import GnssStateSpace
from src.common_log import get_logger
from src.data_types.date import Epoch
from src.errors import PVTComputationFail
from src.io.config import config_dict
from src.io.config.enums import *
from src.models.observation.geometry import SystemGeometry
from src import constants

np.set_printoptions(linewidth=np.inf)


def get_weight(system_geometry, sat):
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

    def _set_solver_info(self, config):

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
        COMBINED_OBS_MODEL = {}
        MODEL = {}
        CODES = {}
        for const in ["GPS", "GAL"]:
            TROPO[const] = EnumTropo(config.get("model", const, "troposphere"))  # 0 - no model, 1 - Saastamoinen
            IONO[const] = EnumIono(config.get("model", const, "ionosphere"))  # 0 - no model, 1 - A priori

            COMBINED_OBS_MODEL[const] = EnumCombined(config.get("model", const, "combined"))  # false - uncombined,
            # true - combined

            # check if the model is single frequency or dual frequency
            code_types = self.obs_data.get_code_types(const)
            if len(code_types) > 1 and COMBINED_OBS_MODEL[const] == EnumCombined.UNCOMBINED_MODEL:
                # Dual-Frequency Model
                self.log.info(f"Dual-Frequency Model for constellation {const} and {COMBINED_OBS_MODEL[const]}")
                MODEL[const] = EnumModel.DUAL_FREQ
                CODES[const] = [code_types[0], code_types[1]]  # setting main and second code type
            else:
                # Single-Frequency Model
                self.log.info(f"Single-Frequency Model for constellation {const} and {COMBINED_OBS_MODEL[const]}")
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
            "COMBINED_OBS_MODEL": COMBINED_OBS_MODEL,
            "TX_TIME_ALG": EnumTransmissionTime(TX_TIME_ALG),
            "INITIAL_POS": INITIAL_POS,
            "INITIAL_CLOCK_BIAS": INITIAL_CLOCK_BIAS
        }

    def solve(self):
        """
        Output dict:
            receiver_pos (src.data_types.containers.TimeSeries.TimeSeries) : receiver position output timeseries store
            receiver_bias (src.data_types.containers.TimeSeries.TimeSeries) : receiver clock bias output timeseries
            estimated_iono (src.data_types.containers.TimeSeries.TimeSeries) : estimated ionosphere (for dual-frequency)
            prefit_residuals (src.data_types.containers.TimeSeries.TimeSeries) : prefit residuals output
            postfit_residuals (src.data_types.containers.TimeSeries.TimeSeries) : postfit residuals output
            DOPs (src.data_types.data_types.DOP.DOP) : DOPs output timeseries
            sat_info (src.data_types.containers.TimeSeries.TimeSeries) : satellite info time series
        """

        # available epochs
        epochs = self.obs_data.get_epochs()

        # initialize receiver_position
        initial_pos = self._info["INITIAL_POS"][0:3]
        initial_clock = self._info["INITIAL_CLOCK_BIAS"][0]
        state = GnssStateSpace(position=initial_pos,
                               clock_bias=initial_clock)

        # iterate over all available epochs
        for epoch in epochs:
            # fetch observation data for this epoch
            epoch_data = self.obs_data.get_epoch_data(epoch)

            # initialize solve-for variables (receiver position and bias) for the present epoch
            state = GnssStateSpace(position=state.position,
                                   clock_bias=state.clock_bias,
                                   date=epoch)
            # call lower level of solve
            _solver_info = self._solve(epoch, epoch_data, state)

            """ Re-add this later!!!
            if success:
                # add solution to Output timeseries
                self.log.info(f"Successfully solved positioning for epoch {epoch.to_time_stamp()} with "
                              f"RMS = {RMS} [m]")

                # store data for this epoch
                receiver_pos.set_data(epoch, state.receiver_position)
                receiver_bias.set_data(epoch, state.receiver_clock)
                if state.iono:
                    estimated_iono.set_data(epoch, state.iono)
                sat_info.set_data(epoch, _debug_info.get("geometry", None))
                DOPs.set_dop(epoch, "matrix", _debug_info.get("DOP", None))
                prefit_residuals.set_data(epoch, _debug_info.get("prefit", None))
                postfit_residuals.set_data(epoch, _debug_info.get("postfit", None))
            else:
                self.log.warning(f"PVT failed to converge for epoch {epoch.to_time_stamp()}. "
                                 f"No solution will be computed for this epoch.")"""

        self.log.info("Successfully ending module GNSS Positioning Solver...")

    @staticmethod
    def _stop(rms_old, rms_new, stop_criteria):
        return abs((rms_old - rms_new) / rms_old) <= stop_criteria

    def _solve(self, epoch, epoch_data, state):
        # begin iterative process for this epoch
        _solver_info = {}
        iteration = 0
        success = False
        RMS_prev = 1
        DOP = prefit_residuals = postfit_residuals = None
        gps_model = self._info['MODEL']['GPS']
        gps_obs_model = self._info['COMBINED_OBS_MODEL']['GPS']
        gps_codes = self._info['CODES']['GPS']

        # build system geometry for this epoch
        system_geometry = SystemGeometry(self.nav_data, epoch_data)

        # check which model to use (Single Frequency / Dual Frequency / not enough data)
        data_available = self._check_model_availability(system_geometry, epoch_data, epoch)

        if not data_available:
            return False, 0  # not enough data to process this epoch

        self.log.info(f"Processing {str(epoch)} with model GPS {gps_model} - {gps_obs_model}. Available Satellites: "
                      f"{system_geometry.get_satellites()}")

        # Iterated Least-Squares algorithm
        # Note: in the eventuality of implementing other filters, e.g.: Kalman Filter, then we need to branch
        # the solver from here...
        while iteration < self._info["MAX_ITER"]:

            # compute geometry-related data for each satellite link
            system_geometry.compute(epoch, state.position, state.clock_bias, self._info["TX_TIME_ALG"],
                                    gps_codes[0], self._info["REL_CORRECTION"])

            # solve the Least Squares
            try:
                if model == 0:
                    # Single Frequency Algorithm
                    postfit_residuals, DOP, prefit_residuals = self._solve_sf_LS(
                        system_geometry, epoch_data, state, epoch)
                else:
                    # Dual Frequency Algorithm
                    postfit_residuals, DOP, prefit_residuals = self._solve_df_LS(
                        system_geometry, epoch_data, state, epoch)

            except PVTComputationFail as e:
                self.log.warning(f"Least Squares failed for {str(epoch)} on iteration {iteration}."
                                 f"Reason: {e}")
                break

            # update RMS for this iteration
            RMS = np.linalg.norm(postfit_residuals)

            # check stop condition
            if self._stop(RMS_prev, RMS, self._info["STOP_CRITERIA"]):
                self.log.info(f"Least Squares was successful. Reached convergence at iteration {iteration}")
                success = True
                break

            # increase iteration counter
            RMS_prev = RMS
            iteration += 1

        # save debug_info
        _solver_info["geometry"] = system_geometry
        _solver_info["DOP"] = DOP
        _solver_info["prefit"] = prefit_residuals
        _solver_info["postfit"] = postfit_residuals
        _solver_info["success"] = success
        return _solver_info

    def _solve_sf_LS(self, system_geometry, epoch_data, state, nav_header, epoch):

        # initializations
        satellite_list = system_geometry.get_satellites()
        obs_length = len(satellite_list)  # number of available satellites
        state_length = len(state)  # number of solve for parameters to estimate

        # least squares arrays
        y = np.zeros(obs_length)  # observation vector <=> prefit residuals
        G = np.ones((obs_length, state_length))  # geometry matrix
        W = np.eye(obs_length)  # diagonal weight matrix

        observation_rec = ObservationReconstruction(system_geometry, self._info["MAIN_CODE"],
                                                    tropo=self._info["TROPO"] == 1,
                                                    iono=self._info["IONO"] == 1, true_range=True,
                                                    relativistic_correction=self._info["REL_CORRECTION"] == 1,
                                                    satellite_clock=True)

        i = 0
        for sat in satellite_list:
            # get observable
            obs = epoch_data.get_observable(sat, self._info["MAIN_CODE"])

            # fetch valid navigation message (closest to the current epoch)
            nav_message = self.nav_data.get_sat_data_for_epoch(sat, epoch)

            # compute predicted observation
            predicted_obs = observation_rec.compute(nav_message, nav_header, sat, epoch)

            # prefit residuals (measured observation - predicted observation)
            prefit_residuals = obs - predicted_obs

            # get LOS vector w.r.t. ECEF frame (column in geometry matrix)
            LOS = system_geometry.get_unit_line_of_sight(sat)

            # filling the corresponding entry of data vectors for the Least Squares quantities
            y[i] = prefit_residuals.value
            G[i][0:3] = LOS

            # Weight matrix -> sigma = 1 / e^{-elevation}
            W[i][i] = get_weight(system_geometry, sat)

            i += 1

        # solve LS problem for this iteration
        try:
            solver = WeightedLeastSquares(y, G, W=W)
            solver.solve()

            DOP = np.linalg.inv(G.T @ G)  # Dilution of precision matrix (without Weights)

        except (AttributeError, np.linalg.LinAlgError) as e:
            # possible error in the numpy.linalg.inv() function -> solution not possible
            raise PVTComputationFail(e)

        dX = solver.get_solution()

        # update state vector with incremental dX
        state.receiver_position.form = "cartesian"
        state.receiver_position += dX[0:3]
        state.receiver_clock = dX[3] / Constant.SPEED_OF_LIGHT  # receiver clock in seconds

        # get post-fit residuals
        post_fit = y - G[:, 0:3] @ dX[0:3]

        return post_fit, DOP, y

    def _solve_df_LS(self, system_geometry, epoch_data, state, nav_header, epoch):

        # initializations
        satellite_list = system_geometry.get_satellites()
        obs_length = len(satellite_list)  # number of available satellites

        # least squares arrays
        y = np.zeros(obs_length * 2)  # observation vector <=> prefit residuals
        G = np.ones((obs_length, 4))  # geometry matrix (state + clock)
        W = np.eye(obs_length * 2)  # diagonal weight matrix
        ionoMatrix1 = np.eye(obs_length) * (f1.freq_value / self._info["MAIN_CODE"].freq.freq_value) ** 2
        ionoMatrix2 = np.eye(obs_length) * (f1.freq_value / self._info["SECOND_CODE"].freq.freq_value) ** 2

        # fill in the LS arrays
        iFreq = 0
        for code in [self._info["MAIN_CODE"], self._info["SECOND_CODE"]]:
            observation_rec = ObservationReconstruction(system_geometry, code,
                                                        tropo=self._info["TROPO"] == 1,
                                                        iono=False, true_range=True,
                                                        relativistic_correction=self._info["REL_CORRECTION"] == 1,
                                                        satellite_clock=True)
            iSat = 0
            for sat in satellite_list:
                # get observable
                obs = epoch_data.get_observable(sat, code)

                # fetch valid navigation message (closest to the current epoch)
                nav_message = self.nav_data.get_sat_data_for_epoch(sat, epoch)

                # compute predicted observation
                predicted_obs = observation_rec.compute(nav_message, nav_header, sat, epoch)

                # prefit residuals (measured observation - predicted observation)
                prefit_residuals = obs - predicted_obs

                if iFreq == 0:
                    # get LOS vector w.r.t. ECEF frame (column in geometry matrix)
                    LOS = system_geometry.get_unit_line_of_sight(sat)
                    G[iSat][0:3] = LOS

                # filling the corresponding entry of data vectors for the Least Squares quantities
                y[iFreq * obs_length + iSat] = prefit_residuals.value

                # Weight matrix -> sigma = 1 / e^{-elevation}
                W[iFreq * obs_length + iSat][iFreq * obs_length + iSat] = get_weight(system_geometry, sat)

                iSat += 1
            iFreq += 1

        # concatenate G matrices and iono matrices
        system_matrix = np.block([[G, ionoMatrix1], [G, ionoMatrix2]])

        # solve LS problem for this iteration
        try:
            solver = WeightedLeastSquares(y, system_matrix, W=W)
            solver.solve()

            DOP = np.linalg.inv(G.T @ G)  # Dilution of precision matrix (without Weights)

        except (AttributeError, np.linalg.LinAlgError) as e:
            # possible error in the numpy.linalg.inv() function -> solution not possible
            raise PVTComputationFail(e)

        dX = solver.get_solution()

        # update state vector with incremental dX
        state.receiver_position.form = "cartesian"
        state.receiver_position += dX[0:3]
        state.receiver_clock = dX[3] / constants.SPEED_OF_LIGHT  # receiver clock in seconds
        state.iono = [(satellite_list[i], dX[i + 4]) for i in range(obs_length)]

        # get post-fit residuals
        post_fit = y[0:obs_length] - G[:, 0:3] @ dX[0:3]

        return post_fit, DOP, y

    def _check_model_availability(self, system_geometry, epoch_data, epoch):
        """
        TODO: this need to be revised
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
        model = self._info["MODEL"]["GPS"]
        codes = self._info["CODES"]["GPS"]
        MIN_SAT = 4

        sat_list = epoch_data.get_sats_for_datatypes(codes)
        if len(sat_list) < MIN_SAT:
            self.log.warning(
                f"Not enough satellites to compute the {model} PVT positioning at {str(epoch)}."
                f"Available satellites with data for codes {codes}: "
                f"{sat_list}. Minimum number of satellites is {MIN_SAT}. Aborting the solution for this epoch...")
            return False

        # Dual Frequency model
        if model == EnumModel.DUAL_FREQ:
            # removing satellites with data for only 1 frequency
            removed = []
            for sat in epoch_data.get_satellites():
                if sat not in sat_list:
                    system_geometry.remove(sat)
                    removed.append(sat)
            if removed:
                self.log.info(f"Removing satellites {removed} at epoch {str(epoch)} due to "
                              f"inconsistencies in code data for both frequencies")
        return True

    def _elevation_filter(self, system_geometry, iteration):
        # only apply the filter after iteration 3
        if iteration > 3:

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
                self.log.debug(f"Removing satellites {sats_to_remove} in iteration {iteration} due to "
                               f"elevation filter. ")
