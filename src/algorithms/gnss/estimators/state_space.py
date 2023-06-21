from numpy.linalg import norm
import numpy as np

from PositioningSolver.src import get_logger
from PositioningSolver.src.algorithms.estimators.state_space import SPPStateSpace
from PositioningSolver.src.algorithms.estimators.weighted_ls import WeightedLeastSquares
from PositioningSolver.src.gnss.observation_models.geometry_obs import SystemGeometry
from PositioningSolver.src.gnss.observation_models.observation_reconstructor import ObservationReconstruction
from PositioningSolver.src.math_utils.Constants import Constant
from PositioningSolver.src.utils.errors import ConfigError, PVTComputationFail
from PositioningSolver.src.data_types.basics.DataType import DataType, DataTypeFactory
from PositioningSolver.src.gnss.observation_models import clock_obs

np.set_printoptions(linewidth=np.inf)


C1 = DataTypeFactory("C1")
f1 = C1.freq


class GPSSolver:
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


        Iterated Least Squares PVT - High Level Algorithm
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
    MODEL = {0: "Single Frequency",
             1: "Dual Frequency"}

    SOLVER = {0: "Least Squares",
              1: "Weighted Least Squares"}

    def __init__(self, obs_data, nav_data, config):
        """

        Args:
            obs_data (src.data_types.containers.ObservationData.ObservationData) : observation data
            nav_data (src.data_types.containers.NavigationData.NavigationDataMap) : navigation data
            config (src.config.Config) : user configurations
        """
        self.obs_data = obs_data
        self.nav_data = nav_data

        self.log = get_logger("gps_solver")
        self.log.info("#########################################################")
        self.log.info("######### Starting module 'GPS PVT Solver' ... ##########")

        # user configurations  #
        self._info, self.compute_TX_time = self._set_solver_info(config, obs_data.get_types())

    def _set_solver_info(self, config, datatypes):

        # Fetching user options
        MAX_ITER = config["gps_solver"]["iterations"]["select"]  # maximum number of iterations
        STOP_CRITERIA = config["gps_solver"]["stop_criteria"]["select"]  # RMS threshold for stop criteria
        MODEL = config["model"]["obs_combination"]["select"]  # 0 - SF, 1 - DF
        SOLVER = config["gps_solver"]["solution_solver"]["select"]  # 0 - LS, 1 - WLS
        TROPO = config["model"]["troposphere"]["select"]  # 0 - no model, 1 - Saastamoinen
        IONO = config["model"]["ionosphere"]["select"]  # 0 - no model, 1 - Klobuchar, 2 - Iono Free Combination
        REL_CORRECTION = config["model"]["relativistic_corrections"]["select"]  # 0 disable, 1 enable

        SIGNAL_STRENGTH_FILTER = config["gps_solver"]["signal_strength_filter"]["select"]  # threshold (in dBHz)
        ELEVATION_FILTER = config["gps_solver"]["elevation_filter"]["select"]
        SATELLITE_STATUS_FILTER = config["gps_solver"]["satellite_status"]

        # Checking Additional information

        # find number of necessary observations per epoch
        NR_EQS = 4  # base number of unidimensional equations/observations for each epoch

        # algorithm to compute transmission time
        if config["gps_solver"]["transmission_time_alg"]["select"] == 0:
            compute_TX_time = clock_obs.compute_TX_time_geometric
        else:
            compute_TX_time = clock_obs.compute_TX_time_pseudorange

        # fetch main and second data code types
        code_types = sorted(DataType.get_code_datatypes(datatypes))

        # dual frequency
        if MODEL == 1:
            if len(code_types) < 2:
                self.log.warning(f"User selected dual frequency PVT but GPS Solver only has one 1 code type available,"
                                 f" {code_types}. Resorting to Single Frequency PVT")
                MODEL = 0

            MAIN_CODE = code_types[0]
            SECOND_CODE = code_types[1]

        else:
            if len(code_types) == 0:
                raise ConfigError(f"No available code datatypes to perform PVT. Exiting.")
            MAIN_CODE = code_types[0]
            SECOND_CODE = None
        self.log.info(f"Main code for PVT: {MAIN_CODE}, second code: {SECOND_CODE}")
        self.log.info(f"SPP Algorithm - {GPSSolver.MODEL[MODEL]}. Solver - {GPSSolver.SOLVER[SOLVER]}")

        # add more info if necessary
        # ...

        # fill info dict
        _info = {
            "MAX_ITER": MAX_ITER,
            "STOP_CRITERIA": STOP_CRITERIA,
            "MODEL": MODEL,
            "SOLVER": SOLVER,
            "TROPO": TROPO,
            "IONO": IONO,
            "REL_CORRECTION": REL_CORRECTION,
            "SIGNAL_STRENGTH_FILTER": SIGNAL_STRENGTH_FILTER,
            "ELEVATION_FILTER": ELEVATION_FILTER,
            "SATELLITE_STATUS_FILTER": SATELLITE_STATUS_FILTER,
            "NR_EQS": NR_EQS,
            "MAIN_CODE": MAIN_CODE,
            "SECOND_CODE": SECOND_CODE
        }

        return _info, compute_TX_time

    def solve(self, receiver_pos, receiver_bias, prefit_residuals, estimated_iono,
              postfit_residuals, DOPs, sat_info):
        """

        Args:
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
        previous_state = SPPStateSpace()

        # iterate over all available epochs
        for epoch in epochs:

            # fetch observation data for this epoch
            epoch_data = self.obs_data.get_epoch_data(epoch)

            # initialize solve-for variables (receiver position and bias) for the present epoch
            state = SPPStateSpace(receiver_position=previous_state.receiver_position.copy(),
                                  receiver_clock=previous_state.receiver_clock)
            state.receiver_position.date = epoch

            # fetch closest navigation message header
            nav_header = self.nav_data.get_header_data(epoch)

            # call lower level of solve
            _debug_info = {}
            success, RMS = self._solve(epoch, epoch_data, state, nav_header, _debug_info)

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
                                 f"No solution will be computed for this epoch.")

            previous_state = state

        self.log.info("########## End of module 'GPS PVT Solver' ... ###########\n")

    @staticmethod
    def _stop(RMS_old, RMS_new, STOP_CRITERIA):
        return abs((RMS_old - RMS_new) / RMS_old) <= STOP_CRITERIA

    def _solve(self, epoch, epoch_data, state, nav_header, _debug_info):
        # begin iterative process for this epoch
        iteration = 0
        success = False
        RMS_prev = RMS = 1
        DOP = prefit_residuals = postfit_residuals = None

        # URA, Satellite health filters. Flagged satellites are removed
        self._initial_satellite_validation(epoch, epoch_data)

        # system geometry manager
        system_geometry = SystemGeometry(self.nav_data, nav_header, epoch_data)

        # check which model to use (Single Frequency / Dual Frequency / no model -> not enough data)
        control, model = self._check_model_availability(system_geometry, epoch_data, epoch)

        if not control:
            return False, 0  # not enough data to process this epoch

        # self.log.info(f"Processing {epoch.to_time_stamp()} with model {GPSSolver.MODEL[model]}. Available Satellites:"
        #              f"{system_geometry.get_satellites()}")

        # Iterated Least Squares algorithm (Note: in the eventuality of implementing other filters, e.g.: Kalman Filter,
        # then we need to branch the solver from here..)
        while iteration < self._info["MAX_ITER"]:

            # compute geometry-related data for each satellite link
            system_geometry.compute(epoch, state.receiver_position, state.receiver_clock, self.compute_TX_time,
                                    self._info["MAIN_CODE"], self._info["REL_CORRECTION"])

            # apply elevation filter
            self._elevation_filter(system_geometry, iteration)

            # solve the Least Squares
            try:
                if model == 0:
                    # Single Frequency Algorithm
                    postfit_residuals, DOP, prefit_residuals = self._solve_sf_LS(
                        system_geometry, epoch_data, state, nav_header, epoch)
                else:
                    # Dual Frequency Algorithm
                    postfit_residuals, DOP, prefit_residuals = self._solve_df_LS(
                        system_geometry, epoch_data, state, nav_header, epoch)

            except PVTComputationFail as e:
                self.log.warning(f"Least Squares failed for {epoch.to_time_stamp()} on iteration {iteration}\n"
                                 f"Reason: {e}")
                break

            # update RMS for this iteration
            RMS = np.linalg.norm(postfit_residuals)

            # check stop condition
            if self._stop(RMS_prev, RMS, self._info["STOP_CRITERIA"]):
                # self.log.info(f"Least Squares was successful. Reached convergence at iteration {iteration}")
                success = True
                break

            # increase iteration counter
            RMS_prev = RMS
            iteration += 1

        # save debug_info
        _debug_info["geometry"] = system_geometry
        _debug_info["DOP"] = DOP
        _debug_info["prefit"] = prefit_residuals
        _debug_info["postfit"] = postfit_residuals
        return success, RMS

    def get_weight(self, system_geometry, sat):
        sigma_elevation = np.e ** (-system_geometry.get("el", sat))
        w = (1 / sigma_elevation) ** 2

        return w

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
            W[i][i] = self.get_weight(system_geometry, sat)

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
        y = np.zeros(obs_length*2)  # observation vector <=> prefit residuals
        G = np.ones((obs_length, 4))  # geometry matrix (state + clock)
        W = np.eye(obs_length*2)  # diagonal weight matrix
        ionoMatrix1 = np.eye(obs_length) * (f1.freq_value / self._info["MAIN_CODE"].freq.freq_value)**2
        ionoMatrix2 = np.eye(obs_length) * (f1.freq_value / self._info["SECOND_CODE"].freq.freq_value)**2

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
                W[iFreq * obs_length + iSat][iFreq * obs_length + iSat] = self.get_weight(system_geometry, sat)

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
        state.receiver_clock = dX[3] / Constant.SPEED_OF_LIGHT  # receiver clock in seconds
        state.iono = [(satellite_list[i], dX[i+4]) for i in range(obs_length)]

        # get post-fit residuals
        post_fit = y[0:obs_length] - G[:, 0:3] @ dX[0:3]

        return post_fit, DOP, y

    def _check_model_availability(self, system_geometry, epoch_data, epoch):
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
        # user selected model (Dual Frequency or Single Frequency)
        model = self._info["MODEL"]
        MIN_SAT = self._info["NR_EQS"]

        # Dual Frequency model
        if model == 1:
            available_sat_list = \
                epoch_data.get_satellites_for_datatypes(self._info["MAIN_CODE"], self._info["SECOND_CODE"])

            if len(available_sat_list) < MIN_SAT:
                self.log.warning(
                    f"Not enough satellites to compute the Dual Frequency PVT positioning at {epoch.to_time_stamp()}."
                    f"Available satellites with data for codes {self._info['MAIN_CODE']} and "
                    f"{self._info['SECOND_CODE']} -- {available_sat_list}. Minimum number of satellites is {MIN_SAT}.")
                return False, model
            else:
                # removing satellites with data for only 1 frequency
                removed = []
                for sat in epoch_data.get_satellites():
                    if sat not in available_sat_list:
                        system_geometry.remove(sat)
                        removed.append(sat)
                if removed:
                    self.log.info(f"Removing satellites {removed} at epoch {epoch.to_time_stamp()} due to "
                                  f"inconsistencies in code data for both frequencies")

        elif model == 0:
            available_sat_list = epoch_data.get_satellites_for_datatypes(self._info["MAIN_CODE"])

            if len(available_sat_list) < MIN_SAT:
                self.log.warning(
                    f"Not enough satellites to compute the PVT positioning at {epoch.to_time_stamp()}.\n"
                    f"Reason: available satellites -- {available_sat_list}\n"
                    f"Minimum number of satellites is {MIN_SAT}. Aborting the solution for this epoch")
                return False, model

        return True, model

    def _elevation_filter(self, system_geometry, iteration):
        # only apply the filter after iteration 3
        if iteration > 3:

            # get elevation threshold in radians
            sats_to_remove = []
            threshold = self._info["ELEVATION_FILTER"]
            for sat, sat_info in system_geometry.items():

                if sat_info.el * Constant.RAD2DEG < threshold:
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

    def _initial_satellite_validation(self, epoch, epoch_data):
        """
        evaluates if this satellite can be used in PVT solver:
            * checks SV accuracy index (SV_URA). **REF[3]** advises to only use satellites for URA < 6144 [m]
            * checks SV_health. if SV_health != 0, then discard satellite
        """

        vSats = epoch_data.get_satellites()
        sats_to_remove = []
        for sat in vSats:

            # fetch valid navigation message (closest to the current epoch)
            nav_message = self.nav_data.get_sat_data_for_epoch(sat, epoch)

            # URA filter
            if self._info["SATELLITE_STATUS_FILTER"]["SV_URA"]:
                URA = nav_message.SV_URA
                URA_threshold = self._info["SATELLITE_STATUS_FILTER"]["SV_minimum_URA"]

                if URA > URA_threshold:
                    self.log.warning(f"Satellite {sat} is being discarded at epoch {epoch.to_time_stamp()} due to high "
                                     f"URA value ({URA}) compared to threshold {URA_threshold}")
                    sats_to_remove.append(sat)
                    continue

            # Satellite Health filter
            if self._info["SATELLITE_STATUS_FILTER"]["SV_health"]:
                SV_health = nav_message.SV_health

                if SV_health != 0:
                    self.log.warning(f"Satellite {sat} is being discarded at epoch {epoch.to_time_stamp()} due to bad "
                                     f"health flag. SV_health = {SV_health} in the navigation message")
                    sats_to_remove.append(sat)
                    continue

        # remove flagged satellites
        for sat in sats_to_remove:
            epoch_data.remove_satellite(sat)
