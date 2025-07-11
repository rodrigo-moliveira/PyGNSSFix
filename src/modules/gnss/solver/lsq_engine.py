""" Module with Least-Squares engines for position and velocity estimators """
import numpy as np

from src import constants
from src.constants import SPEED_OF_LIGHT
from src.errors import SolverError
from src.io.config import EnumSolver
from src.modules.estimators.weighted_ls import WeightedLeastSquares
from src.modules.gnss.solver import PseudorangeReconstructor, RangeRateReconstructor, CarrierPhaseReconstructor
from src.data_types.gnss import DataType


class LSQ_Engine:
    """
    Base class for the LSQ Engine.

    Two child classes are derived:
        * :py:class:`LSQ_Engine_Position` for the estimation of position and clock bias (with pseudorange measurements)
        * :py:class:`LSQ_Engine_Velocity` for the estimation of velocity and clock drift (with range rate measurements)

    Consider the linear system y = G @ x, where
        * y - observation vector
        * x - state vector
        * G - state matrix (relates state vector with the observation vector)

    Given the vector of observations y, the unconstrained estimated state x^{hat} is given by (Normal Equation):
        x^{hat} = (G.T @ W @ G)^-1 @ G.T @ W @ y
    Given the vector of observations y, the initial state vector x0, the initial covariance matrix P0_inv,
    the constrained estimated state x^{hat} is given by:
        x^{hat} = (P0_inv + G.T @ W @ G)^-1 @ (G.T @ W @ y + P0_inv @ x0)
    where W is the weight matrix.

    The `LSQ_Engine` class is responsible for filling the LSQ observation vector y and state matrices G and W.
    Then it calls an :py:class:`src.modules.estimators.WeightedLeastSquares` instance to solve the system and update
    the GNSS state vector object (:py:class:`src.data_mng.gnss.state_space.GnssStateSpace`) accordingly. This process
    is performed iteratively (Iterated Weighted Least-Squares).
    
            
    Attributes:
        datatypes(dict): dict with all the datatypes to be used
        y_vec(numpy.ndarray): observation vector
        design_mat(numpy.ndarray): state matrix
        weight_mat(numpy.ndarray): weight matrix
        sat_list(dict): dict with available constellations as keys and available satellites to be used as values
    """

    def __init__(self, datatypes, satellite_list, epoch, obs_data, state, solver_enum):
        """
        Constructor of the LSQ_Engine instance.

        Args:
            datatypes(dict): dict with all the datatypes to be used (constellation as keys, list of datatypes as values)
            satellite_list(list): list with satellites to be processed
            epoch(src.data_types.date.Epoch): epoch to be solved
            obs_data (src.data_mng.gnss.observation_data.EpochData) : instance of `EpochData` (GNSS observable database
                for a single epoch)
            state(src.data_mng.gnss.state_space.GnssStateSpace): GNSS state vector object
            solver_enum(src.io.config.enums.EnumSolver): enumeration for the type of solver (LS or WLS)
        """
        self.datatypes = datatypes

        self.y_vec = None       # observation vector
        self.design_mat = None  # design matrix
        self.weight_mat = None  # weight matrix

        self.sat_list = dict()
        for sat in satellite_list:
            if sat.sat_system not in self.sat_list:
                self.sat_list[sat.sat_system] = list()
            self.sat_list[sat.sat_system].append(sat)

        self._initialize_matrices(state, obs_data.get_number_of_observations(datatypes))
        self._build_lsq(epoch, obs_data, state, solver_enum)

    def _initialize_matrices(self, state, n_obs):
        """
        Initializes the observation vector y (attribute `y_vec`), the state matrix G (attribute `design_mat`) 
        and the weight matrix W (attribute `weight_mat`) to numpy objects of the correct size and shape 
        (initialized to zeros), according to the user configurations

        The size of the arrays depend on the number of observations (satellites and frequencies) available, the number
        of states to be estimated, etc.

        Args:
            state(src.data_mng.gnss.state_space.GnssStateSpace): GNSS state vector object
            n_obs(int): number of observations to be processed (sum per satellite and datatype)
        """
        pass

    def _build_lsq(self, epoch, obs_data, state, solver_enum):
        """
        Fills the observation vector y (attribute `y_vec`), the state matrix G (attribute `design_mat`) and the weight
        matrix W (attribute `weight_mat`) with the appropriate data.

        Args:
            epoch(src.data_types.date.Epoch): epoch to be solved
            obs_data (src.data_mng.gnss.observation_data.EpochData) : instance of `EpochData` (GNSS observable database
                for a single epoch)
            state(src.data_mng.gnss.state_space.GnssStateSpace): GNSS state vector object
            solver_enum(src.io.config.enums.EnumSolver): enumeration for the type of solver (LS or WLS)
        """
        pass

    def _build_init_state_cov(self, state):
        """ Builds the initial state vector and covariance matrix for the constrained WLS problem.

        Args:
            state(src.data_mng.gnss.state_space.GnssStateSpace): GNSS state vector object
        Returns
            tuple[numpy.ndarray, numpy.ndarray]: Tuple with the initial state vector and the inverse of the initial
                                                covariance matrix
        """
        return None, None

    def compute_residual_los(self, sat, epoch, datatype, obs_data):
        """
        Computes the prefit residuals (observed minus computed observation) and the line of sight vector w.r.t. 
            ECEF frame.

        Args:
            sat(src.data_types.gnss.Satellite): satellite to compute the residual and LOS
            epoch(src.data_types.date.Epoch): epoch to make the computation
            datatype(src.data_types.gnss.DataType): datatype (frequency band) to compute the residual
            obs_data (src.data_mng.gnss.observation_data.EpochData) : instance of `EpochData` (GNSS observable database
                for a single epoch)

        Returns:
             tuple[float, numpy.ndarray] : Tuple with the computed residual (in the units of the datatype)
                                        and the line of sight for [x, y, z] axis of ECEF frame
        """
        pass

    def solve_ls(self, state, constrained=False):
        """
        Solves the LS problem for this iteration.

        Args:
            state(src.data_mng.gnss.state_space.GnssStateSpace): GNSS state vector object
            constrained(bool): if True, the LS problem is solved with constraints (using the initial state vector and
                covariance matrix). If False, the LS problem is solved without constraints (using only the
                observation vector and the design matrix)

        Returns:
            tuple[dict, dict, numpy.ndarray, float]:
                * first element is a dict with the prefit residuals (indexed by constellation, satellite and datatype)
                * second element is a dict with the postfit residuals (indexed by constellation, satellite and datatype)
                * third element is a numpy array object with the DOP matrix, DOP = (G.T * G)^-1
                * fourth element is the norm of the postfit residuals vector
        
        Raises:
            SolverError: an exception is raised if the LS problem fails to be solved
        """
        try:

            if constrained:
                X0, P0_inv = self._build_init_state_cov(state)
            else:
                X0, P0_inv = None, None
            solver = WeightedLeastSquares(self.y_vec, self.design_mat, W=self.weight_mat, P0_inv=P0_inv,
                                          X0=X0)
            solver.solve()
            dop_matrix = np.linalg.inv(self.design_mat.T @ self.weight_mat @ self.design_mat)

        except (AttributeError, np.linalg.LinAlgError) as e:
            # failed to solve the LS: solution not possible
            raise SolverError(e)

        # get estimated state x^hat and covariance
        x_hat = solver.get_solution()
        cov = solver.get_cov()

        # update the state vector object with x^hat
        self._update_state(state, x_hat, cov)

        # build the postfit residuals vector
        post_fit = self.y_vec - self.design_mat @ x_hat
        norm = np.linalg.norm(post_fit)

        # build prefit and postfit residual dicts for output
        pre_fit_dict = self.get_residuals(self.y_vec)
        post_fit_dict = self.get_residuals(post_fit)

        return pre_fit_dict, post_fit_dict, dop_matrix, norm

    def _update_state(self, state, dX, cov):
        """
        Applies computed corrections and updates the GNSS state vector accordingly.

        Args:
            state(src.data_mng.gnss.state_space.GnssStateSpace): GNSS state vector object
            dX(numpy.ndarray): solution array `x^hat` of the LSQ Solver
            cov(numpy.ndarray): covariance matrix `S = (G.T @ W @ G)^-1` of the LSQ Solver
        """
        pass

    def get_residuals(self, residual_vec):
        """
        Rearrange the residual vector (from a numpy.ndarray vector) into a more readable format.

        The output is a dictionary keyed by constellation, satellite, and datatype, such that:
            res_dict[constellation][sat][datatype] = residual

        Args:
            residual_vec (numpy.ndarray): The input residual vector to be rearranged.

        Returns:
            dict: A dictionary with the rearranged residuals.
        """

        res_dict = dict()
        for const in self.sat_list.keys():
            n_sats = len(self.sat_list[const])
            res_dict[const] = dict()

            for iSat, sat in enumerate(self.sat_list[const]):
                res_dict[const][sat] = dict()

                for iFreq, datatype in enumerate(self.datatypes[const]):
                    res_dict[const][sat][datatype] = residual_vec[iFreq * n_sats + iSat]
        return res_dict


class LSQ_Engine_Position(LSQ_Engine):
    """
    LSQ_Engine_Position class inherits from :py:class:`LSQ_Engine`.

    The `LSQ_Engine_Position` class performs estimation of position and clock for a GNSS system using pseudorange
    and carrier phase (optional) observations from GPS and GAL constellations. The system can handle single or dual
    frequencies per constellation.

    Since the pseudorange and carrier phase equations are non-linear with respect to the user position, they have to be
    linearized.

    The linearized PR and CP observation equation adopted in the model are:

        ΔPR = -LOS * Δr + ΔT + c * Δdt_r + ΔI + c * ΔISB
        ΔCP = -LOS * Δr + ΔT + c * Δdt_r + c * Δδ_r - ΔI + c * ΔISB + λ * ΔN

    where:
        - ΔPR: is the prefit residual pseudorange observation (true minus computed PR observation)
        - ΔCP: is the prefit residual carrier phase observation (true minus computed CP observation)
        - LOS: Line of sight vector w.r.t. ECEF frame
        - Δr: Change in receiver position
        - ΔT: Change in tropospheric wet delay (when tropo estimation is enabled). ΔT = Δzwd * map_wet
        - c: Speed of light
        - Δdt_r: Change in receiver clock bias
        - Δδ_r: Change in receiver phase bias
        - ΔI: Change in ionospheric delay (when iono estimation is enabled). ΔI = mu * dI
        - ΔISB: Change in Inter System Bias (only enabled for the slave constellations)
        - λ: Wavelength of the carrier phase observation
        - ΔN: Change in ambiguity (only enabled for the carrier phase observation)

    In the model implemented, the ambiguity is estimated as floating point number (not fixed). The fixing process, if
    enabled, is implemented a-posteriori. The estimated state is then refined with the fixed ambiguity value.

    Note that the residual ΔPR and ΔCP are defined as the difference between the observed and computed observations.
        * The computed pseudorange is obtained from the reconstructed pseudorange observation.
            See the modelled effects in :py:class:`src.modules.gnss.solver.obs_reconstructor.PseudorangeReconstructor`.
        * The computed carrier phase is obtained from the reconstructed carrier phase observation.
            See the modelled effects in :py:class:`src.modules.gnss.solver.obs_reconstructor.CarrierPhaseReconstructor`.

    The terms in the right-hand side of the equation above form the state vector of the LSQ system. To evidence them,
    the equation can be re-written in matrix form:

        ΔPR = [-LOS map_wet c mu c] * [Δr Δzwd Δdt_r dI ΔISB]^T
        ΔCP = [-LOS map_wet c c -mu c λ] * [Δr Δzwd Δdt_r Δδ_r dI ΔISB ΔN]^T

    where [Δr Δzwd Δdt_r Δdt_r dI ΔISB ΔN] is the state vector (x^hat) of the LSQ problem.
    Depending on the user configuration, some states may be enabled or disabled (Δr and Δdt_r are mandatory).

    After the Normal Equation is solved, the GNSS state vector is updated accordingly

        * r -> r + Δr
        * dt_r -> dt_r + Δdt_r
        * δ_r -> δ_r + Δδ_r (only for CP)
        * I -> I + dI
        * T -> T + Δzwd (only zwd is estimated)
        * ISB -> ISB + ΔISB
        * N -> N + ΔN (only for CP)

    These new quantities are then used in the next iteration of the LSQ, for the new computation of the predicted
    observations.
    """

    def __init__(self, system_geometry, metadata, epoch, obs_data, state, init_state, trace_data):
        self.cp_based = metadata["CP_BASED"]
        pr_datatypes = metadata["CODES"]
        cp_datatypes = metadata["PHASES"]
        self._initial_state = init_state.clone()
        self._initial_state.build_index_map(system_geometry.get_satellites())  # this updates potential new sat states

        self.reconstructor = dict()
        self.reconstructor["PR"] = PseudorangeReconstructor(system_geometry, metadata, state, trace_data)

        if self.cp_based:
            datatypes = dict()
            for const in pr_datatypes.keys():
                datatypes[const] = dict()
                datatypes[const] = pr_datatypes[const] + cp_datatypes[const]
            self.reconstructor["CP"] = CarrierPhaseReconstructor(system_geometry, metadata, state, trace_data)
        else:
            datatypes = pr_datatypes
        super().__init__(datatypes, system_geometry.get_satellites(), epoch, obs_data, state, metadata["SOLVER"])

    def _initialize_matrices(self, state, n_obs):
        """
        Initializes the observation vector y (attribute `y_vec`), the state matrix G (attribute `design_mat`)
        and the weight matrix W (attribute `weight_mat`) to numpy objects of the correct size and shape
        (initialized to zeros), according to the user configurations

        The size of the arrays depend on the number of observations (satellites and datatypes) available, the number
        of states to be estimated, etc.

        Example for Pseudorange only estimation:
            Single Frequency Single Constellation
                -> States are position and clock
                -> Design matrix has shape [m_S,4]
                -> Observation vector has shape [m_S,1]
            Dual Frequency Single Constellation
                -> States are position, clock and iono
                -> Design matrix has shape [2*m_S,4+m_S]
                -> Observation vector has shape [2*m_S,1]
            Single Frequency Dual Constellation
                -> States are position, clock and isb
                -> Design matrix has shape [m_S1+m_S2,4+1]
                -> Observation vector has shape [m_S1+m_S2,1]
            Dual Frequency Dual Constellation
                -> States are position, clock, iono and isb
                -> Design matrix has shape [2*(m_S1+m_S2), 4+m_S1+m_S2+1]
                -> Observation vector has shape [2*(m_S1+m_S2),1]
            Dual constellation (constellation 1 with dual frequency and constellation 2 with single frequency)
                -> States: position, clock, iono and isb
                -> Design matrix is Design matrix has shape [2*m_S1+m_S2, 4+m_S1+1]
                -> Observation vector has shape [2*m_S1+m_S2,1]
        m_Si is the number of available satellites for the constellation i
        When tropo estimation is enabled, a new state is added to all cases above.

        See ref [1] for more details about the different cases and setups of the LS Normal Equations.

        Reference:
            [1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
                Springer Cham, 2017
        """
        n_states = state.index_map["total_states"]

        self.y_vec = np.zeros(n_obs)
        self.design_mat = np.zeros((n_obs, n_states))
        self.weight_mat = np.eye(n_obs)

    def compute_residual_los(self, sat, epoch, datatype, obs_data):
        """
        Computes the prefit residuals (observed minus computed observation) and the line of sight vector w.r.t.
        ECEF frame.

        Args:
            sat(src.data_types.gnss.Satellite): satellite to compute the residual and LOS
            epoch(src.data_types.date.Epoch): epoch to make the computation
            datatype(src.data_types.gnss.DataType): datatype (frequency band) to compute the residual
            obs_data (src.data_mng.gnss.observation_data.EpochData) : instance of `EpochData` (GNSS observable database
                for a single epoch)
        """
        # get observable and compute predicted observable
        obs = float(obs_data.get_observable(sat, datatype))
        if DataType.is_code(datatype):
            reconstructor = self.reconstructor["PR"]
        elif DataType.is_carrier(datatype):
            reconstructor = self.reconstructor["CP"]
        else:
            raise SolverError(f"Unknown datatype {datatype} for satellite {sat} at epoch {epoch}.")
        predicted_obs = reconstructor.compute(sat, epoch, datatype)

        # prefit residuals (observed minus computed)
        prefit_residuals = obs - predicted_obs

        # get LOS vector w.r.t. ECEF frame (column in geometry matrix)
        line_sight = reconstructor.get_unit_line_of_sight(sat)

        return prefit_residuals, line_sight

    def _build_init_state_cov(self, state):
        n_observables, n_states = np.shape(self.design_mat)
        index_map = state.index_map

        initial_state = self._initial_state
        P0 = np.zeros((n_states, n_states))
        X0 = np.zeros(n_states)
        X0_prev = np.zeros(n_states)

        # position
        idx_pos = index_map["position"]
        for i in range(3):
            P0[idx_pos + i, idx_pos + i] = initial_state.cov_position[idx_pos + i, idx_pos + i]
        X0[idx_pos:idx_pos+3] = initial_state.position
        X0_prev[idx_pos:idx_pos+3] = state.position

        # clock bias (in meters)
        idx_clock = index_map["clock_bias"]
        P0[idx_clock, idx_clock] = initial_state.cov_clock_bias
        X0[idx_clock] = initial_state.clock_bias
        X0_prev[idx_clock] = state.clock_bias

        # tropo
        if "tropo_wet" in index_map:
            idx_tropo = index_map["tropo_wet"]
            P0[idx_tropo, idx_tropo] = initial_state.cov_tropo_wet
            X0[idx_tropo] = initial_state.tropo_wet
            X0_prev[idx_tropo] = state.tropo_wet

        for iConst, const in enumerate(self.sat_list.keys()):

            # iono
            for iSat, sat in enumerate(self.sat_list[const]):
                if "iono" in index_map and sat in index_map["iono"]:
                    idx_iono = index_map["iono"][sat]
                    cov_iono = initial_state.cov_iono[sat]
                    P0[idx_iono, idx_iono] = cov_iono
                    X0[idx_iono] = initial_state.iono[sat]
                    X0_prev[idx_iono] = state.iono[sat]

            # ISB
            if iConst > 0 and "isb" in index_map:
                idx_isb = index_map["isb"]
                P0[idx_isb, idx_isb] = initial_state.cov_isb
                X0[idx_isb] = initial_state.isb
                X0_prev[idx_isb] = state.isb

        if "ambiguity" in index_map:
            pivot_dict = state.get_additional_info("pivot")
            for sat, cp_types in index_map["ambiguity"].items():
                if pivot_dict[sat.sat_system] != sat:
                    for cp_type in cp_types:
                        if not initial_state.ambiguity[sat][cp_type].fixed:
                            idx_amb = cp_types[cp_type]
                            P0[idx_amb, idx_amb] = initial_state.ambiguity[sat][cp_type].cov
                            X0[idx_amb] = initial_state.ambiguity[sat][cp_type].val
                            X0_prev[idx_amb] = state.ambiguity[sat][cp_type].val

        if "phase_bias" in index_map:
            for const, cp_types in index_map["phase_bias"].items():
                for cp_type in cp_types:
                    idx_phase_bias = cp_types[cp_type]
                    P0[idx_phase_bias, idx_phase_bias] = initial_state.cov_phase_bias[const][cp_type]
                    X0[idx_phase_bias] = initial_state.phase_bias[const][cp_type]
                    X0_prev[idx_phase_bias] = state.phase_bias[const][cp_type]

        return X0 - X0_prev, np.linalg.inv(P0)

    def _build_lsq(self, epoch, obs_data, state, solver_enum):
        obs_offset = 0
        index_map = state.index_map

        for iConst, const in enumerate(self.sat_list.keys()):

            n_sats = len(self.sat_list[const])
            for iFreq, datatype in enumerate(self.datatypes[const]):

                if DataType.is_code(datatype):
                    reconstructor = self.reconstructor["PR"]
                elif DataType.is_carrier(datatype):
                    reconstructor = self.reconstructor["CP"]
                else:
                    raise SolverError(f"Unknown datatype {datatype} at epoch {epoch}.")

                for iSat, sat in enumerate(self.sat_list[const]):
                    residual, los = self.compute_residual_los(sat, epoch, datatype, obs_data)

                    # filling the LS matrices
                    self.y_vec[obs_offset + iSat] = residual

                    # position
                    idx_pos = index_map["position"]
                    self.design_mat[obs_offset + iSat][idx_pos:idx_pos+3] = los

                    # clock
                    idx_clock = index_map["clock_bias"]
                    self.design_mat[obs_offset + iSat, idx_clock] = 1.0

                    # tropo
                    if "tropo_wet" in index_map:
                        map_wet = reconstructor._system_geometry.get("tropo_map_wet", sat)
                        idx_tropo = index_map["tropo_wet"]
                        self.design_mat[obs_offset + iSat, idx_tropo] = map_wet

                    # iono
                    if "iono" in index_map and sat in index_map["iono"]:
                        factor = (self.datatypes[const][0].freq.freq_value / datatype.freq.freq_value) ** 2
                        idx_iono = index_map["iono"][sat]
                        if DataType.is_carrier(datatype):
                            self.design_mat[obs_offset + iSat, idx_iono] = -1.0 * factor  # iono
                        else:
                            self.design_mat[obs_offset + iSat, idx_iono] = 1.0 * factor  # iono

                    # ISB
                    if iConst > 0 and "isb" in index_map:
                        idx_isb = index_map["isb"]
                        self.design_mat[obs_offset + iSat, idx_isb] = 1.0

                    if DataType.is_carrier(datatype):
                        # ambiguity
                        if "ambiguity" in index_map and sat in index_map["ambiguity"]:
                            pivot_dict = state.get_additional_info("pivot")
                            if pivot_dict[sat.sat_system] != sat:
                                if not state.ambiguity[sat][datatype].fixed:
                                    idx_amb = index_map["ambiguity"][sat][datatype]
                                    wavelength = constants.SPEED_OF_LIGHT / datatype.freq.freq_value
                                    self.design_mat[obs_offset + iSat, idx_amb] = wavelength

                        # phase bias
                        if "phase_bias" in index_map and const in index_map["phase_bias"]:
                            idx_phase_bias = index_map["phase_bias"][const][datatype]
                            self.design_mat[obs_offset + iSat, idx_phase_bias] = 1.0

                    # Weight matrix -> as 1/(obs_std^2)
                    if solver_enum == EnumSolver.WLS:
                        std = reconstructor.get_obs_std(sat, datatype)
                        obs_data.get_observable(sat, datatype).set_std(std)
                        self.weight_mat[obs_offset + iSat, obs_offset + iSat] = \
                            1 / (std ** 2)

                obs_offset += n_sats

    def _update_state(self, state, dX, cov):

        index_map = state.index_map

        # Perform Ambiguity Resolution (if enabled)
        if "ambiguity" in index_map and state.ambiguity.amb_resolution_enable:
            dX, cov = state.ambiguity.main_fix(index_map, state, dX, cov)
        idx_pos = index_map["position"]
        idx_clock = index_map["clock_bias"]

        state.position = state.position + dX[idx_pos:idx_pos+3]
        state.clock_bias = state.clock_bias + dX[idx_clock]

        # if iono is estimated
        for const in self.sat_list.keys():
            if "iono" in index_map:
                for iSat, sat in enumerate(self.sat_list[const]):
                    if sat in index_map["iono"]:
                        idx_iono = index_map["iono"][sat]
                        state.iono[sat] += float(dX[idx_iono])
                        state.cov_iono[sat] = cov[idx_iono, idx_iono]

        # ISB
        if "isb" in index_map:
            idx_isb = index_map["isb"]
            state.isb += dX[idx_isb]
            state.cov_isb = float(cov[idx_isb, idx_isb])

        # tropo
        if "tropo_wet" in index_map:
            idx_tropo = index_map["tropo_wet"]
            state.tropo_wet += dX[idx_tropo]
            state.cov_tropo_wet = cov[idx_tropo, idx_tropo]

        if "ambiguity" in index_map:
            pivot_dict = state.get_additional_info("pivot")
            for sat, cp_types in index_map["ambiguity"].items():
                if pivot_dict[sat.sat_system] != sat:
                    for cp_type in cp_types:
                        idx_amb = cp_types[cp_type]
                        if not state.ambiguity[sat][cp_type].fixed:
                            state.ambiguity[sat][cp_type].val += dX[idx_amb]
                            state.ambiguity[sat][cp_type].cov = cov[idx_amb, idx_amb]

        if "phase_bias" in index_map:
            for const, cp_types in index_map["phase_bias"].items():
                for cp_type in cp_types:
                    idx_phase_bias = cp_types[cp_type]
                    state.phase_bias[const][cp_type] += dX[idx_phase_bias]
                    state.cov_phase_bias[const][cp_type] = cov[idx_phase_bias, idx_phase_bias]

        # unpack covariance matrices
        state.cov_position = np.array(cov[idx_pos:idx_pos+3, idx_pos:idx_pos+3])
        state.cov_clock_bias = float(cov[idx_clock, idx_clock])


class LSQ_Engine_Velocity(LSQ_Engine):
    """
    LSQ_Engine_Velocity class inherits from :py:class`LSQ_Engine`.

    The `LSQ_Engine_Velocity` class performs estimation of velocity and clock drift for a GNSS system using doppler
    observations from GPS and GAL constellations. Currently, the system can only handle single frequency data
    per constellation.

    Since the pseudorange rate equation is linear with respect to the user position, therefore no linearization is
    required and this process resorts to a standard WLS without the need to perform iterations.

    The range rate model is:
        PR_Rate = (v^sat - v_rec).LOS + c(clock_rate_rec - clock_rate_sat - rel_clock_rate_sat)

    where:
        - v^sat: inertial satellite velocity vector written in ECEF frame components
        - v_rec: inertial receiver velocity vector written in ECEF frame components
        - LOS: Line of sight vector w.r.t. ECEF frame
        - clock_rate_rec: receiver clock drift
        - clock_rate_sat: satellite clock drift
        - rel_clock_rate_sat: satellite relativistic clock drift correction

    Rearranging the terms in the equation above, we evidence in the right-hand side only the terms to be estimated
        PR_Rate - v^sat.LOS + c * (clock_rate_sat + rel_clock_rate_sat) = -v_rec.LOS + c(clock_rate_rec)

    The state vector (x^hat) of the LSQ Normal Equation is thus [v_rec, clock_rate_rec].

    NOTE: the LSQ System estimates the receiver inertial velocity, so it is required to remove the Earth rotation
    in order to get the ECEF velocity, cf. Eq. (21.29) of [1].
        v_rec(ECEF) = v_rec(ECI) - OMEGA_EARTH x r_rec

    See :py:class:`src.modules.gnss.solver.obs_reconstructor.RangeRateReconstructor` for the reconstruction of the
    pseudorange rate observation.

    Reference
        [1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
            Springer Cham, 2017
    """

    def __init__(self, system_geometry, metadata, epoch, obs_data, state, trace_data):
        # currently we always use one observation per constellation (dual-frequency is not possible)
        datatypes = metadata["DOPPLER"]
        for key in datatypes:
            datatypes[key] = datatypes[key][:1]
        self.reconstructor = RangeRateReconstructor(system_geometry, metadata, state, trace_data)
        super().__init__(datatypes, system_geometry.get_satellites(), epoch, obs_data, state, metadata["SOLVER"])

    def _initialize_matrices(self, state, n_obs):
        """
        Initializes the observation vector y (attribute `y_vec`), the state matrix G (attribute `design_mat`)
        and the weight matrix W (attribute `weight_mat`) to numpy objects of the correct size and shape
        (initialized to zeros), according to the user configurations

        The size of the arrays depend on the number of observations (satellites and frequencies) available, the number
        of states to be estimated, etc.

        Different cases:
            Single Frequency Single Constellation
                -> States are velocity and clock drift
                -> Design matrix has shape [m_S,3+1]
                -> Observation vector has shape [m_S,1]
            Single Frequency Dual Constellation
                -> States are velocity and clock drift per constellation
                -> Design matrix has shape [m_S1+m_S2,3+1+1]
                -> Observation vector has shape [m_S1+m_S2,1]
        """
        n_observables = 0  # number of rows
        n_states = 3  # number of columns (default is 3 - velocity)

        for const in self.sat_list.keys():
            n_sats = len(self.sat_list[const])
            n_observables += n_sats
            n_states += 1  # receiver clock rate for each constellation

        self.y_vec = np.zeros(n_observables)
        self.design_mat = np.zeros((n_observables, n_states))
        self.weight_mat = np.eye(n_observables)

    def compute_residual_los(self, sat, epoch, doppler_datatype, obs_data):

        # get observable and compute predicted observable
        obs = float(obs_data.get_observable(sat, doppler_datatype))  # in Hz

        # transform Doppler to pseudorange rate
        wavelength = SPEED_OF_LIGHT / doppler_datatype.freq_value  # in meters
        obs_range_rate = -wavelength * obs  # in m/s

        predicted_obs = self.reconstructor.compute(sat, epoch, doppler_datatype)

        # prefit residuals (observed minus computed)
        prefit_residuals = obs_range_rate - predicted_obs

        los = -self.reconstructor.get_unit_line_of_sight(sat)

        return prefit_residuals, los

    def _build_lsq(self, epoch, obs_data, state, solver_enum):
        obs_offset = 0
        const_offset = 0
        for iConst, const in enumerate(self.sat_list.keys()):

            n_sats = len(self.sat_list[const])
            datatype = self.datatypes[const][0]

            for iSat, sat in enumerate(self.sat_list[const]):
                residual, los = self.compute_residual_los(sat, epoch, datatype, obs_data)

                # filling the LS matrices
                self.y_vec[obs_offset + iSat] = residual
                self.design_mat[obs_offset + iSat][0:3] = -los  # velocity
                self.design_mat[obs_offset + iSat][3 + const_offset] = 1

                if solver_enum == EnumSolver.WLS:
                    # Weight matrix -> as 1/(obs_std^2)
                    std = self.reconstructor.get_obs_std(sat, datatype)
                    obs_data.get_observable(sat, datatype).set_std(std)
                    self.weight_mat[obs_offset + iSat, obs_offset + iSat] = \
                        1 / (std ** 2)

            obs_offset += n_sats
            const_offset += 1

    def _update_state(self, state, dX, cov):
        vel = np.array(dX[0:3]) - np.cross(constants.EARTH_ANGULAR_RATE, state.position)
        state.velocity = vel  # in m/s
        state.cov_velocity = cov[0:3, 0:3]  # in (m/s)^2

        for iConst, const in enumerate(self.sat_list.keys()):
            # receiver clock drift [m/m]
            state.clock_bias_rate[const] = dX[iConst + 3]
            state.cov_clock_bias_rate[const] = float(cov[iConst + 3, iConst + 3])

    def _build_init_state_cov(self, state):
        n_observables, n_states = np.shape(self.design_mat)

        P0 = np.zeros((n_states, n_states))
        X0 = np.zeros(n_states)

        # position
        for i in range(3):
            P0[i, i] = state.cov_velocity[i, i]
        X0[0:3] = state.velocity + np.cross(constants.EARTH_ANGULAR_RATE, state.position)

        for iConst, const in enumerate(self.sat_list.keys()):
            # receiver clock drift [m/m]
            X0[3+iConst] = state.clock_bias_rate[const]
            P0[3 + iConst, 3 + iConst] = state.cov_clock_bias_rate[const]
        return X0, np.linalg.inv(P0)
