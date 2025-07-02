import numpy as np

from src import constants
from src.data_types.gnss import DataType
from src.errors import SolverError
from src.modules.estimators.EKF import EKF
from src.modules.gnss.solver import PseudorangeReconstructor, CarrierPhaseReconstructor


def delete_state(x, P, index):
    """
    Delete the state at the given index from covariance matrix P.
    Removes both the corresponding row and column.
    """
    x_new = np.delete(x, index)
    P_new = np.delete(P, index, axis=0)  # Delete row
    P_new = np.delete(P_new, index, axis=1)  # Delete column
    return x_new, P_new


def add_state(x, P, index, new_x, new_var):
    """
    Add a new state value at the given index in state vector x and covariance matrix P.
    The new state is assumed uncorrelated with existing ones.

    Parameters:
        x       : state vector, shape (n,)
        P       : covariance matrix, shape (n, n)
        index   : position to insert the new state
        new_x   : value of the new state
        new_var : variance of the new state (scalar)

    Returns:
        x_out   : extended state vector, shape (n+1,)
        P_out   : extended covariance matrix, shape (n+1, n+1)
    """
    # Insert new value in state vector
    x_out = np.insert(x, index, new_x)

    # Create new covariance matrix with zeros
    n = P.shape[0]
    P_out = np.zeros((n + 1, n + 1))

    # Copy existing blocks
    P_out[:index, :index] = P[:index, :index]
    P_out[:index, index+1:] = P[:index, index:]
    P_out[index+1:, :index] = P[index:, :index]
    P_out[index+1:, index+1:] = P[index:, index:]

    # Set the variance of the new state
    P_out[index, index] = new_var

    return x_out, P_out

class EKF_Engine:
    def __init__(self, init_epoch, init_state, metadata):
        self._solver = EKF()
        # self._epoch = init_epoch
        self._state = init_state
        self._state.epoch = init_epoch
        self._init = False

        self._x = None
        self._P = None

        self.cp_based = metadata["CP_BASED"]
        self.pr_datatypes = metadata["CODES"]
        self.cp_datatypes = metadata["PHASES"]
        self.metadata = metadata



    @property
    def epoch(self):
        return self.state.epoch

    @property
    def state(self):
        return self._state

    def _build_init_state_cov(self, sat_list):
        # master_constellation = self._state.get_additional_info("clock_master")
        slave_constellation = self._state.get_additional_info("clock_slave")
        index_map = self._state.index_map
        n_states = index_map["total_states"]

        P0 = np.zeros((n_states, n_states))
        X0 = np.zeros(n_states)

        # position
        idx_pos = index_map["position"]
        for i in range(3):
            P0[idx_pos + i, idx_pos + i] = self._state.cov_position[idx_pos + i, idx_pos + i]
        X0[idx_pos:idx_pos+3] = self._state.position

        # clock bias (in meters)
        idx_clock = index_map["clock_bias"]
        P0[idx_clock, idx_clock] = self._state.cov_clock_bias
        X0[idx_clock] = self._state.clock_bias

        # tropo
        if "tropo_wet" in index_map:
            idx_tropo = index_map["tropo_wet"]
            P0[idx_tropo, idx_tropo] = self._state.cov_tropo_wet
            X0[idx_tropo] = self._state.tropo_wet

        # iono
        if "iono" in index_map:
            for sat in sat_list:
                if sat in index_map["iono"]:
                    idx_iono = index_map["iono"][sat]
                    cov_iono = self._state.cov_iono[sat]
                    P0[idx_iono, idx_iono] = cov_iono
                    X0[idx_iono] = self._state.iono[sat]

        # ISB
        if slave_constellation is not None and "isb" in index_map:
            idx_isb = index_map["isb"]
            P0[idx_isb, idx_isb] = self._state.cov_isb
            X0[idx_isb] = self._state.isb

        if "ambiguity" in index_map:
            pivot_dict = self._state.get_additional_info("pivot")
            for sat, cp_types in index_map["ambiguity"].items():
                if pivot_dict[sat.sat_system] != sat:
                    for cp_type in cp_types:
                        if not self._state.ambiguity[sat][cp_type].fixed:
                            idx_amb = cp_types[cp_type]
                            P0[idx_amb, idx_amb] = self._state.ambiguity[sat][cp_type].cov
                            X0[idx_amb] = self._state.ambiguity[sat][cp_type].val

        if "phase_bias" in index_map:
            for const, cp_types in index_map["phase_bias"].items():
                for cp_type in cp_types:
                    idx_phase_bias = cp_types[cp_type]
                    P0[idx_phase_bias, idx_phase_bias] = self._state.cov_phase_bias[const][cp_type]
                    X0[idx_phase_bias] = self._state.phase_bias[const][cp_type]

        return X0, P0

    def _build_state_cov(self, new_index_map, prev_index_map):
        x_out = self._x
        P_out = self._P

        if new_index_map["total_states"] != prev_index_map["total_states"]:
            # need to update the internal state and covariance dimensions (new states)

            # iono
            if "iono" in new_index_map:

                for sat in new_index_map["iono"]:
                    if sat not in prev_index_map["iono"]:
                        # add a new satellite (new available)
                        idx = new_index_map["iono"][sat]
                        x_out, P_out = add_state(x_out, P_out, idx, 0, 100)

                # update state and cov for removed satellites
                for sat in prev_index_map["iono"]:
                    if sat not in new_index_map["iono"]:
                        # remove this satellite (no longer available)
                        idx = prev_index_map["iono"][sat]
                        x_out, P_out = delete_state(x_out, P_out, idx)

        return x_out, P_out

    def _build_state_cov_2(self, sat_list, new_index_map, prev_index_map):
        # master_constellation = self._state.get_additional_info("clock_master")
        slave_constellation = self._state.get_additional_info("clock_slave")
        n_states = new_index_map["total_states"]

        X = np.zeros(n_states)
        P = self._state.get_additional_info("P")

        # position
        idx_pos = new_index_map["position"]
        X[idx_pos:idx_pos+3] = self._state.position

        # clock bias (in meters)
        idx_clock = new_index_map["clock_bias"]
        X[idx_clock] = self._state.clock_bias

        # tropo
        if "tropo_wet" in new_index_map:
            idx_tropo = new_index_map["tropo_wet"]
            X[idx_tropo] = self._state.tropo_wet

        # iono
        if "iono" in new_index_map:
            for sat in sat_list:
                if sat in new_index_map["iono"]:
                    idx_iono = new_index_map["iono"][sat]
                    X[idx_iono] = self._state.iono[sat]

        # ISB
        if slave_constellation is not None and "isb" in new_index_map:
            idx_isb = new_index_map["isb"]
            X[idx_isb] = self._state.isb

        if "ambiguity" in new_index_map:
            pivot_dict = self._state.get_additional_info("pivot")
            for sat, cp_types in new_index_map["ambiguity"].items():
                if pivot_dict[sat.sat_system] != sat:
                    for cp_type in cp_types:
                        if not self._state.ambiguity[sat][cp_type].fixed:
                            idx_amb = cp_types[cp_type]
                            X[idx_amb] = self._state.ambiguity[sat][cp_type].val

        if "phase_bias" in new_index_map:
            for const, cp_types in new_index_map["phase_bias"].items():
                for cp_type in cp_types:
                    idx_phase_bias = cp_types[cp_type]
                    X[idx_phase_bias] = self._state.phase_bias[const][cp_type]

        return X

    def _build_stm_process_noise(self, sat_list):
        # master_constellation = self._state.get_additional_info("clock_master")
        slave_constellation = self._state.get_additional_info("clock_slave")
        index_map = self._state.index_map
        n_states = index_map["total_states"]

        F = np.zeros((n_states, n_states))  # State Transition Matrix
        Q_c = np.zeros((n_states, n_states))  # Continuous-Time Process Noise Covariance Matrix

        # TODO: make configuration for random walk, white noise or gauss markov process

        process_pos = 100
        process_clock = 100
        process_iono = 100
        process_tropo = 100
        process_ambiguity = 100
        process_phase_bias = 100
        process_isb = 100

        # position
        idx_pos = index_map["position"]
        for i in range(3):
            Q_c[idx_pos + i, idx_pos + i] = process_pos
            F[idx_pos + i, idx_pos + i] = 1

        # clock bias (in meters)
        idx_clock = index_map["clock_bias"]
        Q_c[idx_clock, idx_clock] = process_clock
        F[idx_clock, idx_clock] = 1

        # tropo
        if "tropo_wet" in index_map:
            idx_tropo = index_map["tropo_wet"]
            Q_c[idx_tropo, idx_tropo] = process_tropo
            F[idx_tropo, idx_tropo] = 1

        # iono
        if "iono" in index_map:
            for sat in sat_list:
                if sat in index_map["iono"]:
                    idx_iono = index_map["iono"][sat]
                    Q_c[idx_iono, idx_iono] = process_iono
                    F[idx_iono, idx_iono] = 1

        # ISB
        if slave_constellation is not None and "isb" in index_map:
            idx_isb = index_map["isb"]
            Q_c[idx_isb, idx_isb] = process_isb
            F[idx_isb, idx_isb] = 1

        if "ambiguity" in index_map:
            pivot_dict = self._state.get_additional_info("pivot")
            for sat, cp_types in index_map["ambiguity"].items():
                if pivot_dict[sat.sat_system] != sat:
                    for cp_type in cp_types:
                        if not self._state.ambiguity[sat][cp_type].fixed:
                            idx_amb = cp_types[cp_type]
                            Q_c[idx_amb, idx_amb] = process_ambiguity
                            F[idx_amb, idx_amb] = 1

        if "phase_bias" in index_map:
            for const, cp_types in index_map["phase_bias"].items():
                for cp_type in cp_types:
                    idx_phase_bias = cp_types[cp_type]
                    Q_c[idx_phase_bias, idx_phase_bias] = process_phase_bias
                    F[idx_phase_bias, idx_phase_bias] = 1
        return F, Q_c

    def compute_residual_los(self, sat, epoch, datatype, obs_data, reconstructor):
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

        predicted_obs = reconstructor.compute(sat, epoch, datatype)

        # prefit residuals (observed minus computed)
        prefit_residuals = obs - predicted_obs

        # get LOS vector w.r.t. ECEF frame (column in geometry matrix)
        line_sight = reconstructor.get_unit_line_of_sight(sat)

        return prefit_residuals, line_sight

    def _build_obs_matrix(self, epoch, obs_data, datatypes, state, reconstructor_dict, sat_list):
        obs_offset = 0
        index_map = self._state.index_map
        slave_constellation = self._state.get_additional_info("clock_slave")
        n_states = index_map["total_states"]
        n_obs = obs_data.get_number_of_observations(datatypes)

        y_vec = np.zeros(n_obs)
        design_mat = np.zeros((n_obs, n_states))
        R = np.eye(n_obs)

        for const in datatypes.keys():
            for datatype in datatypes[const]:

                if DataType.is_code(datatype):
                    reconstructor = reconstructor_dict["PR"]
                elif DataType.is_carrier(datatype):
                    reconstructor = reconstructor_dict["CP"]
                else:
                    raise SolverError(f"Unknown datatype {datatype} at epoch {epoch}.")

                iSat = 0
                for sat in sat_list:
                    if sat.sat_system != const:
                        continue

                    residual, los = self.compute_residual_los(sat, epoch, datatype, obs_data, reconstructor)

                    # filling the LS matrices
                    y_vec[obs_offset + iSat] = residual

                    # position
                    idx_pos = index_map["position"]
                    design_mat[obs_offset + iSat][idx_pos:idx_pos+3] = los

                    # clock
                    idx_clock = index_map["clock_bias"]
                    design_mat[obs_offset + iSat, idx_clock] = 1.0

                    # tropo
                    if "tropo_wet" in index_map:
                        map_wet = reconstructor._system_geometry.get("tropo_map_wet", sat)
                        idx_tropo = index_map["tropo_wet"]
                        design_mat[obs_offset + iSat, idx_tropo] = map_wet

                    # iono
                    if "iono" in index_map and sat in index_map["iono"]:
                        factor = (datatypes[const][0].freq.freq_value / datatype.freq.freq_value) ** 2
                        idx_iono = index_map["iono"][sat]
                        if DataType.is_carrier(datatype):
                            design_mat[obs_offset + iSat, idx_iono] = -1.0 * factor  # iono
                        else:
                            design_mat[obs_offset + iSat, idx_iono] = 1.0 * factor  # iono

                    # ISB

                    if "isb" in index_map and sat.sat_system == slave_constellation:
                        idx_isb = index_map["isb"]
                        design_mat[obs_offset + iSat, idx_isb] = 1.0

                    if DataType.is_carrier(datatype):
                        # ambiguity
                        if "ambiguity" in index_map and sat in index_map["ambiguity"]:
                            pivot_dict = state.get_additional_info("pivot")
                            if pivot_dict[sat.sat_system] != sat:
                                if not state.ambiguity[sat][datatype].fixed:
                                    idx_amb = index_map["ambiguity"][sat][datatype]
                                    wavelength = constants.SPEED_OF_LIGHT / datatype.freq.freq_value
                                    design_mat[obs_offset + iSat, idx_amb] = wavelength

                        # phase bias
                        if "phase_bias" in index_map and const in index_map["phase_bias"]:
                            idx_phase_bias = index_map["phase_bias"][const][datatype]
                            design_mat[obs_offset + iSat, idx_phase_bias] = 1.0

                    # Weight matrix -> as 1/(obs_std^2)
                    std = reconstructor.get_obs_std(sat, datatype)
                    obs_data.get_observable(sat, datatype).set_std(std)
                    R[obs_offset + iSat, obs_offset + iSat] = std ** 2

                    iSat += 1
                obs_offset += iSat  # TODO check + or - here
        return R, design_mat, y_vec

    def _update_state(self, x_out, P_out, sat_list):

        index_map = self._state.index_map

        # Perform Ambiguity Resolution (if enabled)
        #if "ambiguity" in index_map and self._state.ambiguity.amb_resolution_enable:
        #    dX, cov = self._state.ambiguity.main_fix(index_map, self._state, dX, cov)
        idx_pos = index_map["position"]
        idx_clock = index_map["clock_bias"]

        self._state.position = x_out[idx_pos:idx_pos+3]
        self._state.clock_bias = x_out[idx_clock]

        # if iono is estimated
        if "iono" in index_map:
            for sat in sat_list:
                if sat in index_map["iono"]:
                    idx_iono = index_map["iono"][sat]
                    self._state.iono[sat] = x_out[idx_iono]
                    self._state.cov_iono[sat] = P_out[idx_iono, idx_iono]

        # ISB
        if "isb" in index_map:
            idx_isb = index_map["isb"]
            self._state.isb = x_out[idx_isb]
            self._state.cov_isb = P_out[idx_isb, idx_isb]

        # tropo
        if "tropo_wet" in index_map:
            idx_tropo = index_map["tropo_wet"]
            self._state.tropo_wet = x_out[idx_tropo]
            self._state.cov_tropo_wet = P_out[idx_tropo, idx_tropo]

        #if "ambiguity" in index_map:
        #    pivot_dict = self._state.get_additional_info("pivot")
        #    for sat, cp_types in index_map["ambiguity"].items():
        #        if pivot_dict[sat.sat_system] != sat:
        #            for cp_type in cp_types:
        #                idx_amb = cp_types[cp_type]
        #                if not self._state.ambiguity[sat][cp_type].fixed:
        #                    self._state.ambiguity[sat][cp_type].val += dX[idx_amb]
        #                    self._state.ambiguity[sat][cp_type].cov = cov[idx_amb, idx_amb]

        #if "phase_bias" in index_map:
        #    for const, cp_types in index_map["phase_bias"].items():
        #        for cp_type in cp_types:
        #            idx_phase_bias = cp_types[cp_type]
        #            self._state.phase_bias[const][cp_type] += dX[idx_phase_bias]
        #            self._state.cov_phase_bias[const][cp_type] = cov[idx_phase_bias, idx_phase_bias]

        # unpack covariance matrices
        self._state.cov_position = np.array(P_out[idx_pos:idx_pos+3, idx_pos:idx_pos+3])
        self._state.cov_clock_bias = P_out[idx_clock, idx_clock]

    def estimate(self, epoch, system_geometry, obs_for_epoch):

        reconstructor = dict()
        reconstructor["PR"] = PseudorangeReconstructor(system_geometry, self.metadata, self._state, None)

        if self.cp_based:
            datatypes = dict()
            for const in self.pr_datatypes.keys():
                datatypes[const] = dict()
                datatypes[const] = self.pr_datatypes[const] + self.cp_datatypes[const]
            reconstructor["CP"] = CarrierPhaseReconstructor(system_geometry, self.metadata, self._state, None)
        else:
            datatypes = self.pr_datatypes

        sat_list = system_geometry.get_satellites()

        if self._init is False:
            self._state.build_index_map(sat_list)
            x_in, P_in = self._build_init_state_cov(sat_list)
            self._init = True
            #self._x = x_in
            #self._P = P_in
        else:
            prev_index_map = self._state.index_map
            self._state.build_index_map(sat_list)
            new_index_map = self._state.index_map
            x_in, P_in = self._build_state_cov(new_index_map, prev_index_map)
            # P_in = self._state.get_additional_info("P")

        # build state vector and covariance matrix (the dimension may be variable)
        time_step = (epoch - self.epoch).total_seconds()
        F, Q_c = self._build_stm_process_noise(sat_list)
        R, H, obs_vector = self._build_obs_matrix(epoch, obs_for_epoch, datatypes, self._state, reconstructor, sat_list)

        # perform predict step from `self._epoch` to `epoch`
        x_pred, P_pred = self._solver.predict(x_in, P_in, time_step, F, Q_c)

        # perform update step
        x_out, P_out = self._solver.update(obs_vector, P_pred, x_pred, H, R)

        # update the state vector object with x^hat
        self._update_state(x_out, P_out, sat_list)

        self._x = x_out
        self._P = P_out
        # self._state.add_additional_info("P", P_out)
        self._state.epoch = epoch

        # build the postfit residuals vector
        #post_fit = self.y_vec - self.design_mat @ x_hat
        #norm = np.linalg.norm(post_fit)

        # build prefit and postfit residual dicts for output
        #pre_fit_dict = self.get_residuals(self.y_vec)
        #post_fit_dict = self.get_residuals(post_fit)
