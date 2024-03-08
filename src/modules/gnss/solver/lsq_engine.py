import numpy as np

from src import constants
from src.modules.estimators import WeightedLeastSquares
from src.constants import SPEED_OF_LIGHT
from src.errors import SolverError
from src.io.config.enums import EnumModel


class LSQ_Engine:
    def __init__(self, satellite_list, metadata, epoch, obs_data, reconstructor, nav_data):
        self.y_vec = None  # observation vector
        self.design_mat = None  # design matrix
        self.weight_mat = None  # weight matrix

        self._metadata = metadata
        self.constellations = metadata["CONSTELLATIONS"]  # master constellation is the first in the list
        self.sat_list = dict()
        for sat in satellite_list:
            if sat.sat_system not in self.sat_list:
                self.sat_list[sat.sat_system] = list()
            self.sat_list[sat.sat_system].append(sat)

        self._initialize_matrices()
        self._build_lsq(epoch, obs_data, reconstructor, nav_data)

    def _initialize_matrices(self):

        # different cases:
        #   * Single Frequency Single Constellation: States: position, clock
        #       -> Design matrix is [m_S,4]
        #       -> Observation vector is [m_S,1]
        #   * Dual Frequency Single Constellation: States: position, clock, iono
        #       -> Design matrix is [2*m_S,4+m_S]
        #       -> Observation vector is [2*m_S,1]
        #   * Single Frequency Dual Constellation: States: position, clock, isb
        #       -> Design matrix is [m_S1+m_S2,4+1]
        #       -> Observation vector is [m_S1+m_S2,1]
        #   * Dual Frequency Dual Constellation: States: position, clock, iono, isb
        #       -> Design matrix is [2*(m_S1+m_S2), 4+m_S1+m_S2+1]
        #       -> Observation vector is [2*(m_S1+m_S2),1]
        #   * (Single + Dual) Frequency Dual Constellation:
        #       -> States: position, clock, iono, isb. Design matrix is Design matrix is [2*m_S1+m_S2, 4+m_S1+1]
        #       -> Observation vector is [2*m_S1+m_S2,1]
        #  When tropo is added, a new state is incremented, that is 4 -> 4 + 1 in the dimensions above

        n_observables = 0  # number of rows
        n_states = 3  # number of columns (default is 3 - position)

        for const in self.constellations:
            n_states += 1  # clock/isb variable

            n_sats = len(self.sat_list[const])

            # add iono states
            if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
                n_states += n_sats
                n_observables += 2 * n_sats
            else:
                n_observables += n_sats

        if self._metadata["TROPO"].estimate_tropo():
            n_states += 1

        self.y_vec = np.zeros(n_observables)
        self.design_mat = np.zeros((n_observables, n_states))
        self.weight_mat = np.eye(n_observables)

    @staticmethod
    def compute_residual_los(nav_data, sat, epoch, datatype, obs_data, reconstructor):
        # fetch valid navigation message (closest to the current epoch)
        nav_message = nav_data.get_closest_message(sat, epoch)

        # get observable and compute predicted observable
        obs = obs_data.get_observable(sat, datatype)
        predicted_obs = reconstructor.compute(nav_message, sat, epoch, datatype)

        # prefit residuals (measured gnss_models - predicted gnss_models)
        prefit_residuals = obs - predicted_obs

        # get LOS vector w.r.t. ECEF frame (column in geometry matrix)
        line_sight = reconstructor.get_unit_line_of_sight(sat)

        return prefit_residuals.value, line_sight

    def _build_lsq(self, epoch, obs_data, reconstructor, nav_data):
        """build the LS matrices y_vec, design_mat, weight_mat"""

        iono_offset = 0
        obs_offset = 0
        estimate_tropo = self._metadata["TROPO"].estimate_tropo()
        tropo_offset = 1 if estimate_tropo else 0

        for iConst, const in enumerate(self.constellations):

            n_sats = len(self.sat_list[const])
            for iFreq, datatype in enumerate(self._metadata["CODES"][const]):

                for iSat, sat in enumerate(self.sat_list[const]):
                    residual, los = self.compute_residual_los(nav_data, sat, epoch, datatype, obs_data, reconstructor)

                    # filling the LS matrices
                    self.y_vec[obs_offset + iSat] = residual
                    self.design_mat[obs_offset + iSat][0:3] = los  # position
                    self.design_mat[obs_offset + iSat, 3] = 1.0  # clock
                    if estimate_tropo:
                        map_wet = reconstructor._system_geometry.get("tropo_map_wet", sat)
                        self.design_mat[obs_offset + iSat, 4] = map_wet
                    if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
                        factor = (self._metadata["CODES"][const][0].freq.freq_value / datatype.freq.freq_value) ** 2
                        self.design_mat[obs_offset + iSat, 4 + tropo_offset + iono_offset + iSat] = 1.0 * factor  # iono
                    if iConst > 0:
                        self.design_mat[obs_offset + iSat, -1] = 1.0  # ISB

                    # Weight matrix -> as 1/(obs_std^2)
                    self.weight_mat[obs_offset + iSat, obs_offset + iSat] = \
                        1 / (reconstructor.get_obs_std(sat, datatype) ** 2)

                obs_offset += n_sats
            if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
                iono_offset += n_sats

    def solve_ls(self, state):
        """solves the LS problem for this iteration"""
        try:
            solver = WeightedLeastSquares(self.y_vec, self.design_mat, W=self.weight_mat)
            solver.solve()
            dop_matrix = np.linalg.inv(self.design_mat.T @ self.design_mat)

        except (AttributeError, np.linalg.LinAlgError) as e:

            # possible error in the numpy.linalg.inv() function -> solution not possible
            raise SolverError(e)

        # update state vector with incremental dX
        dX = solver.get_solution()
        cov = solver.get_cov()

        # apply dX corrections to the state vector
        self.apply_corrections(state, dX, cov)

        # form the post-fit residuals
        post_fit = self.y_vec - self.design_mat[:, 0:3] @ dX[0:3]
        norm = np.linalg.norm(post_fit)

        # form residual dicts for output
        pre_fit_dict = self.get_residuals(self.y_vec)
        post_fit_dict = self.get_residuals(post_fit)

        return pre_fit_dict, post_fit_dict, dop_matrix, norm

    def apply_corrections(self, state, dX, cov):
        """applies corrections to the state vector"""

        estimate_tropo = self._metadata["TROPO"].estimate_tropo()
        tropo_offset = 1 if estimate_tropo else 0

        state.position += dX[0:3]
        state.clock_bias += dX[3] / constants.SPEED_OF_LIGHT  # receiver clock in seconds

        # if iono is estimated
        iono_offset = 0
        for const in self.constellations:
            if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
                for iSat, sat in enumerate(self.sat_list[const]):
                    state.iono[sat] += float(dX[iono_offset + iSat + 4 + tropo_offset])
                    state.cov_iono[sat] = cov[
                        iono_offset + iSat + 4 + tropo_offset, iono_offset + iSat + 4 + tropo_offset]
                iono_offset += len(self.sat_list[const])

        # ISB
        if len(self.constellations) > 1:
            state.isb += dX[-1] / constants.SPEED_OF_LIGHT  # ISB between master and slave constellations
            state.cov_isb = float(cov[-1, -1]) / (constants.SPEED_OF_LIGHT ** 2)  # in seconds^2

        # tropo
        if estimate_tropo:
            state.tropo_wet += dX[4]
            state.cov_tropo_wet = cov[4, 4]
            print("tropo", state.tropo_wet)

        # unpack covariance matrices
        state.cov_position = np.array(cov[0:3, 0:3])
        state.cov_clock_bias = float(cov[3, 3]) / (constants.SPEED_OF_LIGHT ** 2)  # in seconds^2

    def get_residuals(self, residual_vec):

        res_dict = dict()
        for const in self.constellations:
            n_sats = len(self.sat_list[const])
            res_dict[const] = dict()

            for iSat, sat in enumerate(self.sat_list[const]):
                res_dict[const][sat] = dict()

                for iFreq, datatype in enumerate(self._metadata["CODES"][const]):
                    res_dict[const][sat][datatype] = residual_vec[iFreq * n_sats + iSat]
        return res_dict


# TODO: criar uma class LSQ_Engine mãe com as funções high level , e depois fazer duas filhas para pos e vel

class LSQ_Engine_Vel:
    def __init__(self, satellite_list, metadata, epoch, obs_data, reconstructor):
        self.y_vec = None  # observation vector
        self.design_mat = None  # design matrix
        self.weight_mat = None  # weight matrix

        self._metadata = metadata
        self.doppler_measurements = self._metadata["DOPPLER"]
        self.constellations = metadata["CONSTELLATIONS"]
        self.sat_list = dict()
        for sat in satellite_list:
            if sat.sat_system not in self.sat_list:
                self.sat_list[sat.sat_system] = list()
            self.sat_list[sat.sat_system].append(sat)

        self._initialize_matrices()
        self._build_lsq(epoch, obs_data, reconstructor)

    def _initialize_matrices(self):
        # TODO: por agora apenas single constellation GPS+GAL
        # different cases:
        #   * Single Frequency Single Constellation: States: position, clock
        #       -> Design matrix is [m_S,4]
        #       -> Observation vector is [m_S,1]
        #   * Dual Frequency Single Constellation: States: position, clock, iono
        #       -> Design matrix is [2*m_S,4+m_S]
        #       -> Observation vector is [2*m_S,1]
        #   * Single Frequency Dual Constellation: States: position, clock, isb
        #       -> Design matrix is [m_S1+m_S2,4+1]
        #       -> Observation vector is [m_S1+m_S2,1]
        #   * Dual Frequency Dual Constellation: States: position, clock, iono, isb
        #       -> Design matrix is [2*(m_S1+m_S2), 4+m_S1+m_S2+1]
        #       -> Observation vector is [2*(m_S1+m_S2),1]
        #   * (Single + Dual) Frequency Dual Constellation:
        #       -> States: position, clock, iono, isb. Design matrix is Design matrix is [2*m_S1+m_S2, 4+m_S1+1]
        #       -> Observation vector is [2*m_S1+m_S2,1]
        #  When tropo is added, a new state is incremented, that is 4 -> 4 + 1 in the dimensions above

        n_observables = 0  # number of rows
        n_states = 3  # number of columns (default is 3 - velocity)

        for const in self.constellations:
            n_sats = len(self.sat_list[const])
            n_observables += n_sats
            n_states += 1  # receiver clock rate for each constellation

        self.y_vec = np.zeros(n_observables)
        self.design_mat = np.zeros((n_observables, n_states))
        self.weight_mat = np.eye(n_observables)

    @staticmethod
    def compute_residual(sat, epoch, doppler_datatype, obs_data, reconstructor):

        # get observable and compute predicted observable
        obs = obs_data.get_observable(sat, doppler_datatype)

        # transform Doppler to pseudorange rate
        wavelength = SPEED_OF_LIGHT / doppler_datatype.freq_value  # in meters
        obs_range_rate = -wavelength * float(obs)  # in m/s

        predicted_obs = reconstructor.compute(sat, epoch, doppler_datatype)

        # prefit residuals (measured gnss_models - predicted gnss_models)
        prefit_residuals = obs_range_rate - predicted_obs
        return prefit_residuals

    def _build_lsq(self, epoch, obs_data, reconstructor):
        """build the LS matrices y_vec, design_mat, weight_mat"""

        obs_offset = 0
        const_offset = 0
        for iConst, const in enumerate(self.constellations):

            n_sats = len(self.sat_list[const])
            doppler_datatype = self.doppler_measurements[const][0]

            for iSat, sat in enumerate(self.sat_list[const]):
                los = -reconstructor.get_unit_line_of_sight(sat)
                residual = self.compute_residual(sat, epoch, doppler_datatype, obs_data, reconstructor)

                # filling the LS matrices
                self.y_vec[obs_offset + iSat] = residual
                self.design_mat[obs_offset + iSat][0:3] = -los  # velocity
                self.design_mat[obs_offset + iSat][3 + const_offset] = 1

                # Weight matrix -> as 1/(obs_std^2)
                self.weight_mat[obs_offset + iSat, obs_offset + iSat] = 1.0  # / (reconstructor.get_obs_std(sat, datatype)**2)
            obs_offset += n_sats
            const_offset += 1

    def solve_ls(self, state):
        """solves the LS problem for this iteration"""
        try:
            solver = WeightedLeastSquares(self.y_vec, self.design_mat, W=self.weight_mat)
            solver.solve()
        except (AttributeError, np.linalg.LinAlgError) as e:

            # possible error in the numpy.linalg.inv() function -> solution not possible
            raise SolverError(e)

        # update state vector with incremental dX
        dX = solver.get_solution()
        cov = solver.get_cov()

        # apply dX corrections to the state vector
        self.apply_corrections(state, dX, cov)

        # form the post-fit residuals
        post_fit = self.y_vec - self.design_mat[:, 0:3] @ dX[0:3]
        norm = np.linalg.norm(post_fit)

        # form residual dicts for output
        pre_fit_dict = self.get_residuals(self.y_vec)
        post_fit_dict = self.get_residuals(post_fit)

        return pre_fit_dict, post_fit_dict, norm

    def apply_corrections(self, state, dX, cov):
        """applies corrections to the state vector"""

        state.velocity += dX[0:3]  # in m/s
        state.cov_velocity = cov[0:3, 0:3]  # in (m/s)^2

        for iConst, const in enumerate(self.constellations):
            # receiver clock drift [dimensionless]
            state.clock_bias_rate[const] += dX[iConst + 3] / constants.SPEED_OF_LIGHT
            state.cov_clock_bias_rate[const] = float(cov[iConst + 3, iConst + 3]) / (constants.SPEED_OF_LIGHT ** 2)

    def get_residuals(self, residual_vec):

        res_dict = dict()
        iSat = 0
        for const in self.constellations:
            res_dict[const] = dict()

            for sat in self.sat_list[const]:
                res_dict[const][sat] = residual_vec[iSat]
                iSat += 1

        return res_dict
