import numpy as np

from src import constants
from src.algorithms.gnss.estimators.weighted_ls import WeightedLeastSquares
from src.errors import PVTComputationFail
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

        # add tropo...
        # pass

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

        # prefit residuals (measured gnss_obs - predicted gnss_obs)
        prefit_residuals = obs - predicted_obs

        # get LOS vector w.r.t. ECEF frame (column in geometry matrix)
        line_sight = reconstructor.get_unit_line_of_sight(sat)

        return prefit_residuals.value, line_sight

    def _build_lsq(self, epoch, obs_data, reconstructor, nav_data):
        """build the LS matrices y_vec, design_mat, weight_mat"""

        iono_offset = 0
        obs_offset = 0
        for iConst, const in enumerate(self.constellations):

            n_sats = len(self.sat_list[const])
            for iFreq, datatype in enumerate(self._metadata["CODES"][const]):

                for iSat, sat in enumerate(self.sat_list[const]):
                    residual, los = self.compute_residual_los(nav_data, sat, epoch, datatype, obs_data, reconstructor)

                    # filling the LS matrices
                    self.y_vec[obs_offset + iSat] = residual
                    self.design_mat[obs_offset + iSat][0:3] = los  # position
                    self.design_mat[obs_offset + iSat, 3] = 1.0  # clock
                    if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
                        factor = (self._metadata["CODES"][const][0].freq.freq_value / datatype.freq.freq_value) ** 2
                        self.design_mat[obs_offset + iSat, 4 + iono_offset + iSat] = 1.0 * factor  # iono
                    if iConst > 0:
                        self.design_mat[obs_offset + iSat, -1] = 1.0  # ISB

                    # Weight matrix -> as 1/(obs_std^2)
                    self.weight_mat[obs_offset + iSat, obs_offset + iSat] = \
                        1 / (reconstructor.get_obs_std(sat, datatype)**2)

                obs_offset += n_sats
            iono_offset += n_sats

    def solve_ls(self, state):
        """solves the LS problem for this iteration"""
        try:
            solver = WeightedLeastSquares(self.y_vec, self.design_mat, W=self.weight_mat)
            solver.solve()
            dop_matrix = np.linalg.inv(self.design_mat.T @ self.design_mat)

        except (AttributeError, np.linalg.LinAlgError) as e:

            # possible error in the numpy.linalg.inv() function -> solution not possible
            raise PVTComputationFail(e)
        print("solved with success")
        exit()
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
        const = self.constellations[0]

        state.position += dX[0:3]
        state.clock_bias += dX[3] / constants.SPEED_OF_LIGHT  # receiver clock in seconds

        # if iono is estimated
        if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
            for iSat, sat in enumerate(self.sat_list[const]):
                state.iono[sat] += float(dX[iSat + 4])

        # if isb is estimated
        # ...

        # unpack covariance matrices
        state.cov_position = np.array(cov[0:3, 0:3])
        state.cov_clock_bias = float(cov[3, 3]) / (constants.SPEED_OF_LIGHT ** 2)  # in seconds^2

        # if iono is estimated
        if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
            for iSat, sat in enumerate(self.sat_list[const]):
                state.cov_iono[sat] = cov[iSat + 4, iSat + 4]

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


"""if self.n_consts == 1:
const = self._metadata["CONSTELLATIONS"][0]

if self._metadata["MODEL"][const] == EnumModel.SINGLE_FREQ:
    self.y_vec = np.zeros(self.n_sats)  # gnss_obs vector <=> prefit residuals
    self.design_mat = np.ones((self.n_sats, 4))  # geometry matrix
    self.weight_mat = np.eye(self.n_sats)  # diagonal weight matrix
else:
    self.y_vec = np.zeros(self.n_sats * 2)  # gnss_obs vector <=> prefit residuals
    self.weight_mat = np.eye(self.n_sats * 2)  # diagonal weight matrix
    geometry = np.ones((self.n_sats, 4))  # geometry matrix (state + clock)
    ionoMatrix1 = np.eye(self.n_sats)  # iono matrices for freq 1 and freq 2
    factor = (self._metadata["CODES"][const][0].freq.freq_value /
              self._metadata["CODES"][const][1].freq.freq_value) ** 2
    ionoMatrix2 = np.eye(self.n_sats) * factor

    self.design_mat = np.block([
        [geometry, ionoMatrix1],
        [geometry, ionoMatrix2]]
    )
"""