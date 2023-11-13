import numpy as np

from src import constants
from src.algorithms.gnss.estimators.weighted_ls import WeightedLeastSquares
from src.errors import PVTComputationFail
from src.io.config.enums import EnumModel


class LSQ_Engine:
    def __init__(self, satellite_list, metadata):
        self.y_vec = None  # observation vector
        self.design_mat = None  # design matrix
        self.weight_mat = None  # weight matrix

        self._metadata = metadata
        self.satellite_list = satellite_list

        self.n_sats = len(satellite_list)
        self.n_consts = len(self._metadata["CONSTELLATIONS"])

        self._initialize_ls()

    def _initialize_ls(self):
        if self.n_consts == 1:
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
        else:
            raise PVTComputationFail("Currently only 1 constellation is possible!")

    def build_lsq(self, epoch, obs_data, reconstructor, nav_data):
        """build the LS matrices"""
        if self.n_consts != 1:
            raise PVTComputationFail("Currently only 1 constellation is possible!")

        const = self._metadata["CONSTELLATIONS"][0]
        for iFreq, datatype in enumerate(self._metadata["CODES"][const]):
            for iSat, sat in enumerate(self.satellite_list):
                # get observable
                obs = obs_data.get_observable(sat, datatype)

                # fetch valid navigation message (closest to the current epoch)
                nav_message = nav_data.get_closest_message(sat, epoch)

                # compute predicted gnss_obs
                predicted_obs = reconstructor.compute(nav_message, sat, epoch, datatype)

                # prefit residuals (measured gnss_obs - predicted gnss_obs)
                prefit_residuals = obs - predicted_obs

                # get LOS vector w.r.t. ECEF frame (column in geometry matrix)
                line_sight = reconstructor.get_unit_line_of_sight(sat)

                # filling the corresponding entry of data vectors for the Least Squares quantities
                self.y_vec[iFreq * self.n_sats + iSat] = prefit_residuals.value
                self.design_mat[iFreq * self.n_sats + iSat][0:3] = line_sight

                # Weight matrix -> sigma = 1 / e^{-elevation}
                self.weight_mat[iFreq * self.n_sats + iSat][iFreq * self.n_sats + iSat] = reconstructor.get_weight(sat)

    def solve_ls(self, state):
        """solves the LS problem for this iteration"""
        try:
            solver = WeightedLeastSquares(self.y_vec, self.design_mat, W=self.weight_mat)
            solver.solve()
            dop_matrix = np.linalg.inv(self.design_mat.T @ self.design_mat)

        except (AttributeError, np.linalg.LinAlgError) as e:

            # possible error in the numpy.linalg.inv() function -> solution not possible
            raise PVTComputationFail(e)

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
        constellations = self._metadata["CONSTELLATIONS"]
        const = constellations[0]

        state.position += dX[0:3]
        state.clock_bias += dX[3] / constants.SPEED_OF_LIGHT  # receiver clock in seconds

        # if iono is estimated
        if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
            for iSat, sat in enumerate(self.satellite_list):
                state.iono[const][sat] += float(dX[iSat + 4])

        # if isb is estimated
        # ...

        # unpack covariance matrices
        state.cov_position = np.array(cov[0:3, 0:3])
        state.cov_clock_bias = float(cov[3, 3]) / (constants.SPEED_OF_LIGHT ** 2)  # in seconds^2

        # if iono is estimated
        if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
            for iSat, sat in enumerate(self.satellite_list):
                state.cov_iono[const][sat] = cov[iSat + 4, iSat + 4]

    def get_residuals(self, residual_vec):
        constellations = self._metadata["CONSTELLATIONS"]
        n_sats = len(self.satellite_list)

        res_dict = dict()
        for const in constellations:
            res_dict[const] = dict()

            for iSat, sat in enumerate(self.satellite_list):
                res_dict[const][sat] = dict()

                for iFreq, datatype in enumerate(self._metadata["CODES"][const]):
                    res_dict[const][sat][datatype] = residual_vec[iFreq * n_sats + iSat]
        return res_dict
