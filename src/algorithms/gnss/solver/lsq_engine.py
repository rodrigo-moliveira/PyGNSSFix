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
        self._initialize_ls(len(satellite_list))
        self.satellite_list = satellite_list

    def _initialize_ls(self, n_sats):
        nr_consts = len(self._metadata["CONSTELLATIONS"])

        if nr_consts == 1:
            const = self._metadata["CONSTELLATIONS"][0]

            if self._metadata["MODEL"][const] == EnumModel.SINGLE_FREQ:
                self.y_vec = np.zeros(n_sats)  # gnss_obs vector <=> prefit residuals
                self.design_mat = np.ones((n_sats, 4))  # geometry matrix
                self.weight_mat = np.eye(n_sats)  # diagonal weight matrix
            else:
                self.y_vec = np.zeros(n_sats * 2)  # gnss_obs vector <=> prefit residuals
                self.weight_mat = np.eye(n_sats * 2)  # diagonal weight matrix
                geometry = np.ones((n_sats, 4))  # geometry matrix (state + clock)
                ionoMatrix1 = np.eye(n_sats)  # iono matrices for freq 1 and freq 2
                factor = (self._metadata["CODES"][const][0].freq.freq_value /
                          self._metadata["CODES"][const][1].freq.freq_value) ** 2
                ionoMatrix2 = np.eye(n_sats) * factor

                self.design_mat = np.block([
                    [geometry, ionoMatrix1],
                    [geometry, ionoMatrix2]]
                )
        else:
            raise PVTComputationFail("Currently only 1 constellation is possible!")

    def build_lsq(self, epoch, obs_data, reconstructor, nav_data):
        n_sats = len(self.satellite_list)
        constellations = self._metadata["CONSTELLATIONS"]
        nr_consts = len(constellations)

        if nr_consts != 1:
            raise PVTComputationFail("Currently only 1 constellation is possible!")

        iFreq = 0
        const = constellations[0]
        for datatype in self._metadata["CODES"][const]:
            if datatype is None:  # TODO or datatype is not pseudorange
                continue

            iSat = 0
            for sat in self.satellite_list:
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
                self.y_vec[iFreq * n_sats + iSat] = prefit_residuals.value
                self.design_mat[iFreq * n_sats + iSat][0:3] = line_sight

                # Weight matrix -> sigma = 1 / e^{-elevation}
                self.weight_mat[iFreq * n_sats + iSat][iFreq * n_sats + iSat] = reconstructor.get_weight(sat)

                iSat += 1
            iFreq += 1

    def solve_ls(self, state):
        # TODO this function needs to be cleaned
        constellations = self._metadata["CONSTELLATIONS"]
        const = constellations[0]

        # solve LS problem for this iteration
        try:
            solver = WeightedLeastSquares(self.y_vec, self.design_mat, W=self.weight_mat)
            solver.solve()

            # covariance matrix of the LS estimator
            cov = np.linalg.inv(self.design_mat.T @ self.weight_mat @ self.design_mat)
            dop_matrix = np.linalg.inv(self.design_mat.T @ self.design_mat)

        except (AttributeError, np.linalg.LinAlgError) as e:

            # possible error in the numpy.linalg.inv() function -> solution not possible
            raise PVTComputationFail(e)

        # update state vector with incremental dX
        dX = solver.get_solution()
        state.position += dX[0:3]
        state.clock_bias = dX[3] / constants.SPEED_OF_LIGHT  # receiver clock in seconds

        # if iono is estimated
        if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
            for iSat, sat in enumerate(self.satellite_list):
                state.iono[const][sat] = dX[iSat + 4]

        # if isb is estimated
        # ...

        # get post-fit residuals
        post_fit = self.y_vec - self.design_mat[:, 0:3] @ dX[0:3]
        norm = np.linalg.norm(post_fit)

        # form residual dicts for output
        pre_fit_dict = self.get_residuals(self.y_vec)
        post_fit_dict = self.get_residuals(post_fit)

        # unpack covariance matrices
        state.cov_position = cov[0:3, 0:3]
        state.cov_clock_bias = cov[3, 3] / (constants.SPEED_OF_LIGHT ** 2)  # in seconds^2
        # if iono is estimated
        if self._metadata["MODEL"][const] == EnumModel.DUAL_FREQ:
            for iSat, sat in enumerate(self.satellite_list):
                state.cov_iono[const][sat] = dX[iSat + 4]

        return pre_fit_dict, post_fit_dict, dop_matrix, norm

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
