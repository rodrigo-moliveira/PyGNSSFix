import numpy as np

from src import constants
from src.algorithms.gnss.estimators.weighted_ls import WeightedLeastSquares
from src.errors import PVTComputationFail
from src.io.config.enums import EnumModel


class LSQ_Engine:
    def __init__(self, satellite_list, model_config, services):
        self.y_vec = None       # observation vector
        self.design_mat = None  # design matrix
        self.weight_mat = None  # weight matrix
        self._available = []

        self._initialize_ls(len(satellite_list), model_config, services)
        self.services = services
        self.model_config = model_config
        self.satellite_list = satellite_list

    def _initialize_ls(self, n_sats, model_config, services):
        # least squares arrays

        # TODO: Currently this is GPS specific
        if model_config["GPS"] == EnumModel.SINGLE_FREQ:
            self.y_vec = np.zeros(n_sats)           # gnss_obs vector <=> prefit residuals
            self.design_mat = np.ones((n_sats, 4))  # geometry matrix
            self.weight_mat = np.eye(n_sats)        # diagonal weight matrix
        else:
            self.y_vec = np.zeros(n_sats * 2)       # gnss_obs vector <=> prefit residuals
            self.weight_mat = np.eye(n_sats * 2)    # diagonal weight matrix
            geometry = np.ones((n_sats, 4))         # geometry matrix (state + clock)
            ionoMatrix1 = np.eye(n_sats)            # iono matrices for freq 1 and freq 2
            ionoMatrix2 = np.eye(n_sats) * (
                    services["GPS"][0].freq.freq_value / services["GPS"][1].freq.freq_value) ** 2

            self.design_mat = np.block([
                [geometry, ionoMatrix1],
                [geometry, ionoMatrix2]]
            )

    def build_lsq(self, epoch, obs_data, reconstructor, system_geometry, nav_data):
        n_sats = len(self.satellite_list)

        iFreq = 0
        for datatype in self.services["GPS"]:
            if datatype is None:
                continue

            iSat = 0
            for sat in self.satellite_list:
                # get observable
                obs = obs_data.get_observable(sat, datatype)

                # fetch valid navigation message (closest to the current epoch)
                nav_message = nav_data.get_sat_data_for_epoch(sat, epoch)

                # compute predicted gnss_obs
                predicted_obs = reconstructor.compute(nav_message, sat, epoch, datatype)

                # prefit residuals (measured gnss_obs - predicted gnss_obs)
                prefit_residuals = obs - predicted_obs

                # get LOS vector w.r.t. ECEF frame (column in geometry matrix)
                line_sight = system_geometry.get_unit_line_of_sight(sat)

                # filling the corresponding entry of data vectors for the Least Squares quantities
                self.y_vec[iFreq * n_sats + iSat] = prefit_residuals.value
                self.design_mat[iFreq * n_sats + iSat][0:3] = line_sight

                # Weight matrix -> sigma = 1 / e^{-elevation}
                self.weight_mat[iFreq * n_sats + iSat][iFreq * n_sats + iSat] = self.get_weight(system_geometry, sat)

                iSat += 1
            iFreq += 1

    def solve_ls(self, state):
        n_sats = len(self.satellite_list)

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

        dX = solver.get_solution()

        # update state vector with incremental dX
        state.position += dX[0:3]
        state.clock_bias = dX[3] / constants.SPEED_OF_LIGHT  # receiver clock in seconds

        # if iono is estimated
        if self.model_config["GPS"] == EnumModel.DUAL_FREQ:
            state.iono = [dX[i + 4] for i in range(n_sats)]

        # if isb is estimated
        # ...

        # get post-fit residuals
        post_fit = self.y_vec - self.design_mat[:, 0:3] @ dX[0:3]

        return post_fit, cov, dop_matrix

    @classmethod
    def get_weight(cls, system_geometry, sat):
        # TODO: need to add here the user defined sigmas as a multiplication factor
        # "obs_std", and can add the possibility of this mask as well.
        sigma_elevation = np.e ** (-system_geometry.get("el", sat))
        w = (1 / sigma_elevation) ** 2

        return w
