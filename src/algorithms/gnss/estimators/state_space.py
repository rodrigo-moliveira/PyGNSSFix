import numpy as np

from src import constants
from src.data_mng.container import Container


class GnssStateSpace(Container):
    __slots__ = ["position", "velocity", "clock_bias", "iono", "isb", "date", "_info"]

    def __init__(self, **kwargs):
        super().__init__()

        # mandatory fields (initialized to 0)
        self.position = kwargs.get("position", np.array([0, 0, 0]))  # with fallback
        self.clock_bias = kwargs.get("clock_bias", 0)  # with fallback

        # optional solve-for variables (initialized to None)
        self.isb = kwargs.get("ISB", None)  # fallback to None (not estimated)
        self.iono = kwargs.get("iono", None)  # fallback to None (not estimated)
        self.velocity = kwargs.get("velocity", None)  # fallback to None (not estimated)

        # time tag
        self.date = kwargs.get("date", None)

        self._info = dict()

    def __str__(self):
        _str = f"position = {self.position}, " \
               f"clock bias = {self.clock_bias}"
        if self.velocity:
            _str += f", velocity = {self.velocity}"
        if self.isb:
            _str += f", ISB = {self.isb}"
        if self.iono:
            _str += f", iono = {self.iono}"

        return f'{type(self).__name__}[{str(self.date)}]({_str})'

    def __repr__(self):
        return str(self)

    def __len__(self):

        # n = receiver position (3) + clock bias (1)
        n = len(self.position) + 1

        if self.iono is not None:
            n += len(self.iono)
        if self.isb is not None:
            n += 1
        if self.velocity is not None:
            n += len(self.velocity)

        return n

    def add_solver_info(self, arg, val):
        self._info[arg] = val

    def get_solver_info(self, arg):
        if arg in self._info:
            return self._info[arg]
        return None

    def get_estimables(self):
        estimables = set()
        for est in self.__slots__:
            if est != "_info" and est != "date":
                if getattr(self, est, None) is not None:
                    estimables.add(est)

            if est == "_info":
                info = getattr(self, est, None)
                if info is not None:
                    if "geometry" in info.keys():
                        estimables.add("satellite_azel")  # save satellite azimuth and elevation
                    if "dop" in info.keys():
                        estimables.add("dop")
                    if "prefit_residuals" in info.keys():
                        estimables.add("prefit_residuals")
                    if "postfit_residuals" in info.keys():
                        estimables.add("postfit_residuals")
        return estimables

    def get_header(self, estimable):
        if estimable == "position":
            return "Epoch,X_ECEF[m],Y_ECEF[m],Z_ECEF[m],cov_XX[m^2],cov_YY[m^2],cov_ZZ[m^2],cov_XY[m^2],cov_XZ[m^2]," \
                   "cov_YZ[m^2]"
        if estimable == "clock_bias":
            return "Epoch,clock_bias[s],cov[s^2]"
        elif estimable == "prefit_residuals":
            return "Epoch,sat,prefit_residuals_i[m^2]"
        elif estimable == "postfit_residuals":
            return "Epoch,sat,postfit_residuals_i[m^2]"
        elif estimable == "satellite_azel":
            return "Epoch,sat,azimuth[deg],elevation[deg]"
        elif estimable == "dop":
            return "Epoch,DOP_X[m],DOP_Y[m],DOP_Z[m],DOP_T[m]"
        else:
            raise ValueError(f"Undefined header due to unknown estimable {estimable}")

    def export_to_file(self, estimable):
        if estimable == "position":
            cov = self._info["cov"]
            cov_xx = cov[0, 0]
            cov_yy = cov[1, 1]
            cov_zz = cov[2, 2]
            cov_xy = cov[0, 1]
            cov_xz = cov[0, 2]
            cov_yz = cov[1, 2]
            return f"{self.position[0]},{self.position[1]},{self.position[2]},{cov_xx},{cov_yy},{cov_zz},{cov_xy}," \
                   f"{cov_xz},{cov_yz}"

        if estimable == "clock_bias":
            cov = self._info["cov"]
            cov_t = cov[3, 3] / (constants.SPEED_OF_LIGHT**2)  # in seconds^2
            return f"{self.clock_bias},{cov_t}"

        elif estimable == "prefit_residuals" or estimable == "postfit_residuals":
            residuals = self._info[estimable]
            sat_list = self._info["sat_list"]
            data = []
            if len(sat_list) > 0:
                for sat, res in zip(sat_list, residuals):
                    data.append(f"{sat},{res}")
            return data

        elif estimable == "satellite_azel":
            geometry = self._info["geometry"]
            sat_list = self._info["sat_list"]
            data = []
            if len(sat_list) > 0:
                for sat in sat_list:
                    el = geometry.get("el", sat) * constants.RAD2DEG
                    az = geometry.get("az", sat) * constants.RAD2DEG
                    data.append(f"{sat},{az},{el}")
            return data

        elif estimable == "dop":
            dop = self._info["dop"]
            return f"{dop[0, 0]},{dop[1, 1]},{dop[2, 2]},{dop[3, 3]}"

        else:
            raise ValueError(f"Undefined header due to unknown estimable {estimable}")
