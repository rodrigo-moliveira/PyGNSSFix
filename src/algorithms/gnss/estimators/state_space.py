import numpy as np

from src.data_mng.container import Container


class GnssStateSpace(Container):
    __slots__ = ["position", "velocity", "clock_bias", "iono", "ISB", "date", "_info"]

    def __init__(self, **kwargs):
        super().__init__()

        # mandatory fields (initialized to 0)
        self.position = kwargs.get("position", np.array([0, 0, 0]))  # with fallback
        self.clock_bias = kwargs.get("clock_bias", 0)  # with fallback

        # optional solve-for variables (initialized to None)
        self.ISB = kwargs.get("ISB", None)  # fallback to None (not estimated)
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
        if self.ISB:
            _str += f", ISB = {self.ISB}"
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
        if self.ISB is not None:
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
                    if "DOP" in info.keys():
                        estimables.add("DOP")
                    if "prefit_residuals" in info.keys():
                        estimables.add("residuals")
                    if "postfit_residuals" in info.keys():
                        estimables.add("residuals")
                    if "rms" in info.keys():
                        estimables.add("rms")
        return estimables

    def get_header(self, estimable):
        if estimable == "position":
            return "Epoch,X_ECEF[m],Y_ECEF[m],Z_ECEF[m],cov_XX[m^2],cov_YY[m^2],cov_ZZ[m^2],cov_XZ[m^2],cov_YZ[m^2],cov_ZZ[m^2]"
        if estimable == "position":
            return "Epoch,clock_bias[s],cov[s^2]"
        elif estimable == "residuals":
            return "Epoch,prefit_residuals[TBC],postfit_residuals[TBC]"
        elif estimable == "rms":
            return "Epoch,rms[TBC]"
        elif estimable == "satellite_azel":
            return "Epoch,sat,azimuth[deg],elevation[deg]"
        else:
            raise ValueError(f"Undefined header due to unknown estimable {estimable}")

    def export_to_file(self, directory):
        pass
