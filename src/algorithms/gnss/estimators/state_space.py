import numpy as np

from src.data_mng.container import Container


class GnssStateSpace(Container):
    __states__ = ["position", "velocity", "clock_bias", "iono", "isb"]
    __slots__ = __states__ + ["epoch", "_info"]

    def __init__(self, **kwargs):
        super().__init__()

        # initialize system solve-for states
        self.position = np.array(kwargs.get("position")) if "position" in kwargs else None
        self.clock_bias = kwargs.get("clock_bias", None)
        self.isb = kwargs.get("isb", None)
        self.iono = kwargs.get("iono", None)
        self.velocity = np.array(kwargs.get("velocity")) if "velocity" in kwargs else None

        # epoch
        self.epoch = kwargs.get("epoch", None)

        # additional info
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
        # return length of state vector (number of solve-for parameters)
        n = 0

        if self.position is not None:
            n += len(self.position)
        if self.clock_bias is not None:
            n += 1
        if self.iono is not None:
            n += len(self.iono)
        if self.isb is not None:
            n += 1
        if self.velocity is not None:
            n += len(self.velocity)

        return n

    def add_additional_info(self, arg, val):
        self._info[arg] = val

    def get_additional_info(self, arg):
        if arg in self._info:
            return self._info[arg]
        return None

    def get_exportable_lst(self):
        exportable_lst = list()
        for ext in self.__states__:
            if getattr(self, ext, None) is not None:
                exportable_lst.append(ext)

        if self.epoch is not None:
            exportable_lst.append("time")

        # additional info
        if "geometry" in self._info.keys():
            exportable_lst.append("satellite_azel")
        if "dop_ecef" in self._info.keys():
            exportable_lst.append("dop_ecef")
        if "dop_local" in self._info.keys():
            exportable_lst.append("dop_local")
        if "prefit_residuals" in self._info.keys():
            exportable_lst.append("prefit_residuals")
        if "postfit_residuals" in self._info.keys():
            exportable_lst.append("postfit_residuals")
        return exportable_lst
