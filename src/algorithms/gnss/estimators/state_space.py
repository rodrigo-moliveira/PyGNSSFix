import numpy as np

from src.data_mng.container import Container
from src.io.config.enums import EnumModel


class GnssStateSpace(Container):
    __states__ = ["position", "clock_bias", "iono", "isb"]
    __slots__ = __states__ + ["epoch", "_info"]

    def __init__(self, metadata, **kwargs):
        super().__init__()

        # dict to store state information
        self._info = dict()

        # initialize system solve-for states
        self._init_states(metadata, kwargs)

        # epoch
        self.epoch = kwargs.get("epoch", None)

    def _init_states(self, metadata, kwargs):
        _states = ["position", "clock_bias"]  # mandatory states
        sat_list = kwargs.get("sat_list", [])

        # position
        self.position = np.array(kwargs.get("position", [0, 0, 0]))

        # clock
        self.clock_bias = kwargs.get("clock_bias", 0)

        # iono (optional -> in case there are 2 frequencies for this constellation)
        self.iono = dict()
        for constellation in metadata["CONSTELLATIONS"]:
            if metadata["MODEL"][constellation] == EnumModel.DUAL_FREQ:
                self.iono[constellation] = dict()
                for sat in sat_list:
                    if sat.sat_system == constellation:
                        self.iono[constellation][sat] = 0  # initialize iono for this satellite

                # remove this constellation if no satellite was actually inserted
                if len(self.iono[constellation]) == 0:
                    self.iono.pop(constellation)
        if len(self.iono) >= 1:
            _states.append("iono")
        else:
            self.iono = None

        # isb (optional -> in case there are 2 constellations)
        if len(metadata["CONSTELLATIONS"]) > 1:
            self.isb = kwargs.get("isb", 0)
            _states.append("isb")
        else:
            self.isb = None

        self.add_additional_info("states", _states)

    def __str__(self):
        _str = f"position = {self.position}, " \
               f"clock bias = {self.clock_bias}"
        if "isb" in self.get_additional_info("states"):
            _str += f", ISB = {self.isb}"
        if "iono" in self.get_additional_info("states"):
            _str += f", iono = {self.iono}"

        return f'{type(self).__name__}[{str(self.date)}]({_str})'

    def __repr__(self):
        return str(self)

    def add_additional_info(self, arg, val):
        self._info[arg] = val

    def get_additional_info(self, arg):
        if arg in self._info:
            return self._info[arg]
        return None

    def get_exportable_lst(self):
        exportable_lst = list(self.get_additional_info("states"))

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
