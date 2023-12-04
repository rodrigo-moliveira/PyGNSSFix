import numpy as np

from src.data_mng.container import Container
from src.io.config.enums import EnumModel
from src.models.gnss_obs.clock_obs import compute_ggto


class GnssStateSpace(Container):
    __states__ = ["position", "clock_bias", "iono", "isb"]
    __covs__ = ["cov_position", "cov_clock_bias", "cov_iono", "cov_isb"]
    __slots__ = __states__ + __covs__ + ["epoch", "_info"]

    def __init__(self, metadata=None, position=None, clock_bias=None, epoch=None, sat_list=None):
        super().__init__()

        # dict to store state information
        self._info = dict()

        # initialize system solve-for states
        self._init_states(metadata, position, clock_bias, sat_list)

        # epoch
        self.epoch = epoch

    def clone(self):
        _states = self.get_additional_info("states")
        state = GnssStateSpace(position=self.position, clock_bias=self.clock_bias,
                               epoch=self.epoch)

        if "iono" in _states:
            for sat, iono in self.iono.items():
                state.iono[sat] = iono

        if "isb" in _states:
            state.isb = self.isb

        state.add_additional_info("states", _states)
        state.add_additional_info("clock_master", self.get_additional_info("clock_master"))
        state.add_additional_info("clock_slave", self.get_additional_info("clock_slave"))

        return state

    def _init_states(self, metadata, position, clock_bias, sat_list):
        _states = ["position", "clock_bias"]  # mandatory states

        # position (with default to [0, 0, 0])
        self.position = np.array(position if position is not None else [0, 0, 0])
        self.cov_position = np.zeros((3, 3))

        # clock with default to 0
        self.clock_bias = clock_bias if clock_bias is not None else 0.0
        self.cov_clock_bias = 0.0

        # iono dict (if dual-frequency mode is selected) with default to None
        self.iono = dict()
        self.cov_iono = dict()
        if metadata is not None and sat_list is not None:
            for constellation in metadata["CONSTELLATIONS"]:
                if metadata["MODEL"][constellation] == EnumModel.DUAL_FREQ:
                    # initialize iono vector for the available satellites of this constellation
                    for sat in sat_list:
                        if sat.sat_system == constellation:
                            self.iono[sat] = 0.0
                            self.cov_iono[sat] = 0.0
        if len(self.iono) >= 1:
            _states.append("iono")

        # isb (optional -> in case there are 2 constellations)
        self.isb = None
        self.cov_isb = None
        if metadata is not None and len(metadata["CONSTELLATIONS"]) > 1:
            self.isb = 0.0
            self.cov_isb = 0.0
            self.add_additional_info("clock_master", metadata["CONSTELLATIONS"][0])
            self.add_additional_info("clock_slave", metadata["CONSTELLATIONS"][1])
            _states.append("isb")
        else:
            self.add_additional_info("clock_master", None)
            self.add_additional_info("clock_slave", None)

        self.add_additional_info("states", _states)

    def __str__(self):
        _states = self.get_additional_info("states")
        _str = f"position = {self.position}, " \
               f"clock bias = {self.clock_bias}"
        if "isb" in _states:
            _str += f", ISB = {self.isb}"
        if "iono" in _states:
            _str += f", iono = {self.iono}"

        return f'{type(self).__name__}[{str(self.epoch)}]({_str})'

    def __repr__(self):
        return str(self)

    def get_clock_bias(self, constellation, time_correction):
        if "isb" in self.get_additional_info("states"):
            if constellation == self.get_additional_info("clock_master"):
                clock = self.clock_bias
            else:
                isb = self.get_isb(constellation, time_correction)
                clock = self.clock_bias + isb
        else:
            clock = self.clock_bias
        return clock

    def get_isb(self, constellation, time_correction):
        ggto = compute_ggto(time_correction, self.epoch)
        # if estimate_ggto is false
        #       if constellation is slave
        #           if constellation is GPS
        #               ggto = -ggt0
        #       return self.isb - ggto
        # else
        return self.isb

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
