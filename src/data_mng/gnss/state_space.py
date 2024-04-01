import numpy as np

from src.data_mng import Container
from src.io.config import EnumModel, config_dict
from src.models.gnss_models import compute_ggto


class GnssStateSpace(Container):
    """
    This class contains the state space vector with all variables to be estimated in a GNSS run,
    and the corresponding covariance values

    States:
        * epoch (Epoch): Epoch object for this time instant
        * position (numpy.ndarray): cartesian components of ECEF position vector (XYZ)
        * velocity (numpy.ndarray): cartesian components of ECEF velocity vector (XYZ)
        * clock_bias (float): estimated clock bias
        * iono (dict): dictionary with satellites as keys and estimated ionosphere delay as values
        * tropo_wet (float): estimated wet component of troposphere delay
        * isb (float): estimated inter system bias
        * clock_bias_rate (dict): dictionary with constellations as keys and clock bias rates as values

    Some auxiliary information is saved in the `_info` attribute dict, namely:
        * _info["states"] provides a list with all valid states to be estimated, depending on the setup
        * _info["clock_master"] provides the master constellation, when the ISB state is active
        * _info["clock_slave"] provides the slave constellation, when the ISB state is active

    Additionally, the GNSS solver will also save the following information in the _info dict:
        * geometry
        * dop_ecef
        * dop_local
        * prefit_residuals
        * postfit_residuals
    """
    __states__ = ["position", "velocity", "clock_bias", "iono", "tropo_wet", "isb", "clock_bias_rate"]
    __covs__ = ["cov_position", "cov_velocity", "cov_clock_bias", "cov_iono", "cov_tropo_wet", "cov_isb",
                "cov_clock_bias_rate"]
    __slots__ = __states__ + __covs__ + ["epoch", "_info"]

    def __init__(self, metadata=None, position=None, velocity=None, clock_bias=None, epoch=None,
                 clock_bias_rate=None, sat_list=None):
        super().__init__()

        # dict to store state information
        self._info = dict()

        # initialize system solve-for states
        self._init_states(metadata, position, velocity, clock_bias, clock_bias_rate, sat_list)

        # epoch
        self.epoch = epoch

    def clone(self):
        """Returns a copy of this GnssStateSpace object. All state variables are deep-copied."""
        _states = self.get_additional_info("states")
        state = GnssStateSpace(position=self.position, velocity=self.velocity, clock_bias=self.clock_bias,
                               epoch=self.epoch, clock_bias_rate=self.clock_bias_rate)

        if "clock_bias_rate" in _states:
            state.clock_bias_rate = self.clock_bias_rate

        if "iono" in _states:
            for sat, iono in self.iono.items():
                state.iono[sat] = iono

        if "isb" in _states:
            state.isb = self.isb

        if "tropo_wet" in _states:
            state.tropo_wet = self.tropo_wet

        state.add_additional_info("states", _states)
        state.add_additional_info("clock_master", self.get_additional_info("clock_master"))
        state.add_additional_info("clock_slave", self.get_additional_info("clock_slave"))

        return state

    def _init_states(self, metadata, position, velocity, clock_bias, clock_bias_rate, sat_list):
        _states = ["position", "clock_bias"]  # mandatory states

        # position (with default to [0, 0, 0])
        self.position = np.array(position if position is not None else [0, 0, 0])
        self.cov_position = np.zeros((3, 3))

        # velocity is an additional state, that is estimated when set by the user
        self.velocity = None
        self.cov_velocity = None
        self.clock_bias_rate = None
        self.cov_clock_bias_rate = None
        if metadata is not None and metadata["VELOCITY_EST"] or velocity is not None:
            _states += ["velocity", "clock_bias_rate"]
            self.velocity = np.array(velocity if velocity is not None else [0, 0, 0])
            self.cov_velocity = np.zeros((3, 3))
            self.clock_bias_rate = dict()
            self.cov_clock_bias_rate = dict()

            if metadata is not None:
                _constellations = metadata["CONSTELLATIONS"]
            elif clock_bias_rate is not None:
                _constellations = list(clock_bias_rate.keys())
            else:
                _constellations = None
            if _constellations:
                for constellation in _constellations:
                    self.clock_bias_rate[constellation] = clock_bias_rate[constellation] if clock_bias_rate else 0.0
                    self.cov_clock_bias_rate[constellation] = 0.0

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

        # tropo wet delay (optional -> in case the user defined it)
        self.tropo_wet = None
        self.cov_tropo_wet = None
        if metadata is not None and metadata["TROPO"].estimate_tropo():
            self.tropo_wet = 0.0
            self.cov_tropo_wet = 0.0
            _states.append("tropo_wet")

        self.add_additional_info("states", _states)

    def __str__(self):
        _states = self.get_additional_info("states")
        _str = f"position = {self.position}, " \
               f"clock bias = {self.clock_bias}"
        if "velocity" in _states:
            _str += f", velocity = {self.velocity}"
        if "isb" in _states:
            _str += f", ISB = {self.isb}"
        if "iono" in _states:
            _str += f", iono = {self.iono}"
        if "tropo_wet" in _states:
            _str += f", tropo_wet = {self.tropo_wet}"
        if "clock_bias_rate" in _states:
            _str += f", clock_bias_rate = {self.clock_bias_rate}"

        return f'{type(self).__name__}[{str(self.epoch)}]({_str})'

    def __repr__(self):
        return str(self)

    def get_clock_bias(self, constellation, time_correction):
        """
        Computes the clock bias for the required constellation, applying the ISB correction
        in case the slave constellation is selected

        Args:
             constellation(str)
             time_correction(dict): dict with the GGTO data from the NavigationData header
        """
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
        estimate_ggto = config_dict.get("model", "estimate_ggto")

        ggto = 0.0
        if estimate_ggto is False:
            week, sow = self.epoch.gnss_time
            # compute GGTO from broadcast message
            ggto = compute_ggto(time_correction, week, sow)
            if constellation == self.get_additional_info("clock_slave"):
                if constellation == "GPS":
                    # in case the slave constellation is GPS, we need to fix the GGTO
                    ggto = -ggto
        # else: GGTO is estimated together with the ISB
        return self.isb - ggto

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
        if "geometry" in self._info:
            exportable_lst.append("satellite_azel")
        if "dop_ecef" in self._info:
            exportable_lst.append("dop_ecef")
        if "dop_local" in self._info:
            exportable_lst.append("dop_local")
        if "prefit_residuals" in self._info:
            exportable_lst.append("prefit_residuals")
        if "postfit_residuals" in self._info:
            exportable_lst.append("postfit_residuals")
        return exportable_lst
