""" GNSS Sate Space Module
"""
import numpy as np

from src.data_mng import Container
from src import constants
from src.io.config import config_dict, EnumFrequencyModel
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

    Additional parameters:
        * epoch (Epoch): Epoch object for this time instant
        * _info (dict): auxiliary information (see paragraph below)
        * initial_state (GnssStateSpace): initial state space object with the initial states and covariances
            for the iteration process of the GNSS solver

    Some auxiliary information is saved in the `_info` attribute dict, namely:
        * _info["states"] provides a list with all valid states to be estimated, depending on the configuration
        * _info["clock_master"] provides the master constellation, when the ISB state is active
        * _info["clock_slave"] provides the slave constellation, when the ISB state is active

    Additionally, the GNSS solver will also save the following information in the `_info` dict parameter:
        * geometry
        * dop_ecef
        * dop_local
        * pr_prefit_residuals
        * pr_postfit_residuals
        * pr_rate_prefit_residuals
        * pr_rate_postfit_residuals
    """
    __states__ = ["position", "velocity", "clock_bias", "iono", "tropo_wet", "isb", "clock_bias_rate"]
    __covs__ = ["cov_position", "cov_velocity", "cov_clock_bias", "cov_iono", "cov_tropo_wet", "cov_isb",
                "cov_clock_bias_rate"]
    __slots__ = __states__ + __covs__ + ["epoch", "_info", "initial_state"]

    def __init__(self, metadata=None, epoch=None, sat_list=None):
        """
        GnssStateSpace constructor.

        Builds a GNSS state space object with the initial states and covariances defined in the `metadata` dict, defined
        in :py:class:`src.modules.gnss_solver.gnss_solver.GnssSolver`.
        If the metadata is not provided, all the state variables and covariances are initialized as None.

        Args:
            metadata(dict): dictionary with the metadata information on initial states and covariances
            epoch(src.data_types.date.Epoch): Epoch object for this time instant
            sat_list(list[src.data_types.gnss.Satellite]): list of available satellites (for iono initialization)
        """
        super().__init__()

        # dict to store state information
        self._info = dict()
        self.initial_state = None

        # initialize state variables
        if metadata is not None and sat_list is not None:
            self._init_states(metadata, sat_list)
        else:
            for state in self.__states__:
                setattr(self, state, None)
            for cov in self.__covs__:
                setattr(self, cov, None)

        # epoch
        self.epoch = epoch

    def clone(self):
        """ Returns a copy of this GnssStateSpace object. All state variables are deep-copied.

        Returns:
            GnssStateSpace : cloned GnssStateSpace instance
        """
        _states = self.get_additional_info("states")
        state = GnssStateSpace()
        state.epoch = self.epoch

        # deep copy position
        state.position = np.copy(self.position)
        state.cov_position = np.copy(self.cov_position)

        # deep copy velocity
        state.velocity = np.copy(self.velocity)
        state.cov_velocity = np.copy(self.cov_velocity)

        # deep copy clock bias
        state.clock_bias = self.clock_bias
        state.cov_clock_bias = self.cov_clock_bias

        # deep copy clock bias rate (if active)
        if "clock_bias_rate" in _states:
            state.clock_bias_rate = dict()
            state.cov_clock_bias_rate = dict()
            for const in self.clock_bias_rate.keys():
                state.clock_bias_rate[const] = self.clock_bias_rate[const]
                state.cov_clock_bias_rate[const] = self.cov_clock_bias_rate[const]

        if "iono" in _states:
            state.iono = dict()
            state.cov_iono = dict()
            for sat, iono in self.iono.items():
                state.iono[sat] = iono
                state.cov_iono[sat] = self.cov_iono[sat]

        if "isb" in _states:
            state.isb = self.isb
            state.cov_isb = self.cov_isb

        if "tropo_wet" in _states:
            state.tropo_wet = self.tropo_wet
            state.cov_tropo_wet = self.cov_tropo_wet

        state.add_additional_info("states", _states)
        state.add_additional_info("clock_master", self.get_additional_info("clock_master"))
        state.add_additional_info("clock_slave", self.get_additional_info("clock_slave"))

        return state

    def _init_states(self, metadata, sat_list):
        """ Initialize the state variables and covariances with the values provided in the metadata dict """
        _states = ["position", "clock_bias"]  # mandatory states

        # initialize position
        self.position = np.array(metadata["INITIAL_STATES"]["pos"][0:3], dtype=np.float64)
        self.cov_position = np.diag(np.array(metadata["INITIAL_STATES"]["pos"][3:6]))

        # initialize clock bias (convert input units from seconds to meters)
        self.clock_bias = list(metadata["INITIAL_STATES"].get("clock"))[0] * constants.SPEED_OF_LIGHT
        self.cov_clock_bias = list(metadata["INITIAL_STATES"].get("clock"))[1] * constants.SPEED_OF_LIGHT ** 2

        # velocity and clock bias rates are additional states, that are estimated when set by the user
        self.velocity = None
        self.cov_velocity = None
        self.clock_bias_rate = None
        self.cov_clock_bias_rate = None
        if metadata["VELOCITY_EST"]:
            _states += ["velocity", "clock_bias_rate"]
            self.velocity = np.array(metadata["INITIAL_STATES"]["vel"][0:3], dtype=np.float64)
            self.cov_velocity = np.diag(metadata["INITIAL_STATES"]["vel"][3:6])
            self.clock_bias_rate = dict()
            self.cov_clock_bias_rate = dict()

            for constellation in metadata["CONSTELLATIONS"]:
                self.clock_bias_rate[constellation] = list(metadata["INITIAL_STATES"].get("clock_rate"))[0] \
                                                      * constants.SPEED_OF_LIGHT
                self.cov_clock_bias_rate[constellation] = list(metadata["INITIAL_STATES"].get("clock_rate"))[1] \
                                                          * constants.SPEED_OF_LIGHT ** 2

        # iono dict (if dual-frequency mode is selected)
        self.iono = dict()
        self.cov_iono = dict()
        for constellation in metadata["CONSTELLATIONS"]:
            if metadata["MODEL"][constellation] == EnumFrequencyModel.DUAL_FREQ and \
                    metadata["IONO"][constellation].estimate_diono():
                for sat in sat_list:
                    if sat.sat_system == constellation:
                        self.iono[sat] = list(metadata["INITIAL_STATES"].get("iono"))[0]
                        self.cov_iono[sat] = list(metadata["INITIAL_STATES"].get("iono"))[1]
        if len(self.iono) >= 1:
            _states.append("iono")

        # isb (optional -> in case there are 2 constellations)
        self.isb = None
        self.cov_isb = None
        if len(metadata["CONSTELLATIONS"]) > 1:
            # convert input units from seconds to meters
            self.isb = list(metadata["INITIAL_STATES"].get("isb"))[0] * constants.SPEED_OF_LIGHT
            self.cov_isb = list(metadata["INITIAL_STATES"].get("isb"))[1] * constants.SPEED_OF_LIGHT ** 2
            self.add_additional_info("clock_master", metadata["CONSTELLATIONS"][0])
            self.add_additional_info("clock_slave", metadata["CONSTELLATIONS"][1])
            _states.append("isb")
        else:
            self.add_additional_info("clock_master", None)
            self.add_additional_info("clock_slave", None)

        # tropo wet delay (optional -> in case the user defined it)
        self.tropo_wet = None
        self.cov_tropo_wet = None
        if metadata["TROPO"].estimate_tropo():
            self.tropo_wet = list(metadata["INITIAL_STATES"].get("tropo"))[0]
            self.cov_tropo_wet = list(metadata["INITIAL_STATES"].get("tropo"))[1]
            _states.append("tropo_wet")

        self.add_additional_info("states", _states)

    def update_sat_list(self, sat_list):
        """ Update the satellite list of the internal iono state with the one provided as argument.

        Args:
            sat_list(list[src.data_types.gnss.Satellite]) : list of available satellites
        """
        _states = self.get_additional_info("states")
        if "iono" not in _states:
            return

        for sat in sat_list:
            if sat not in self.iono:
                self.iono[sat] = 0.0
                self.cov_iono[sat] = 1.0
                if self.initial_state is not None and sat in self.initial_state.iono:
                    self.iono[sat] = self.initial_state.iono[sat]
                    self.cov_iono[sat] = self.initial_state.cov_iono[sat]

        _to_remove = []
        for sat in self.iono:
            if sat not in sat_list:
                _to_remove.append(sat)  # satellite to be removed
        for sat in _to_remove:
            self.iono.pop(sat)
            self.cov_iono.pop(sat)

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
        in case the slave constellation is selected.

        Args:
             constellation(str): the clock bias is computed with respect to the system time defined by the
                `constellation` parameter
             time_correction(dict or None): dict with the GGTO data from the
                :py:class:`src.data_mng.gnss.navigation_data.NavigationHeader` instance
        Returns:
            float: clock bias for the required constellation in meters
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
        """
        Fetches the ISB (Inter-System bias) for the required constellation. Two alternatives are available:
            * the filter estimates the ISB together with the GGTO (GPS-to-GAL Time Offset), that is
                filter_ISB = (dr(B) - dr(A)) - t^AB
            * the filter applies the GGTO correction from the navigation data and only estimates the ISB, that is
                filter_ISB = (dr(B) - dr(A))

        In any case, the convertion from the clock bias of the master to the slave constellation is as follows:
            dtr(B) = dtr(A) + ISB(AB)
        where:
            ISB = (dr(B) - dr(A)) - t^AB
            dtr(A) and dtr(B) are the master and slave receiver clock biases
            dr(A) and dr(B) are the master and slave receiver hardware biases
            t^AB is the system time offset of B wrt A (t_B - t_A)

        For more information, see the discussion of Eqs. (21.60) and (21.61) of [1].

        Args:
             constellation(str): the clock bias is computed with respect to the system time defined by the
                `constellation` parameter
             time_correction(dict or None): dict with the GGTO data from the
                :py:class:`src.data_mng.gnss.navigation_data.NavigationHeader` instance
        Returns:
            float: ISB (Inter-System bias) for the required constellation in meters

        Reference:
            [1] Springer Handbook of Global Navigation Satellite Systems, Peter J.G. Teunissen, Oliver Montenbruck,
                Springer Cham, 2017
        """
        estimate_ggto = config_dict.get("model", "estimate_ggto")

        ggto = 0.0
        if estimate_ggto is False:
            week, sow = self.epoch.gnss_time
            # compute GGTO from broadcast message
            ggto = compute_ggto(time_correction, week, sow) * constants.SPEED_OF_LIGHT
            if constellation == self.get_additional_info("clock_slave"):
                if constellation == "GPS":
                    # in case the slave constellation is GPS, we need to fix the GGTO
                    ggto = -ggto
        # else: GGTO is estimated together with the ISB
        return self.isb - ggto

    def add_additional_info(self, arg: str, val):
        self._info[arg] = val

    def get_additional_info(self, arg):
        if arg in self._info:
            return self._info[arg]
        return None

    def get_exportable_lst(self):
        """
        Returns:
            list: list with the name of all state variables and additional data to be written in the output files
            See module :py:mod:`src.io.states.export_states`
        """
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
        if "pr_prefit_residuals" in self._info:
            exportable_lst.append("pr_prefit_residuals")
        if "pr_postfit_residuals" in self._info:
            exportable_lst.append("pr_postfit_residuals")
        if "pr_rate_prefit_residuals" in self._info:
            exportable_lst.append("pr_rate_prefit_residuals")
        if "pr_rate_postfit_residuals" in self._info:
            exportable_lst.append("pr_rate_postfit_residuals")
        return exportable_lst
