""" GNSS Sate Space Module
"""
import numpy as np

from src.data_mng import Container
from src import constants
from src.data_types.gnss.ambiguity_mng import AmbiguityManager
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
        * ambiguity (dict): two-level dictionary with satellites as first keys and CP types as second keys
            and estimated ambiguity values as values
        * phase_bias (dict): two-level dictionary with constellations as first keys and CP types as second keys
            and estimated phase bias values as values

    Additional parameters:
        * epoch (Epoch): Epoch object for this time instant
        * _info (dict): auxiliary information (see paragraph below)
        * initial_state (GnssStateSpace): initial state space object with the initial states and covariances
            for the iteration process of the GNSS solver
        * index_map (dict): dictionary with the index of each state variable in the state vector. Map to be used
            in the Normal Equations of the GNSS solver
            (see :py:class:`src.modules.gnss.solver.lsq_engine.LSQ_Engine_Position`)

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
        * sat_list
        * pivot: satellites used as pivot for the ambiguity estimation (one per constellation).
    """
    __states__ = ["position", "velocity", "clock_bias", "iono", "tropo_wet", "isb", "clock_bias_rate", "ambiguity",
                  "phase_bias"]
    __covs__ = ["cov_position", "cov_velocity", "cov_clock_bias", "cov_iono", "cov_tropo_wet", "cov_isb",
                "cov_clock_bias_rate", "cov_phase_bias"]
    __slots__ = __states__ + __covs__ + ["epoch", "_info", "initial_state", "index_map"]

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
        self.index_map = dict()

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

        state.add_additional_info("estimate_iono", self.get_additional_info("estimate_iono"))
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

        if "ambiguity" in _states:
            state.ambiguity = self.ambiguity.copy()

        if "phase_bias" in _states:
            state.phase_bias = dict()
            state.cov_phase_bias = dict()
            for const, cp_types in self.phase_bias.items():
                if const not in state.phase_bias:
                    state.phase_bias[const] = dict()
                    state.cov_phase_bias[const] = dict()
                for cp_type in cp_types:
                    state.phase_bias[const][cp_type] = self.phase_bias[const][cp_type]
                    state.cov_phase_bias[const][cp_type] = self.cov_phase_bias[const][cp_type]

        state.add_additional_info("states", _states)
        state.add_additional_info("clock_master", self.get_additional_info("clock_master"))
        state.add_additional_info("clock_slave", self.get_additional_info("clock_slave"))
        state.add_additional_info("pivot", self.get_additional_info("pivot"))
        state.add_additional_info("sat_list", self.get_additional_info("sat_list"))

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
        estimate_iono = set()
        for constellation in metadata["CONSTELLATIONS"]:
            if metadata["MODEL"][constellation] == EnumFrequencyModel.DUAL_FREQ and \
                    metadata["IONO"][constellation].estimate_diono():
                estimate_iono.add(constellation)
                for sat in sat_list:
                    if sat.sat_system == constellation:
                        self.iono[sat] = list(metadata["INITIAL_STATES"].get("iono"))[0]
                        self.cov_iono[sat] = list(metadata["INITIAL_STATES"].get("iono"))[1]
        if len(self.iono) >= 1:
            _states.append("iono")
        self.add_additional_info("estimate_iono", estimate_iono)

        # isb (optional -> in case there are 2 constellations)
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

        # ambiguity (optional -> in case algorithm is CP-based)
        if metadata["CP_BASED"]:
            cp_types = metadata["PHASES"]
            self.ambiguity = AmbiguityManager(sat_list, metadata["INITIAL_STATES"]["ambiguity"][0],
                                              metadata["INITIAL_STATES"]["ambiguity"][1], cp_types)
            _states.append("ambiguity")

            self.phase_bias = dict()
            self.cov_phase_bias = dict()
            for const, types in cp_types.items():
                if const not in self.phase_bias:
                    self.phase_bias[const] = dict()
                    self.cov_phase_bias[const] = dict()
                for cp_type in types:
                    self.phase_bias[const][cp_type] = metadata["INITIAL_STATES"]["phase_bias"][0] * \
                                                      constants.SPEED_OF_LIGHT
                    self.cov_phase_bias[const][cp_type] = metadata["INITIAL_STATES"]["phase_bias"][1] * \
                                                          constants.SPEED_OF_LIGHT ** 2
            _states.append("phase_bias")

        # Define Pivot Satellite
        if metadata["CP_BASED"]:
            pivot_dict = dict.fromkeys(metadata["CONSTELLATIONS"], None)

            for sat in sat_list:
                if pivot_dict[sat.sat_system] is None:
                    pivot_dict[sat.sat_system] = sat
        else:
            pivot_dict = None

        self.add_additional_info("pivot", pivot_dict)
        self.add_additional_info("states", _states)
        self.add_additional_info("sat_list", sat_list)

    def _update_sat_list(self, new_sat_list):
        """
        Update the internal satellite list (used for the iono and ambiguity states) with the one provided as argument.

        Args:
            new_sat_list(list[src.data_types.gnss.Satellite]) : list of available satellites
        """
        prev_sat_list = self.get_additional_info("sat_list")

        # add the new satellites to the internal data (new visible satellites)
        for sat in new_sat_list:
            if sat not in prev_sat_list:
                self._add_sat(sat)

        # remove satellites that are not in the new list (satellites that are not visible anymore)
        _to_remove = []
        for sat in prev_sat_list:
            if sat not in new_sat_list:
                _to_remove.append(sat)  # satellite to be removed
        for sat in _to_remove:
            self._remove_sat(sat)

        self.add_additional_info("sat_list", new_sat_list)

        # update pivot information
        pivot_dict = self.get_additional_info("pivot")
        for pivot in pivot_dict.values():
            if pivot not in new_sat_list:
                for sat in new_sat_list:
                    if sat.sat_system == pivot.sat_system:
                        pivot_dict[pivot.sat_system] = sat
                        break
        self.add_additional_info("pivot", pivot_dict)

    def _add_sat(self, sat):
        """ Add a new satellite to the state space. """
        _states = self.get_additional_info("states")
        if "iono" in _states:
            self.iono[sat] = 0.0
            self.cov_iono[sat] = 1.0

        if "ambiguity" in _states:
            cp_types = self.phase_bias[sat.sat_system].keys()
            for cp_type in cp_types:
                self.ambiguity.add_ambiguity(sat, cp_type)

        # initialize from initial_state
        if self.initial_state is not None:
            if "iono" in _states:
                if sat in self.initial_state.iono:
                    self.iono[sat] = self.initial_state.iono[sat]
                    self.cov_iono[sat] = self.initial_state.cov_iono[sat]

                else:
                    self.initial_state.iono[sat] = 0.0
                    self.initial_state.cov_iono[sat] = 1.0

            if "ambiguity" in _states:
                cp_types = self.phase_bias[sat.sat_system].keys()
                if sat in self.initial_state.ambiguity:
                    for cp_type in cp_types:
                        self.ambiguity.add_ambiguity(sat, cp_type, self.initial_state.ambiguity[sat][cp_type])
                else:
                    for cp_type in cp_types:
                        self.initial_state.ambiguity.add_ambiguity(sat, cp_type)

    def _remove_sat(self, sat):
        """ Remove a satellite from the state space. """
        _states = self.get_additional_info("states")
        if "iono" in _states:
            if sat in self.iono:
                self.iono.pop(sat)
                self.cov_iono.pop(sat)
        if "ambiguity" in _states:
            if sat in self.ambiguity:
                self.ambiguity.pop(sat)

    def build_index_map(self, sat_list):
        """
        Builds the index map for the state vector. The index map is a dictionary with:
            * keys - state names
            * values - index of the state in the state vector

        Available states:
            * position (dimension 3)
            * clock bias (dimension 1)
            * ISB (dimension 1)
            * iono (dimension 1 per satellite)
            * tropo wet (dimension 1)
            * ambiguity (dimension 1 per satellite and CP type, excluding the pivot satellite)
            * phase bias (dimension 1 per constellation and CP type)

        Args:
            sat_list(list[src.data_types.gnss.Satellite]) : list of available satellites
        """
        self._update_sat_list(sat_list)

        index_map = {}
        state_counter = 0
        states = self.get_additional_info("states")

        # position index
        if "position" in states:
            index_map["position"] = state_counter
            state_counter += 3

        # clock index
        if "clock_bias" in states:
            index_map["clock_bias"] = state_counter
            state_counter += 1

        # check if ISB is present
        if "isb" in states:
            index_map["isb"] = state_counter
            state_counter += 1

        # add iono states
        if "iono" in states:
            index_map["iono"] = dict()
            for sat in self.iono.keys():
                index_map["iono"][sat] = state_counter
                state_counter += 1

        if "tropo_wet" in states:
            index_map["tropo_wet"] = state_counter
            state_counter += 1

        if "ambiguity" in states:
            pivot_dict = self.get_additional_info("pivot")
            index_map["ambiguity"] = dict()
            for sat, cp_types in self.ambiguity:
                if pivot_dict[sat.sat_system] != sat:
                    index_map["ambiguity"][sat] = dict()
                    for cp_type, ambiguity in cp_types.items():
                        index_map["ambiguity"][sat][cp_type] = state_counter
                        state_counter += 1

        if "phase_bias" in states:
            index_map["phase_bias"] = dict()
            for const, cp_types in self.phase_bias.items():
                index_map["phase_bias"][const] = dict()
                for cp_type in cp_types.keys():
                    index_map["phase_bias"][const][cp_type] = state_counter
                    state_counter += 1

        index_map["total_states"] = state_counter

        self.index_map = index_map

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
        if "ambiguity" in _states:
            _str += f", ambiguity = {self.ambiguity}"
        if "phase_bias" in _states:
            _str += f", phase_bias = {self.phase_bias}"

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
        if "pos_prefit_residuals" in self._info:
            exportable_lst.append("pos_prefit_residuals")
        if "pos_postfit_residuals" in self._info:
            exportable_lst.append("pos_postfit_residuals")
        if "vel_prefit_residuals" in self._info:
            exportable_lst.append("vel_prefit_residuals")
        if "vel_postfit_residuals" in self._info:
            exportable_lst.append("vel_postfit_residuals")
        return exportable_lst
