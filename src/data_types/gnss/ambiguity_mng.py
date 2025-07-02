"""Module for managing GNSS ambiguities and performing ambiguity resolution using the LAMBDA algorithm."""
import numpy as np

from src.lambda_alg import LAMBDA
from src.io.config.enums import EnumLambdaMethod, EnumSatelliteBias
from src.io.config.config import config_dict
from src.common_log import get_logger, GNSS_ALG_LOG


class Ambiguity:
    """Class representing a GNSS ambiguity with its value, covariance, and fixed status.

    Attributes:
        val (float): The value of the ambiguity (in cycles).
        cov (float): The covariance of the ambiguity (in cycles^2).
        fixed (bool): Indicates whether the ambiguity is fixed or not.
    """

    def __init__(self, val, cov, fixed=False):
        """ Constructor for the Ambiguity class.

        Args:
            val (float): The value of the ambiguity (in cycles).
            cov (float): The covariance of the ambiguity (in cycles^2).
            fixed (bool): Indicates whether the ambiguity is fixed or not. Defaults to False.
        """
        self.val = val
        self.cov = cov
        self.fixed = fixed

    def reset(self, cov):
        """ Resets the ambiguity to its initial state. """
        self.val = 0.0
        self.cov = cov
        self.fixed = False

    def clone(self):
        """ Creates a clone (deep copy) of the Ambiguity object. """
        obj = Ambiguity(self.val, self.cov, fixed=self.fixed)
        return obj

    def __str__(self):
        return f"Ambiguity(val={self.val}, cov={self.cov}, fixed={self.fixed})"

    def __repr__(self):
        return str(self)


class AmbiguityManager:
    """ Class for managing GNSS ambiguities across multiple satellites and carrier phases.
    Manages the ambiguity resolution algorithm using the LAMBDA algorithm.

    Attributes:
        ambiguities (dict): Dictionary mapping satellites to their ambiguities for different carrier phase types.
        init_val (float): Initial value for ambiguities (in cycles).
        init_cov (float): Initial covariance for ambiguities (in cycles^2).
        cp_types (dict): Dictionary mapping satellite systems to their carrier phase types.
        amb_resolution_enable (bool): Flag to enable ambiguity resolution.
        amb_resolution_model (EnumLambdaMethod): The model used for ambiguity resolution.
        config_P0_PAR (float): The minimum required success rate for method PAR.
        config_P0_RatioTest (float): The fixed failure rate for method ILS + Ratio Test.
    """

    def __init__(self, sat_list, init_val, init_cov, cp_types, bias_enum=None):
        """ Constructor for the AmbiguityManager class.

        Args:
            sat_list (list): List of satellites for which ambiguities are managed.
            init_val (float): Initial value for ambiguities (in cycles).
            init_cov (float): Initial covariance for ambiguities (in cycles^2).
            cp_types (dict): Dictionary mapping satellite systems to their carrier phase types.
            bias_enum (EnumSatelliteBias): Enumeration for satellite biases (e.g., broadcast, OSB, DCB).
        """
        self.ambiguities = dict()
        self._log = get_logger(GNSS_ALG_LOG)
        self.init_val = init_val
        self.init_cov = init_cov
        self.cp_types = cp_types

        # user configurations
        self.amb_resolution_enable = config_dict.get("solver", "ambiguity_resolution", "enabled")
        self.amb_resolution_model = EnumLambdaMethod.init_model(config_dict.get("solver",
                                                                                "ambiguity_resolution", "method"))
        self.config_P0_PAR = config_dict.get("solver", "ambiguity_resolution", "P0_PAR")
        self.config_P0_RatioTest = config_dict.get("solver", "ambiguity_resolution", "P0_RatioTest")

        # Check if OSB products are enabled (AR is not possible with DCBs)
        if bias_enum is not None:
            if bias_enum != EnumSatelliteBias.OSB:
                self._log.warning(f"Ambiguity resolution is not possible with satellite biases as {bias_enum}. "
                                  f"OSB products are required for ambiguity resolution.")
                self.amb_resolution_enable = False

        for sat in sat_list:
            # NOTE: when other types of ambiguity are implemented (NL or WL), this part should be modified
            self.ambiguities[sat] = dict()
            for cp_type in cp_types[sat.sat_system]:
                self.ambiguities[sat][cp_type] = Ambiguity(init_val, init_cov)

    def copy(self):
        """ Creates a shallow copy of the AmbiguityManager object.

        Returns:
            AmbiguityManager: A shallow copy of AmbiguityManager with the same ambiguities.
        """
        return self

    def clone(self):
        """ Creates a deep copy of the AmbiguityManager object.

        Returns:
            AmbiguityManager: A deep copy of the AmbiguityManager with cloned ambiguities.
        """
        obj = AmbiguityManager(
            sat_list=list(self.ambiguities.keys()),
            init_val=self.init_val,
            init_cov=self.init_cov,
            cp_types=self.cp_types
        )
        obj.amb_resolution_enable = self.amb_resolution_enable
        for sat in self.ambiguities:
            for cp_type in self.ambiguities[sat]:
                obj.ambiguities[sat][cp_type] = self.ambiguities[sat][cp_type].clone()
        return obj

    def add_ambiguity(self, sat, cp_type, other_ambiguity=None):
        """ Adds an ambiguity for a given satellite and carrier phase type.

        Args:
            sat (src.data_types.gnss.Satellite): The satellite for which the ambiguity is added.
            cp_type (src.data_types.gnss.DataType): The carrier phase type for the ambiguity.
            other_ambiguity (Ambiguity or None): An existing ambiguity to clone (optional). Defaults to None.
        """
        if sat not in self.ambiguities:
            self.ambiguities[sat] = dict()
        if cp_type not in self.ambiguities[sat]:
            if other_ambiguity is None:
                self.ambiguities[sat][cp_type] = Ambiguity(self.init_val, self.init_cov)
            else:
                self.ambiguities[sat][cp_type] = other_ambiguity.clone()

    def __str__(self):
        return f"AmbiguityManager(ambiguities={self.ambiguities})"

    def __repr__(self):
        return str(self)

    def __iter__(self):
        return iter(self.ambiguities.items())  # enables unpacking as (key, val)

    def __getitem__(self, key):
        return self.ambiguities[key]

    def __setitem__(self, key, value):
        self.ambiguities[key] = value

    def __delitem__(self, key):
        del self.ambiguities[key]

    def __contains__(self, key):
        return key in self.ambiguities

    def pop(self, key, default=None):
        """ Removes and returns an ambiguity for a given key. """
        return self.ambiguities.pop(key, default)

    def unfix_ambiguities(self, constellation: str):
        """ Unfix ambiguities for all available satellites in the given constellation."""
        for sat in self.ambiguities:
            if sat.sat_system == constellation:
                for cp_type in self.ambiguities[sat]:
                    self.ambiguities[sat][cp_type].reset(cov=self.init_cov)

    def main_fix(self, index_map, state, dX, cov):
        """
        Main function to perform ambiguity resolution using the LAMBDA algorithm.

        This function is divided into several steps:
            1. Identify ambiguities and their indices.
            2. Create reduced state vector and covariance matrix.
            3. Call the LAMBDA algorithm to resolve ambiguities.
            4. Perform the acceptance test to assess the quality of the resolved ambiguities.
            5. Update the state vector and covariance with resolved ambiguities.

        This function handles exceptions and logs errors if ambiguity resolution fails.

        Args:
            index_map (dict): Dictionary mapping satellites to their ambiguity indices.
            state (src.data_mng.gnss.state_space.GnssStateSpace): The current state containing additional information.
            dX (np.ndarray): The estimated GNSS state vector.
            cov (np.ndarray): The covariance matrix of the estimated state vector.

        Returns:
            tuple[np.ndarray, np.ndarray]: Updated state vector and covariance matrix after ambiguity resolution.

        If any error occurs during the process, the original state vector and covariance are returned without changes.
        """
        if not self.amb_resolution_enable:
            # If ambiguity resolution is disabled, return the original state vector and covariance
            return dX, cov

        try:
            pivot_dict = state.get_additional_info("pivot")
            lowest_index = -1
            highest_index = -1
            ambiguities = []

            # Find the lowest and highest index of ambiguities
            for sat, cp_types in index_map["ambiguity"].items():
                if pivot_dict[sat.sat_system] != sat:
                    for cp_type in cp_types:
                        idx_amb = cp_types[cp_type]
                        if idx_amb < lowest_index or lowest_index == -1:
                            lowest_index = idx_amb
                        if idx_amb > highest_index:
                            highest_index = idx_amb
                        ambiguities.append(self[sat][cp_type].val + dX[idx_amb])

            # Create reduced state vector and covariance
            b_indices = list(range(0, lowest_index)) + list(range(highest_index + 1, len(dX)))
            a_indices = [i for i in range(len(dX)) if i not in b_indices]
            amb_cov = cov[lowest_index:highest_index + 1, lowest_index:highest_index + 1]

            # Partition the state vector and covariance
            x_B = dX[b_indices]
            # x_A = dX[a_indices]
            P_B = cov[np.ix_(b_indices, b_indices)]
            P_A = cov[np.ix_(a_indices, a_indices)]
            P_AB = cov[np.ix_(b_indices, a_indices)]

            # check if there are ambiguities to resolve
            if len(ambiguities) == 0:
                return dX, cov

            # Perform Ambiguity Resolution (call to LAMBDA)
            p0 = None
            if self.amb_resolution_model == EnumLambdaMethod.PAR:
                p0 = self.config_P0_PAR
            elif self.amb_resolution_model == EnumLambdaMethod.ILS_RATIO_TEST:
                p0 = self.config_P0_RatioTest
            a_fixed, sq_norm, Ps, Qz_hat, Z, n_fixed, mu_out = LAMBDA.main(
                np.array(ambiguities), amb_cov, ncands=2, method=self.amb_resolution_model.value, P0=p0)

            if n_fixed == 0:
                self._log.warning("No ambiguities were fixed during the resolution process (LAMBDA AR failed).")
                return dX, cov
            else:
                pass

            # get fixed integers
            idx_fixed = np.where(a_fixed[:, 0] - a_fixed[:, 0].astype(int) == 0)
            idx_fixed = idx_fixed[0]
            n_fixed = len(idx_fixed)
            if n_fixed == 0:
                self._log.warning("No ambiguities were fixed during the resolution process (LAMBDA AR failed).")
                return dX, cov
            else:
                pass

            # Perform acceptance test...
            # TODO: consider adding an acceptance test here (later, after KF)

            # fix ambiguities in state space
            for sat, cp_types in index_map["ambiguity"].items():
                if pivot_dict[sat.sat_system] != sat:
                    for cp_type in cp_types:
                        idx_amb = cp_types[cp_type]
                        if idx_amb - lowest_index in idx_fixed:
                            self[sat][cp_type].val = a_fixed[idx_amb - lowest_index, 0]
                            self[sat][cp_type].cov = amb_cov[idx_amb - lowest_index, idx_amb - lowest_index]
                            self[sat][cp_type].fixed = True
                            self._log.debug(f"Fixed ambiguity for satellite {sat}, type {cp_type}: "
                                            f"{self[sat][cp_type].val} [cycles], "
                                            f"cov = {self[sat][cp_type].cov} [cycles^2]")

            # Update the state vector and covariance with fixed ambiguities
            P_A_inv = np.linalg.inv(P_A)
            x_B2 = x_B - P_AB @ P_A_inv @ (np.array(ambiguities) - a_fixed[:, 0])
            P_B2 = P_B - P_AB @ P_A_inv @ P_AB.T
            cov_updated = cov.copy()
            cov_updated[np.ix_(b_indices, b_indices)] = P_B2

            dX_updated = np.zeros(len(dX))
            dX_updated[b_indices] = x_B2
            dX_updated[a_indices] = a_fixed[:, 0]

        except (Exception, SystemExit) as e:
            log = get_logger(GNSS_ALG_LOG)
            log.warning(f"Error in ambiguity resolution: {e}. Not fixing ambiguities in this epoch.")
            return dX, cov
        return dX_updated, cov_updated
