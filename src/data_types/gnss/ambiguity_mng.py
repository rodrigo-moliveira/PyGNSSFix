"""Module for managing GNSS ambiguities and performing ambiguity resolution using the LAMBDA algorithm."""
import numpy as np

from src.lambda_alg import LAMBDA
from src.io.config.enums import EnumLambdaMethod
from src.io.config.config import config_dict

# TODO: add documentation

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

    def reset(self):
        """ Resets the ambiguity to its initial state. """
        self.val = 0.0
        self.cov = 0.0
        self.fixed = False

    def clone(self):
        """ Creates a clone of the Ambiguity object. """
        obj = Ambiguity(self.val, self.cov, fixed=self.fixed)
        return obj

    def __str__(self):
        return f"Ambiguity(val={self.val}, cov={self.cov}, fixed={self.fixed})"

    def __repr__(self):
        return str(self)


class AmbiguityManager:
    def __init__(self, sat_list, init_val, init_cov, cp_types):
        self.ambiguities = dict()
        self.init_val = init_val
        self.init_cov = init_cov
        self.cp_types = cp_types

        # user configurations
        self.amb_resolution_enable = config_dict.get("solver", "ambiguity_resolution", "enabled")
        self.amb_resolution_model = EnumLambdaMethod.init_model(config_dict.get("solver",
                                                                                "ambiguity_resolution", "method"))
        self.amb_resolution_p0 = config_dict.get("solver", "ambiguity_resolution", "P0")
        self.amb_resolution_mu = config_dict.get("solver", "ambiguity_resolution", "mu")

        for sat in sat_list:
            # NOTE: when other types of ambiguity are implemented (NL or WL), this part should be modified
            self.ambiguities[sat] = dict()
            for cp_type in cp_types[sat.sat_system]:
                self.ambiguities[sat][cp_type] = Ambiguity(init_val, init_cov)

    def copy(self):
        return self

    def clone(self):
        obj = AmbiguityManager(
            sat_list=list(self.ambiguities.keys()),
            init_val=self.init_val,
            init_cov=self.init_cov,
            cp_types=self.cp_types
        )
        for sat in self.ambiguities:
            for cp_type in self.ambiguities[sat]:
                obj.ambiguities[sat][cp_type] = self.ambiguities[sat][cp_type].clone()
        return obj

    def add_ambiguity(self, sat, cp_type, other_ambiguity=None):
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
        return self.ambiguities.pop(key, default)

    def unfix_ambiguities(self, sat_system):
        """Unfix ambiguities for a given satellite system."""
        for sat in self.ambiguities:
            if sat.sat_system == sat_system:
                for cp_type in self.ambiguities[sat]:
                    self.ambiguities[sat][cp_type].reset()

    def main_fix(self, index_map, state, dX, cov):
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
            a_fixed, sq_norm, Ps, Qz_hat, Z, n_fixed, mu_out = LAMBDA.main(
                np.array(ambiguities), amb_cov, self.amb_resolution_model.value, 2, self.amb_resolution_p0,
                self.amb_resolution_mu)

            # Perform acceptance test...
            # TODO...

            # fix ambiguities in state space
            for sat, cp_types in index_map["ambiguity"].items():
                if pivot_dict[sat.sat_system] != sat:
                    for cp_type in cp_types:
                        idx_amb = cp_types[cp_type]
                        self[sat][cp_type].val = a_fixed[idx_amb - lowest_index, 0]
                        self[sat][cp_type].cov = amb_cov[idx_amb - lowest_index, idx_amb - lowest_index]
                        self[sat][cp_type].fixed = True

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
            print("Error in main_fix:", e)
            return None, None
            # TODO: properly take care of this
        return dX_updated, cov_updated
