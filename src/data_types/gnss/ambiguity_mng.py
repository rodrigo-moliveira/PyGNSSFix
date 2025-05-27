# TODO: add documentation
from src.io.config.enums import EnumLambdaMethod
from src.io.config.config import config_dict


class Ambiguity:
    def __init__(self, val, cov, fixed=False):
        self.val = val
        self.cov = cov
        self.fixed = fixed

    def clone(self):
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
        self.amb_resolution_model = EnumLambdaMethod.init_model(config_dict.get("solver", "ambiguity_resolution", "method"))
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



