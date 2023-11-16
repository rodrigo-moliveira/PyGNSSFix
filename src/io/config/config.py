"""Configuration handler

The configuration is a simple dictionary. See :ref:`configuration` for
details.
"""
from src.data_types.gnss.data_type import data_type_from_rinex
from src.errors import ConfigError
from src.io.config.enums import EnumPositioningMode


class Config(dict):
    """Configuration

    Example:
        .. code-block:: python

            from space.config import config

            print(config['env']['eop_missing_policy'])
            print(config.get('env', 'non-existent-field', fallback=25))

    """

    _instance = None

    def init(self, initial_dict):
        super().__init__()  # Call the parent class (dict) constructor
        self.update(initial_dict)  # Update the dictionary with the values from initial_dict

        # TODO: add config validation
        #   see https://towardsdatascience.com/how-to-use-json-schema-to-validate-json-documents-ae9d8d1db344

        # model initializations
        self["model"]["mode"] = EnumPositioningMode.init_model(self["model"]["mode"])

    def get(self, *keys, fallback=ConfigError):
        """Retrieve a value in the config, if the value is not available
        give the fallback value specified.
        """

        full_keys = list(keys).copy()
        key = None

        section, *keys = keys
        out = super().get(section, fallback)

        while isinstance(out, dict):
            key = keys.pop(0)
            out = out.get(key, fallback)

        if keys and out is not fallback:
            raise ConfigError(
                "Dict structure mismatch : Looked for '{}', stopped at '{}'".format(
                    ".".join(full_keys), key
                )
            )

        if out is ConfigError:
            raise ConfigError(
                "Invalid dict structure: Could not find {} in config".format(".".join(full_keys))
            )

        return out

    def set(self, *args):
        """Set a value in the config dictionary

        The last argument is the value to set

        Example:

        .. code-block:: python

            config.set('aaa', 'bbb', True)
            print(config)
            # {
            #     'aaa': {
            #         'bbb': True
            #     }
            # }
        """

        # split arguments in keys and value
        *first_keys, last_key, value = args

        subdict = self
        for k in first_keys:
            subdict.setdefault(k, {})
            subdict = subdict[k]

        subdict[last_key] = value

    def get_services(self):
        if "services" not in self:
            services = {}
            # iterate over all active constellations
            constellations = self.get("model", "constellations")
            for constellation in constellations:
                const_upper = constellation.upper()
                if const_upper == "GPS" or const_upper == "GAL":
                    services[const_upper] = self.get("model", const_upper, "observations")
                    self.set("services", services)
        return self["services"]

    def get_model(self):
        return self["model"]["mode"]

    def get_obs_std(self):
        if "obs_std" not in self:
            obs_dict = {}
            service_dict = self.get_services()
            for constellation, services in service_dict.items():
                obs_std_list = self.get("model", constellation, "pr_obs_std")
                obs_dict[constellation] = {}
                for index, service in enumerate(services):
                    datatype = data_type_from_rinex(f"C{service}", constellation)
                    obs_dict[constellation][datatype] = obs_std_list[index]
            self.set("obs_std", obs_dict)

        return self["obs_std"]

    def set_obs_std(self, constellation, datatype, std):
        # first build the obs_std dict
        self.get_obs_std()

        if constellation not in self["obs_std"]:
            self["obs_std"][constellation] = {}
        self["obs_std"][constellation][datatype] = std


config_dict = Config()
