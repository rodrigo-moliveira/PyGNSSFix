"""Configuration handler

The configuration is a simple dictionary. See :ref:`configuration` for
details.
"""

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
        self["model"]["mode"] = EnumPositioningMode.init_model(self["model"]["mode"])
        # TODO: add config validation

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
        services = {}

        # iterate over all active constellations
        constellations = self.get("model", "constellations")
        for constellation in constellations:
            if constellation.upper() == "GPS":
                services["GPS"] = self.get("model", "GPS", "observations")
            elif constellation.upper() == "GAL":
                services["GAL"] = self.get("model", "GAL", "observations")

        return services

    def get_model(self):
        return self["model"]["mode"]

    def is_iono_free(self):
        """
        Returns true if user selected combined observation model (compute iono free observations)
        """
        return self["model"]["mode"] == EnumPositioningMode.SPS_IF


config_dict = Config()
