""" Configuration handler module
The configuration is read from the json file and validated by the jsonschema tool.

The default jsonschema files are stored in the ./resources folder
"""
import json
from jsonschema import validate, ValidationError

from src import PROJECT_PATH
from src.data_types.gnss import data_type_from_rinex
from src.errors import ConfigError
from .enums import EnumPositioningMode
from src.common_log import get_logger, IO_LOG

__all__ = ["config_dict"]


class Config(dict):
    """
    Configuration class that inherits from :py:class:`dict`.

    The configuration handler is initialzed from the configuration xml file. Several xml files are available
        (with different content) depending on the selected algorithm.
        The available algorithms are:
            * "gnss": GNSS PNT solver algorithm
            * "post_processing": post-processing / performance evaluation algorithm

        The :py:mod:`jsonschema` module is used for validating the xml files.

    Once the handler is initialized it can be imported everywhere and used as a global variable.
        The initialization and usage of the handler is examplified in the following example below:
        Examples:
            >>> from src.io.config import config_dict
            >>> config_filename = "path/to/config/file.xml"
            >>> try:
            >>>     with open(config_filename) as json_file:
            >>>         data = json.load(json_file)
            >>>         config_dict.init(data, alg="GNSS")
            >>>     except Exception as e:
            >>>         print(f"Error Reading Configuration File.")
    """
    _instance = None

    def __init__(self):
        super().__init__()  # Call the parent class (dict) constructor
        self.alg = None

    def init(self, initial_dict, alg="gnss"):
        """
        Initializes the configuration handler instance.

        Args:
            initial_dict(dict): a dict instance with the content of the loaded xml file. See :py:meth:`json.load`.
            alg(str): the algorithm of the run.

        Raises:
            ConfigError: when the initialization of the configuration handler fails, a ConfigError is raised
        """

        if alg.lower() not in ("gnss", "post_processing"):
            raise ConfigError(f"illegal value for {alg} argument. Available values are 'gnss', 'post_processing'")
        self.alg = alg.lower()

        # validate config file
        self._validate(initial_dict, alg)

        self.update(initial_dict)  # Update the dictionary with the values from initial_dict

        # model initializations
        if alg.lower() == "gnss":
            self["model"]["mode"] = EnumPositioningMode.init_model(self["model"]["mode"])

    def _validate(self, initial_dict, alg):
        # Read the schema from the file
        if alg.lower() == "gnss":
            schema_path = PROJECT_PATH / "src/io/config/resources/gnss_schema.json"
        elif alg.lower() == "post_processing":
            schema_path = PROJECT_PATH / "src/io/config/resources/performance_schema.json"
        else:
            raise ConfigError(f"invalid value for argument alg: {alg}")

        with open(schema_path) as schema_file:
            schema = json.load(schema_file)

        try:
            validate(initial_dict, schema)
        except ValidationError as e:
            raise ConfigError(f"Error validating config file: {e}")

    def get(self, *keys, fallback=ConfigError):
        """
        Retrieve a value in the config, if the value is not available,
        give the fallback value specified.

        Args:
            keys: list of cascated keys for the request configuration field
            fallback: fallback in case the list of keys is not found.
        Raises:
            ConfigError: if the fallback is the `ConfigError` exception, it is raised

        Examples:
            >>> from src.io.config import config_dict
            >>> config_dict.get("a", "b", "c", fallback=10)
            This will attempt to return the `field` {"a":{"b":{"c":field}}}, with a fallback of 10
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
        """ Set a value in the config dictionary

        The last argument is the value to set.

        Examples:
            >>> from src.io.config import config_dict
            >>> config_dict.set('a', 'b', True)
            will set the field {"a":{"b":True}}
        """

        # split arguments in keys and value
        *first_keys, last_key, value = args

        subdict = self
        for k in first_keys:
            subdict.setdefault(k, {})
            subdict = subdict[k]

        subdict[last_key] = value

    def get_services(self) -> dict:
        """
        Utility function for the config object initialized with the gnss algorithm.

        Returns:
            map: returns a map with the services (observation attributes) configured for each
                    constellation, as defined in the configuration field 'observations'.
        Raises:
            ConfigError: if this function is called with for an algorithm other than gnss
        """
        if self.alg != "gnss":
            raise ConfigError("The `get_services` function is only available for Config instances initialized"
                              " with the gnss algorithm.")
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

    def get_obs_std(self) -> dict:
        """
        Utility function for the config object initialized with the gnss algorithm.

        Returns:
            dict: a dict with the observation noise (standard deviation) for each observation.
        Raises:
            ConfigError: if this function is called with for an algorithm other than gnss. The exception is also
                raised if there are inconsistencies between number of signals and correspondent observation noise fields
        """
        if self.alg != "gnss":
            raise ConfigError("The `get_obs_std` function is only available for Config instances initialized"
                              " with the gnss algorithm.")

        if "obs_std" not in self:
            log = get_logger(IO_LOG)
            obs_dict = {}
            service_dict = self.get_services()

            for constellation, services in service_dict.items():
                obs_dict[constellation] = {}

                for obs_std_list, datatype_char in zip(("pr_obs_std", "doppler_obs_std"), ("C", "D")):
                    std_list = self.get("model", constellation, obs_std_list)
                    if len(std_list) < len(services):
                        raise ConfigError(f"Inconsistency between number of signals for {constellation}: {services} and"
                                          f" the correspondent observation noise in field {std_list}. "
                                          f"Please fix the configuration.")

                    for index, service in enumerate(services):
                        datatype = data_type_from_rinex(f"{datatype_char}{service}", constellation)
                        obs_dict[constellation][datatype] = std_list[index]
                        log.info(f"Setting observation noise std for datatype {datatype} and constellation "
                                 f"{constellation} equal to {obs_dict[constellation][datatype]}")
            self.set("obs_std", obs_dict)

        return self["obs_std"]

    def update_obs_std(self, constellation, datatype, std):
        """
        Utility function for the config object initialized with the gnss algorithm.

        This function updates the field variable `obs_std` for the provided datatype.
        Useful when updating the std value for the iono free data types

        Args:
            constellation(str): the constellation of the observation
            datatype(src.data_types.gnss.data_type.DataType): the `DataType` instance of the observation
            std(float): the std value to be updated
        Raises:
            ConfigError: if this function is called with for an algorithm other than gnss. The exception is also
                raised if there are inconsistencies between number of signals and correspondent observation noise fields
        """
        if self.alg != "gnss":
            raise ConfigError("The `update_obs_std` function is only available for Config instances initialized"
                              " with the gnss algorithm.")

        # first build the obs_std dict (if not built before)
        self.get_obs_std()

        if constellation not in self["obs_std"]:
            self["obs_std"][constellation] = {}
        self["obs_std"][constellation][datatype] = std


config_dict = Config()
