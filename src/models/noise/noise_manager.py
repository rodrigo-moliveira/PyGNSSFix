""" Module to manage the use of stochastic processes. """
from src.data_mng import Container
from src.errors import ConfigError
from src.models.noise.noise_configuration import NoiseModel


class GNSSNoiseManager(Container):
    """
    Manages noise models for GNSS Kalman Filter states.

    This class inherits from :py:class:`src.data_mng.Container` and holds one noise model per state.
    Each noise model is configured from a dictionary and supports types like white noise,
    random walk, and Gauss-Markov processes.

    Inherits from:
        Container: A base container class that stores named attributes.

    Attributes:
        position (NoiseModel): Noise model for receiver position (typically random walk).
        clock_bias (NoiseModel): Noise model for receiver clock bias (typically random walk).
        isb (NoiseModel): Inter-system bias model.
        iono (NoiseModel): Ionospheric delay model (e.g., Gauss-Markov or random walk).
        tropo (NoiseModel): Tropospheric delay model (e.g., Gauss-Markov or random walk).
        ambiguity (NoiseModel): Carrier-phase ambiguity model (typically random walk).
        phase_bias (NoiseModel): Phase bias model (usually fixed or random walk).
    """

    __slots__ = ["position", "clock_bias", "isb", "iono", "tropo", "ambiguity", "phase_bias"]

    def __init__(self, config_dict):
        """
        Initialize all GNSS state noise models from configuration.

        For each state in the GNSS filter, reads its noise model type and corresponding
        parameters from the provided configuration dictionary and instantiates the
        appropriate `NoiseModel` object.

        Args:
            config_dict (src.io.config.config.Config):  A nested configuration dictionary that provides noise model
                parameters for each state, including the type (white_noise, random_walk, gauss_markov),
                process noise standard deviation, and optionally the correlation time (for Gauss-Markov).

        Raises:
            ConfigError: If the noise model type is unknown or if required parameters are missing.
        """
        super().__init__()
        self.position = None
        self.clock_bias = None
        self.isb = None
        self.tropo = None
        self.iono = None
        self.ambiguity = None
        self.phase_bias = None

        for state in GNSSNoiseManager.__slots__:
            type_str = config_dict.get("solver", "noise_model", state, "type")

            relative_re_param = config_dict.get("solver", "noise_model", state, "relative_re-parameterization")

            if type_str == "white_noise" or type_str == "random_walk":
                process_noise = config_dict.get("solver", "noise_model", state, "white_noise_random_walk",
                                                "process_noise_std")
                noise_model = NoiseModel(state, type_str, process_noise=process_noise)
            elif type_str == "gauss_markov":
                process_noise = config_dict.get("solver", "noise_model", state, "gauss_markov",
                                                "process_noise_std")
                correlation_time = config_dict.get("solver", "noise_model", state, "gauss_markov",
                                                   "correlation_time")
                noise_model = NoiseModel(state, type_str, process_noise=process_noise,
                                         correlation_time=correlation_time)
            else:
                raise ConfigError(f"Unknown process noise type {type_str} for state {state}."
                                  f"Available options are white_noise, random_walk or gauss_markov.")

            if noise_model:
                noise_model.relative_re_param = relative_re_param
                self.__setattr__(state, noise_model)
            else:
                raise ConfigError(f"Configuration file has missing information for process noise {state}. "
                                  f"Please review it.")
