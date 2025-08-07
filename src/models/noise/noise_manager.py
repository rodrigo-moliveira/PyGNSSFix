""" Module to manage the use of stochastic processes. """
from src.data_mng import Container
from src.errors import ConfigError
from src.io.config import EnumNoiseProcess
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
        velocity (NoiseModel): Noise model for receiver velocity (typically random walk).
        clock_bias (NoiseModel): Noise model for receiver clock bias (typically random walk).
        clock_drift (NoiseModel): Noise model for receiver clock bias rate (typically random walk).
        isb (NoiseModel): Inter-system bias model.
        iono (NoiseModel): Ionospheric delay model (e.g., Gauss-Markov or random walk).
        tropo (NoiseModel): Tropospheric delay model (e.g., Gauss-Markov or random walk).
        ambiguity (NoiseModel): Carrier-phase ambiguity model (typically random walk).
        phase_bias (NoiseModel): Phase bias model (usually fixed or random walk).
        b_pv_model(bool): boolean to enable the PV Model (position and velocity correlated)
        b_clock_model(bool): boolean to enable the Clock Model (clock bias and drift correlated)
    """

    __slots__ = ["position", "velocity", "clock_bias", "clock_drift", "isb", "iono", "tropo", "ambiguity",
                 "phase_bias", "b_pv_model", "b_clock_model"]

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
        self.velocity = None
        self.clock_bias = None
        self.clock_drift = None
        self.isb = None
        self.tropo = None
        self.iono = None
        self.ambiguity = None
        self.phase_bias = None

        self.b_pv_model = config_dict.get("solver", "noise_model", "PV_Model_Enabled")
        self.b_clock_model = config_dict.get("solver", "noise_model", "Clock_Model_Enabled")

        for state in GNSSNoiseManager.__slots__:
            if state in ["b_pv_model", "b_clock_model"]:
                continue
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

        # validate model
        # In PV Model the Position and Velocity must be Random Walks
        if self.b_pv_model and getattr(self, "velocity").process_enum != EnumNoiseProcess.RANDOM_WALK:
            raise ConfigError(f"When the PV Model is active, the velocity model must be set to Random Walk. "
                              f"Please revise the configurations (position model = "
                              f"{getattr(self, 'position').process_enum}). ")

        if self.b_pv_model and getattr(self, "position").process_enum != EnumNoiseProcess.RANDOM_WALK:
            raise ConfigError(f"When the PV Model is active, the position model must be set to Random Walk. "
                              f"Please revise the configurations (position model = "
                              f"{getattr(self, 'position').process_enum}). ")

        # In Clock Model the clock bias and clock drift must be Random Walks
        if self.b_clock_model and getattr(self, "clock_bias").process_enum != EnumNoiseProcess.RANDOM_WALK:
            raise ConfigError(f"When the Clock Model is active, the clock bias model must be set to Random Walk. "
                              f"Please revise the configurations (clock bias model = "
                              f"{getattr(self, 'clock_bias').process_enum}). ")

        if self.b_clock_model and getattr(self, "clock_drift").process_enum != EnumNoiseProcess.RANDOM_WALK:
            raise ConfigError(f"When the Clock Model is active, the clock drift model must be set to Random Walk. "
                              f"Please revise the configurations (clock drift model = "
                              f"{getattr(self, 'clock_drift').process_enum}). ")
