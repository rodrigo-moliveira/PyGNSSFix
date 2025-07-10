from src.data_mng import Container
from src.errors import ConfigError
from src.models.noise.noise_configuration import NoiseModel


class GNSSNoiseManager(Container):
    # TODO: missing docstrings
    __slots__ = ["position", "clock_bias", "isb", "iono", "tropo", "ambiguity", "phase_bias"]

    def __init__(self, config_dict):
        super().__init__()

        for state in GNSSNoiseManager.__slots__:
            type_str = config_dict.get("solver", "noise_model", state, "type")
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
                self.__setattr__(state, noise_model)
            else:
                raise ConfigError(f"Configuration file has missing information for process noise {state}. "
                                  f"Please review it.")
