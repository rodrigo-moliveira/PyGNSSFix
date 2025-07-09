from src.data_mng import Container
from src.errors import ConfigError
from src.models.noise.noise_configuration import NoiseModel


class GNSSNoiseManager(Container):
    # TODO: missing docstrings
    __slots__ = ["position", "velocity", "clock_bias", "clock_drift", "isb", "iono", "tropo", "ambiguity", "phase_bias"]

    def __init__(self, config_dict):
        super().__init__()

        for arg in GNSSNoiseManager.__slots__:
            model = config_dict.get(arg, None)
            if model:
                type_str = model.get("type", None)
                process_noise = model.get("process_noise", None)
                if type_str and process_noise:
                    self.__setattr__(arg, NoiseModel(arg, type_str, process_noise))
                else:
                    raise ConfigError(f"Configuration file has missing information for process noise {arg}: "
                                      f"missing fields 'type' or 'process_noise'. ")
            else:
                raise ConfigError(f"Configuration file does not contain the process noise {arg}.")
