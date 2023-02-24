import os
from src.errors import ConfigError
from src import PROJECT_PATH


def get_field(config_dict, field):
    try:
        _field = config_dict[field]
    except KeyError:
        raise ConfigError(f"Missing field `{field}` in configuration file.") from None

    return _field


class ConfigField(object):
    def __setattr__(self, *args):  # pragma: no cover
        raise TypeError("Cannot modify attributes of immutable object")

    def __delattr__(self, *args):  # pragma: no cover
        raise TypeError("Cannot modify attributes of immutable object")

    def super(self):
        return super()


class Log(ConfigField):
    available_levels = ["DEBUG", "INFO", "WARN", "ERROR", "FATAL"]
    __slots__ = ["log_level"]

    def __init__(self, **kwargs):
        # get & set fields
        super().super().__setattr__("log_level", get_field(kwargs, "minimum_level"))

        # input validation
        if self.log_level not in Log.available_levels:
            raise ConfigError(f"Invalid option in field log.minimum_level ({self.log_level}).\n\tAvailable values are "
                              f"{repr(Log.available_levels)}")


class Inputs(ConfigField):
    __slots__ = ["rinex_obs", "rinex_nav", "rinex_sp3", "rinex_clk", "snr_control", "first_epoch", "last_epoch"]

    def __init__(self, **kwargs):
        # get fields
        rinex_obs = get_field(kwargs, "rinex_obs_dir_path")
        rinex_nav = get_field(kwargs, "rinex_nav_dir_path")
        rinex_clk = get_field(kwargs, "rinex_clk_dir_path")
        rinex_sp3 = get_field(kwargs, "rinex_sp3_dir_path")

        snr_control = get_field(kwargs, "snr_control")

        arc = get_field(kwargs, "arc")

        first_epoch = get_field(arc, "first_epoch")
        last_epoch = get_field(arc, "last_epoch")

        # set fields
        super().super().__setattr__("rinex_obs", rinex_obs)
        super().super().__setattr__("rinex_nav", rinex_nav)
        super().super().__setattr__("rinex_clk", rinex_clk)
        super().super().__setattr__("rinex_sp3", rinex_sp3)
        super().super().__setattr__("snr_control", snr_control)
        super().super().__setattr__("first_epoch", first_epoch)
        super().super().__setattr__("last_epoch", last_epoch)

        # input validation
        print(os.path.exists(f"{PROJECT_PATH}\\{self.rinex_obs}"))
        for name, file in zip(["rinex obs", "rinex nav", "rinex sp3", "rinex clk"],
                              [self.rinex_obs, self.rinex_nav, self.rinex_sp3, self.rinex_clk]):
            if not os.path.exists(f"{PROJECT_PATH}\\{file}"):
                raise ConfigError(f"Input {name} directory path {PROJECT_PATH}\\{file} does not exist")


class ConfigGNSS(ConfigField):
    __slots__ = ["log", "inputs"]

    def __init__(self, **kwargs):
        # get sub-configs
        log_ = get_field(kwargs, "log")
        inputs_ = get_field(kwargs, "inputs")

        super().super().__setattr__("log", Log(**log_))
        super().super().__setattr__("inputs", Inputs(**inputs_))
