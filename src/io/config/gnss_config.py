import os
from src.errors import ConfigError, ConfigTypeError
from src import WORKSPACE_PATH
from src.data_types.gnss.service_utils import AvailableConstellations, Services
from .enums import *


def get_field(config_dict, field, field_type, root: str = None):
    _root = f"{root}.{field}" if root is not None else f"{field}"

    # get the field from the dict
    try:
        _field = config_dict[field]
    except KeyError:
        raise ConfigError(f"Missing field `{_root}` in configuration file.") from None

    # check if it is of the specified type
    # dual behaviour -> field_type can be a single type or a tuple of types
    if type(field_type) == tuple:
        if type(_field) not in field_type:
            raise ConfigTypeError(
                f"Field {_root} ({_field}) is not of the specified type. It is of type {type(_field)} "
                f"but should be one of the following types {field_type}")
    else:
        if type(_field) != field_type:
            raise ConfigTypeError(
                f"Field {_root} ({_field}) is not of the specified type. It is of type {type(_field)} "
                f"but should be of type {field_type}")
    return _field


def get_enum(enum_class, field, root: str = None):
    try:
        enum = enum_class(field)
    except ValueError:
        raise ConfigTypeError(
            f"The selected value in field {root} is illegal. Value is {field}, but available selections are "
            f"{enum_class.show_options()}") from None

    return enum


class ConfigField(object):
    def __setattr__(self, *args):  # pragma: no cover
        raise TypeError("Cannot modify attributes of immutable object")

    def __delattr__(self, *args):  # pragma: no cover
        raise TypeError("Cannot modify attributes of immutable object")

    def super(self):
        return super()


class Log(ConfigField):
    available_levels = ["DEBUG", "INFO", "WARN", "ERROR", "FATAL"]
    _root = "log"
    __slots__ = ["log_level"]

    def __init__(self, **kwargs):
        # get & set fields
        super().super().__setattr__("log_level", get_field(kwargs, "minimum_level", str, Log._root))

        # input validation
        if self.log_level not in Log.available_levels:
            raise ConfigError(f"Invalid option in field {Log._root}.minimum_level ({self.log_level}).\n\tAvailable "
                              f"values are {repr(Log.available_levels)}")


class Inputs(ConfigField):
    __slots__ = ["rinex_obs", "rinex_nav", "rinex_sp3", "rinex_clk", "snr_control", "first_epoch", "last_epoch", "rate"]
    _root = "inputs"

    def __init__(self, **kwargs):
        # get fields
        rinex_obs = get_field(kwargs, "rinex_obs_files", list, Inputs._root)
        rinex_nav = get_field(kwargs, "rinex_nav_files", list, Inputs._root)
        rinex_clk = get_field(kwargs, "rinex_clk_files", list, Inputs._root)
        rinex_sp3 = get_field(kwargs, "rinex_sp3_files", list, Inputs._root)

        snr_control = get_field(kwargs, "snr_control", int, Inputs._root)

        arc = get_field(kwargs, "arc", dict, Inputs._root)

        first_epoch = get_field(arc, "first_epoch", str, Inputs._root + ".arc")
        last_epoch = get_field(arc, "last_epoch", str, Inputs._root + ".arc")

        rate = get_field(kwargs, "rate", int, Inputs._root)

        # set fields
        super().super().__setattr__("rinex_obs", rinex_obs)
        super().super().__setattr__("rinex_nav", rinex_nav)
        super().super().__setattr__("rinex_clk", rinex_clk)
        super().super().__setattr__("rinex_sp3", rinex_sp3)
        super().super().__setattr__("snr_control", snr_control)
        super().super().__setattr__("first_epoch", first_epoch)
        super().super().__setattr__("last_epoch", last_epoch)
        super().super().__setattr__("rate", rate)

        # input validation
        # check rinex paths
        for name, file_list in zip(["rinex obs", "rinex nav", "rinex sp3", "rinex clk"],
                                   [self.rinex_obs, self.rinex_nav, self.rinex_sp3, self.rinex_clk]):
            for file in file_list:
                if not os.path.exists(f"{WORKSPACE_PATH}/{file}"):
                    raise ConfigError(f"Input {name} directory path {WORKSPACE_PATH}/{file} does not exist")

        # check snr control
        if self.snr_control < 1 or self.snr_control > 9:
            raise ConfigError(f"Field {Inputs._root}.snr_control must be an integer between 1 and 9")


class Constellation(ConfigField):
    __slots__ = ["observations", "obs_std", "tropo", "iono", "relativistic_corr"]
    _root = "model"

    def __init__(self, name: str, **kwargs):
        # get fields
        root = Constellation._root + "." + name

        observations = get_field(kwargs, "observations", list, root)
        obs_std = get_field(kwargs, "obs_std", list, root)
        tropo = get_enum(EnumOnOff, get_field(kwargs, "troposphere", int, root), f"{root}.troposphere")
        iono = get_enum(EnumIono, get_field(kwargs, "ionosphere", int, root), f"{root}.ionosphere")
        relativistic_corr = get_enum(EnumOnOff, get_field(kwargs, "relativistic_corrections", int, root),
                                     f"{root}.relativistic_corrections")

        super().super().__setattr__("observations", observations)
        super().super().__setattr__("obs_std", obs_std)
        super().super().__setattr__("tropo", tropo)
        super().super().__setattr__("iono", iono)
        super().super().__setattr__("relativistic_corr", relativistic_corr)

        # input validation
        # observations and obs_std lists must have the same size
        if len(observations) != len(obs_std):
            raise ConfigError(f"Lists in fields {root}.observations and {root}.obs_std must have the same size")

        # check if obs_std are floats
        for obs in obs_std:
            if type(obs) != float and type(obs) != int:
                raise ConfigTypeError(f"Values in field {root}.obs_std must be floats or ints")

        # check if observations are valid
        for obs in observations:
            if obs not in Services[name]:
                raise ConfigError(f"Observation {obs} in field {root}.observations is not valid. Available"
                                  f" observations for constellation {name} are {Services[name]}")


class Model(ConfigField):
    __slots__ = ["constellations", "GPS", "GAL"]
    _root = "model"

    def __init__(self, **kwargs):
        # get fields
        constellations = get_field(kwargs, "constellations", list, Inputs._root)
        gps = get_field(kwargs, "GPS", dict, Inputs._root)
        gal = get_field(kwargs, "GAL", dict, Inputs._root)

        # set fields
        super().super().__setattr__("constellations", constellations)
        super().super().__setattr__("GPS", Constellation("GPS", **gps))
        super().super().__setattr__("GAL", Constellation("GAL", **gal))

        # input validation
        for constellation in constellations:
            if constellation not in AvailableConstellations:
                raise ConfigError(f"Constellation {constellation} in field {Model._root}.constellations is not "
                                  f"known. Available constellations are {AvailableConstellations}")


class SatelliteStatus(ConfigField):
    __slots__ = ["sv_ura", "sv_minimum_ura", "sv_health"]
    _root = "solver.satellite_status"

    def __init__(self, **kwargs):
        # get fields
        sv_ura = get_field(kwargs, "SV_URA", bool, SatelliteStatus._root)
        sv_minimum_ura = get_field(kwargs, "SV_minimum_URA", (int, float), SatelliteStatus._root)
        sv_health = get_field(kwargs, "SV_health", bool, SatelliteStatus._root)

        super().super().__setattr__("sv_ura", sv_ura)
        super().super().__setattr__("sv_minimum_ura", sv_minimum_ura)
        super().super().__setattr__("sv_health", sv_health)


class Solver(ConfigField):
    __slots__ = ["algorithm", "n_iterations", "stop_criteria", "snr_filter", "elevation_filter", "sat_status",
                 "trans_time_alg"]
    _root = "solver"

    def __init__(self, **kwargs):
        # get fields
        algorithm = get_enum(EnumSolver, get_field(kwargs, "algorithm", int, Solver._root), f"{Solver._root}.algorithm")
        n_iterations = get_field(kwargs, "n_iterations", int, Solver._root)
        stop_criteria = get_field(kwargs, "stop_criteria", (float, int), Solver._root)
        snr_filter = get_field(kwargs, "snr_filter", (bool, float, int), Solver._root)
        elevation_filter = get_field(kwargs, "elevation_filter", (bool, float, int), Solver._root)
        sat_status = get_field(kwargs, "satellite_status", dict, Solver._root)
        trans_time_alg = get_enum(EnumTransmissionTime, get_field(kwargs, "transmission_time_alg", int, Solver._root),
                                  f"{Solver._root}.transmission_time_alg")

        super().super().__setattr__("algorithm", algorithm)
        super().super().__setattr__("n_iterations", n_iterations)
        super().super().__setattr__("stop_criteria", stop_criteria)
        super().super().__setattr__("snr_filter", snr_filter)
        super().super().__setattr__("elevation_filter", elevation_filter)
        super().super().__setattr__("trans_time_alg", trans_time_alg)
        super().super().__setattr__("sat_status", SatelliteStatus(**sat_status))

        # input validation
        for field in ["n_iterations", "stop_criteria", "snr_filter", "elevation_filter"]:
            val = getattr(self, field)

            if val <= 0:
                raise ConfigTypeError(f"Number of iterations in field {Solver._root}.{field} must be "
                                      f"positive. Selected value is {val}")


class PerformanceEval(ConfigField):
    __slots__ = ["static", "true_static_position", "show_plots", "output_path"]
    _root = "performance_evaluation"

    def __init__(self, **kwargs):
        # get fields
        static = get_field(kwargs, "static", bool, PerformanceEval._root)
        coordinates = get_field(kwargs, "true_static_position", dict, PerformanceEval._root)
        x = get_field(coordinates, "x_ecef", (float, int), f"{PerformanceEval._root}.x_ecef")
        y = get_field(coordinates, "y_ecef", (float, int), f"{PerformanceEval._root}.y_ecef")
        z = get_field(coordinates, "z_ecef", (float, int), f"{PerformanceEval._root}.z_ecef")

        show_plots = get_field(kwargs, "show_plots", bool, PerformanceEval._root)
        output_path = get_field(kwargs, "output_path", str, PerformanceEval._root)

        super().super().__setattr__("static", static)
        super().super().__setattr__("true_static_position", [x, y, z])
        super().super().__setattr__("show_plots", show_plots)
        super().super().__setattr__("output_path", output_path)


class ConfigGNSS(ConfigField):
    __slots__ = ["log", "inputs", "model", "solver", "performance_evaluation"]

    def __init__(self, **kwargs):
        # get sub-configs
        log_ = get_field(kwargs, "log", dict)
        inputs_ = get_field(kwargs, "inputs", dict)
        model_ = get_field(kwargs, "model", dict)
        solver_ = get_field(kwargs, "solver", dict)
        performance_eval_ = get_field(kwargs, "performance_evaluation", dict)

        super().super().__setattr__("log", Log(**log_))
        super().super().__setattr__("inputs", Inputs(**inputs_))
        super().super().__setattr__("model", Model(**model_))
        super().super().__setattr__("solver", Solver(**solver_))
        super().super().__setattr__("performance_evaluation", PerformanceEval(**performance_eval_))

    def get_services(self):
        services = {}

        # iterate over all active constellations
        constellations = self.model.constellations
        for constellation in constellations:
            if constellation.upper() == "GPS":
                services["GPS"] = self.model.GPS.observations
            elif constellation.upper() == "GAL":
                services["GAL"] = self.model.GAL.observations

        return services
