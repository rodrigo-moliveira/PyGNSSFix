""" Module for importing IMU sensor data from the config json file. """
import json
from jsonschema import validate, ValidationError
import numpy as np

from src import PROJECT_PATH
from src.data_mng.unit_conversions import convert_unit
from src.errors import ConfigError

_SENSORS = ['gyroscope', 'accelerometer']
_ERROR_MODELS = ['misalignment', 'scale_factor', 'bias_constant', 'bias_drift', 'observation_noise']
_STOCHASTIC_PROCESSES = ["gauss_markov", "random_walk", "random_constant", "constant", "white_noise"]
_UNIT_MAP = {
    "gyroscope": {
        "misalignment": ("deg", "rad"),
        "scale_factor": ("ppm", ""),
        "bias_constant": ("deg/hr", "rad/s"),
        "bias_drift": ("deg/hr", "rad/s"),
        "observation_noise": ("deg/sqrt(hr)", "rad/sqrt(s)")
    },
    "accelerometer": {
        "misalignment": ("deg", "rad"),
        "scale_factor": ("ppm", ""),
        "bias_constant": ("mg", "m/s^2"),
        "bias_drift": ("mg", "m/s^2"),
        "observation_noise": ("m/s/sqrt(hr)", "m/s/sqrt(s)")
    },
}


def _get_process_stats(process, process_dict):
    """
    Extract the process statistics from the configuration dictionary.

    The statistics stored in the JSON configuration depend on the
    selected stochastic process. This function converts the
    process-specific configuration into the common dictionary format
    expected by :meth:`IMU.set`.
    """
    stats = {}
    if process == "gauss_markov":
        stats["process_noise"] = process_dict["Prn"]
        stats["tau"] = process_dict["Tau"]

    elif process == "random_walk":
        stats["process_noise"] = process_dict["Prn"]

    elif process == "random_constant" or process == "white_noise":
        stats["process_noise"] = process_dict["std"]

    elif process == "constant":
        stats["process_noise"] = process_dict["val"]

    return stats


def _convert_to_SI_units(sensor, model, stats):
    """
    Convert process statistics from configuration units to SI units.

    The input units are determined from the sensor and error model using
    ``_UNIT_MAP``. The ``process_noise`` values are converted using
    :func:`convert_unit` and the converted values replace the original
    values in ``stats``.

    This conversion is applied before the process statistics are passed
    to the IMU model.

    Examples of conversions include degrees to radians and values
    specified relative to standard gravity (g) to acceleration in
    meters per second squared.
    """
    # get units
    unit_in, unit_out = _UNIT_MAP[sensor][model]

    vec = np.array(stats["process_noise"])
    vec_out = convert_unit(vec, 3 * [unit_in], 3 * [unit_out])

    # update stats dict
    stats["process_noise"] = vec_out


def read_imu_file(path, imu, log=None):
    """
    Read and process an IMU configuration file.

    The configuration file is expected to be a JSON file conforming to
    the IMU JSON schema. For each sensor and error model, the selected
    stochastic process and its associated statistics are read from the
    configuration file, converted to SI units, and assigned to the
    provided :class:`IMU` instance.

    The processing sequence is:

        1. Read the JSON configuration file.
        2. Validate the configuration against the IMU JSON schema.
        3. Iterate over all supported sensors and error models.
        4. Determine the selected stochastic process.
        5. Extract the process statistics.
        6. Convert the statistics to SI units.
        7. Configure the corresponding IMU error model.

    Args:
        path (str or pathlib.Path): Path to the IMU JSON configuration
            file.

        imu (IMU): IMU instance that will be configured using the
            parameters read from the file.

        log (logging.Logger, optional): Logger used to report the
            processing of each sensor error model. If ``None``, no
            logging is performed.

    Returns:
        None: The provided ``imu`` object is configured in place.

    Raises:
        ConfigError: If the JSON configuration does not conform to the
            IMU JSON schema.

        FileNotFoundError: If the configuration file or the IMU schema
            file cannot be found.

        json.JSONDecodeError: If the configuration file or schema is
            not valid JSON.

        KeyError: If a required sensor, error model, process, or
            configuration parameter is missing.

    Notes:
        Process statistics are converted to SI units before being passed
        to :meth:`IMU.set`. The original JSON configuration is not
        modified.
    """
    # Opening JSON file
    with open(path) as json_file:
        json_dict = json.load(json_file)

        schema_path = PROJECT_PATH / "src/io/config/resources/imu_schema.json"
        with open(schema_path) as schema_file:
            schema = json.load(schema_file)

        try:
            validate(json_dict, schema)
        except ValidationError as e:
            raise ConfigError(f"Error validating config file: {e}")

        for sensor in _SENSORS:
            for model in _ERROR_MODELS:

                # get selected model
                process = json_dict[sensor][model]["select_model"]

                # get statistics for this model and process
                stats = _get_process_stats(process, json_dict[sensor][model][process])

                # convert data to SI units (degrees to rads, ref acceleration g to m/s^2)
                _convert_to_SI_units(sensor, model, stats)

                if log is not None:
                    log.info(f"processing sensor {sensor} for model {model} with process {process} with stats {stats} [SI units]")

                # set the imu for this model
                imu.set(sensor, model, process, stats)

