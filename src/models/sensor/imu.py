""" Class to store the IMU Model (stochastic model stats) """

import os
from pathlib import Path
from typing import Any

from src.models.noise.noise_configuration import NoiseModel


__PATH__ = Path(__file__).parent.absolute()
SENSOR_PATH_MAP = {
    "low-end": os.path.abspath(__PATH__/"../../../inputs/sensors/low-end/imu.json"),
    "mid-end": os.path.abspath(__PATH__/"../../../inputs/sensors/mid-end/imu.json"),
    "high-end": os.path.abspath(__PATH__/"../../../inputs/sensors/high-end/imu.json")
}


class IMU:
    """
    Manage stochastic error models for an inertial measurement unit.

    The IMU contains separate stochastic models for the gyroscope and accelerometer. Each sensor can be configured
    with models for misalignment, scale factor, bias, bias drift, and observation noise.

    The available stochastic processes for each error model are:

        gyroscope:
            misalignment:
                - "gauss_markov"
                - "random_walk"
                - "random_constant"
                - "constant"

            scale_factor:
                - "gauss_markov"
                - "random_walk"
                - "random_constant"
                - "constant"

            bias_constant:
                - "constant"
                - "random_constant"

            bias_drift:
                - "gauss_markov"
                - "random_walk"

            observation_noise:
                - "white_noise"

        accelerometer:
            misalignment:
                - "gauss_markov"
                - "random_walk"
                - "random_constant"
                - "constant"

            scale_factor:
                - "gauss_markov"
                - "random_walk"
                - "random_constant"
                - "constant"

            bias_constant:
                - "constant"
                - "random_constant"

            bias_drift:
                - "gauss_markov"
                - "random_walk"

            observation_noise:
                - "white_noise"

    Continuous- and discrete-time process definitions
    --------------------------------------------------
    The stochastic processes "gauss_markov", "random_walk", and "white_noise" are defined using continuous-time
    statistical parameters. These parameters are discretized internally when the process realization is generated,
    using the specified sampling interval.

    In contrast, "constant" and "random_constant" processes are defined directly in discrete time and do not require
    time discretization.

    In particular:

        - Gauss-Markov:
            continuous-time PSD and correlation time are provided.

        - Random walk:
            continuous-time PSD is provided.

        - White noise:
            continuous-time PSD is provided.

        - Constant:
            value is directly specified in discrete time.

        - Random constant:
            standard deviation is directly specified in discrete time
            and a single random value is generated and maintained for
            the entire realization.

    Attributes:
        _data (dict):
            Dictionary containing the configured stochastic models for
            the gyroscope and accelerometer. Each sensor contains the
            following model entries:

                - "misalignment"
                - "scale_factor"
                - "bias_constant"
                - "bias_drift"
                - "observation_noise"
    """

    def __init__(self):
        """Initialize an IMU with empty stochastic error models."""
        # initialize data dict
        self._data: dict[str, dict[str, Any]] = {
            "gyroscope": {
                "misalignment": None,
                "scale_factor": None,
                "bias_constant": None,
                "bias_drift": None,
                "observation_noise": None,
            },
            "accelerometer": {
                "misalignment": None,
                "scale_factor": None,
                "bias_constant": None,
                "bias_drift": None,
                "observation_noise": None,
            },
        }

    def set(self, sensor, model, process, stats):
        """
        Configure a stochastic error model for an IMU sensor.

        Args:
            sensor (str): Sensor to configure. Must be either
                ``"gyroscope"`` or ``"accelerometer"``.

            model (str): Error model to configure. Must be one of:

                - ``"misalignment"``
                - ``"scale_factor"``
                - ``"bias_constant"``
                - ``"bias_drift"``
                - ``"observation_noise"``

            process (str): Stochastic process used to model the error.
                The available processes depend on the selected model:

                    misalignment:
                        ``"gauss_markov"``, ``"random_walk"``,
                        ``"random_constant"``, ``"constant"``

                    scale_factor:
                        ``"gauss_markov"``, ``"random_walk"``,
                        ``"random_constant"``, ``"constant"``

                    bias_constant:
                        ``"constant"``, ``"random_constant"``

                    bias_drift:
                        ``"gauss_markov"``, ``"random_walk"``

                    observation_noise:
                        ``"white_noise"``

            stats (dict): Statistical parameters required by the
                selected stochastic process.

                The following keys are supported:

                    ``process_noise``:
                        Process noise statistic. For continuous-time
                        processes, this represents the continuous-time
                        noise PSD (or the corresponding configured
                        representation). For discrete-time constant
                        processes, it represents the corresponding
                        discrete-time statistic.

                    ``tau``:
                        Correlation time for a Gauss-Markov process.
                        This parameter is not required by other process
                        types.

        Returns:
            None: This method updates the corresponding entry in the
            internal IMU model dictionary.

        Raises:
            KeyError: If ``sensor`` is not ``"gyroscope"`` or
                ``"accelerometer"``.
        """
        self._data[sensor][model] = NoiseModel(f"{sensor}_{model}", process,
                                               process_noise=stats.get("process_noise",None),
                                               correlation_time=stats.get("tau",None))


    def __getitem__(self, sensor):
        """
        Get the configured stochastic models for a sensor.

        Args:
            sensor (str): Sensor to retrieve. Must be either
                ``"gyroscope"`` or ``"accelerometer"``.

        Returns:
            dict: Dictionary containing the configured stochastic models
                for the requested sensor.

        Raises:
            KeyError: If ``sensor`` is not ``"gyroscope"`` or
                ``"accelerometer"``.
        """
        if sensor not in self._data:
            raise KeyError(f"Sensor must be 'gyroscope' or 'accelerometer', and not {sensor}")
        return self._data[sensor]

    def __str__(self):
        """
        Return a string representation of the IMU configuration.

        Returns:
            str: String representation containing the configured
                gyroscope and accelerometer models.
        """
        return f"IMU[{str(self._data)}]"
