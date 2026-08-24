# Class to store the IMU Model (stochastic model stats)

import os
from pathlib import Path
from typing import Any

from src.models.noise.noise_configuration import NoiseModel

"""
Available stochastic models:
    - White Noise
    - Gauss Markov (may resort to Random Walk)
    - Random Constant
    - Constant

Definir o sensors model:
    - scale factor
        * constant or random constant or Gauss Markov
    - misalignment
        * constant or random constant or Gauss Markov
    - observation noise
        * White Noise
    - bias constant + bias stability drift
        * (Random constant or Constant) + ( Gauss Markov or Random Walk)


"""

__PATH__ = Path(__file__).parent.absolute()
SENSOR_PATH_MAP = {
    "low-end": os.path.abspath(__PATH__/"../../../inputs/sensors/low-end/imu.json"),
    "mid-end": os.path.abspath(__PATH__/"../../../inputs/sensors/mid-end/imu.json"),
    "high-end": os.path.abspath(__PATH__/"../../../inputs/sensors/high-end/imu.json")
}


class IMU:

    def __init__(self):

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
        self._data[sensor][model] = NoiseModel(f"{sensor}_{model}", process,
                                               process_noise=stats.get("process_noise",None),
                                               correlation_time=stats.get("process_noise",None))


    def __getitem__(self, sensor):
        if sensor not in self._data:
            raise KeyError(f"Sensor must be 'gyroscope' or 'accelerometer', and not {sensor}")
        return self._data[sensor]

    def __str__(self):
        return f"IMU[{str(self._data)}]"
