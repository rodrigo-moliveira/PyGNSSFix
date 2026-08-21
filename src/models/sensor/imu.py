# Class to store the IMU Model (stochastic model stats)

import os
from pathlib import Path
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


"""class _Process:
    def __init__(self):
        self.process = None
        self.stats = None

    def set(self, process, stats):
        self.process = process
        self.stats = stats

    def __str__(self):
        return f"Process({self.process}, stats={self.stats})"

    def __repr__(self):
        return str(self)

    def get_stochastic_process(self, dim):
        if self.process == "gauss_markov":
            std = self.stats["prn"]
            tau = self.stats["tau"]
            if tau is None or tau < 0:
                # resort to Random Walk
                return st.RandomWalk(dim=dim, std=std, axis=3)
            return st.GaussMarkov(dim=dim, std=std, correlation_time=tau, axis=3)

        elif self.process == "random_constant":
            std = self.stats["std"]
            return st.RandomConstant(dim=dim, std=std, axis=3)

        elif self.process == "constant":
            bias = self.stats["bias"]
            return st.RandomConstant(dim=dim, std=bias, axis=3)

        elif self.process == "white_noise":
            std = self.stats["std"]
            return st.WhiteNoise(dim=dim, std=std, axis=3)

        else:
            # this should never happen...
            raise ValueError(f"Unknown process {self.process}")
"""

class IMU:

    def __init__(self):

        # initialize data dict
        self._data = {
            "gyroscope": {
                "misalignment": None,
                "scale_factor": None,
                "bias_constant": None,
                "bias_drift": None,
                "observation_noise": None
            },
            "accelerometer": {
                "misalignment": None,
                "scale_factor": None,
                "bias_constant": None,
                "bias_drift": None,
                "observation_noise": None
            }
        }

        #if isinstance(accuracy, str):
        #    if accuracy in SENSOR_PATH_MAP.keys():
        #        path = SENSOR_PATH_MAP[accuracy]
        #    else:
        #        path = accuracy

            # read imu file (may raise exceptions)
            #read_imu_file(path, self)
        #else:
         #   raise AttributeError(f"Provided 'accuracy' argument is not valid. Must be a string 'low-end', 'mid-end', "
          #                       f" 'high-end', or, alternatively, a path to a valid imu json file")

    def set(self, sensor, model, process, stats):
        self._data[sensor][model] = NoiseModel()
        self._data[sensor][model].set(process, stats)

    def __getitem__(self, sensor):
        if sensor not in self._data:
            raise KeyError(f"Sensor must be 'gyroscope' or 'accelerometer', and not {sensor}")
        return self._data[sensor]

    def __str__(self):
        return f"IMU[{str(self._data)}]"
