import numpy as np

class NoiseProcess:
    def __init__(self):
        self._name = "General Process"

    def __str__(self):
        return self._name

    def compute(self, shape, *args):
        return np.zeros(shape)


class RandomWalkProcess(NoiseProcess):
    pass


class GaussMarkovProcess(NoiseProcess):
    pass


class WhiteNoiseProcess(NoiseProcess):
    pass


class RandomConstantProcess(NoiseProcess):
    pass
