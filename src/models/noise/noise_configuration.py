import numpy as np

from src.io.config import EnumNoiseProcess
from .noise_gen import *


class NoiseModel:
    # TODO: missing docstrings
    def __init__(self, state_str, process_str, process_noise=None, correlation_time=None):
        self.name = state_str
        self.process_enum = EnumNoiseProcess.init_model(process_str)

        if self.process_enum == EnumNoiseProcess.WHITE_NOISE:
            self.process_gen = WhiteNoiseProcess()
        elif self.process_enum == EnumNoiseProcess.RANDOM_CONSTANT:
            self.process_gen = RandomConstantProcess()
        elif self.process_enum == EnumNoiseProcess.RANDOM_WALK:
            self.process_gen = RandomWalkProcess()
        elif self.process_enum == EnumNoiseProcess.GAUSS_MARKOV:
            self.process_gen = GaussMarkovProcess()
        else:
            raise ValueError(f"Expected EnumNoiseProcess enum, got {type(self.process_enum).__name__}")

        # internal variables
        self.process_noise = process_noise  # units of std (RW/WN - [units/sqrt(s)] or GM[units])
        self.correlation_time = correlation_time  # in [s]

    def gen(self):
        pass

    def get_stm_entry(self, time_step):
        if self.process_enum == EnumNoiseProcess.WHITE_NOISE:
            stm = 0

        elif self.process_enum == EnumNoiseProcess.RANDOM_CONSTANT:
            stm = 1

        elif self.process_enum == EnumNoiseProcess.RANDOM_WALK:
            stm = 1

        elif self.process_enum == EnumNoiseProcess.GAUSS_MARKOV:
            stm = np.exp(-time_step / self.correlation_time)

        else:
            raise ValueError(f"Invalid noise model: {self.process_enum}")

        return stm

    def get_process_noise(self, time_step):
        if self.process_enum == EnumNoiseProcess.WHITE_NOISE:
            # continuous-time Power Spectral Density q in [units^2 / s]
            q = np.multiply(self.process_noise, self.process_noise)
            process_noise = q * time_step

        elif self.process_enum == EnumNoiseProcess.RANDOM_CONSTANT:
            process_noise = 0

        elif self.process_enum == EnumNoiseProcess.RANDOM_WALK:
            # continuous-time sigma^2 of RW noise in [unit^2 / s]
            sigma_sq = np.multiply(self.process_noise, self.process_noise)
            process_noise = sigma_sq * time_step

        elif self.process_enum == EnumNoiseProcess.GAUSS_MARKOV:
            sigma_sq = np.multiply(self.process_noise, self.process_noise)
            alfa = np.exp(-time_step / self.correlation_time)
            process_noise = sigma_sq * (1 - alfa**2)

        else:
            raise ValueError(f"Invalid noise model: {self.process_enum}")

        return process_noise
