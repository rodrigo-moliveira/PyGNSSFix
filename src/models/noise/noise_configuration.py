""" Module with the definition of Noise Models (Stochastic Processes). """
import numpy as np

from src.io.config import EnumNoiseProcess
from .noise_gen import *


class NoiseModel:
    """
    Represents a stochastic process noise model used in Kalman filtering for GNSS states.

    Supports common models such as White Noise, Random Walk, Random Constant, and Gauss-Markov,
    and provides methods to compute both the State Transition Matrix (STM) and the process noise
    covariance contribution over a given time step.

    Attributes:
        name (str): The state name.
        process_enum (EnumNoiseProcess): Enum representation of the noise model type.
        process_gen (object): Instance of the specific process model class (used for future extensions).
        process_noise (float): Standard deviation input (see constructor documentation).
        correlation_time (float): Correlation time tau in seconds for Gauss-Markov models.

    Methods:
        gen():
        get_stm_entry(time_step):
        get_process_noise(time_step):
    """
    def __init__(self, state_str, process_str, process_noise=None, correlation_time=None):
        """
        Constructor of NoiseModel instances.

        Args:
            state_str (str): Name of the associated state (e.g., 'position', 'clock_bias', etc.).
            process_str (str): Type of noise model to use ('white_noise', 'random_walk', etc.).
            process_noise (float, optional): Standard deviation of the noise process.
                - For White Noise and Random Walk: units of `[unit / sqrt(s)]` (continuous-time).
                - For Gauss-Markov: units of `[unit]`.
            correlation_time (float, optional): Time constant tau for Gauss-Markov models in seconds.
                Required only for Gauss-Markov.

        Raises:
            ValueError: If the provided noise model type is not recognized.
        """
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
        """ Placeholder method for generating noise samples (currently not implemented). """
        pass

    def get_stm_entry(self, time_step: float):
        """
        Returns the scalar state transition matrix STM for the process.

        Args:
            time_step(float): time step for the discretization of the continuous-time process

        Raises:
            ValueError: If the provided noise model type is not recognized.
        """
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
        """
        Computes discrete-time process noise covariance Q.

        Args:
            time_step(float): time step for the discretization of the continuous-time process

        Raises:
            ValueError: If the provided noise model type is not recognized.
        """
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
