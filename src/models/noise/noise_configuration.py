""" Module with the definition of Noise Models (Stochastic Processes). """
import numpy as np
from src.io.config import EnumNoiseProcess
from .noise_gen import *


class NoiseModel:
    """
    Represents a stochastic process noise model for an estimated state.

    ``NoiseModel`` provides a common interface for configuring and
    generating several types of stochastic processes used to model
    sensor errors or process noise in estimation algorithms.

    The supported stochastic processes are:

        - ``"white_noise"``
        - ``"random_walk"``
        - ``"random_constant"``
        - ``"gauss_markov"``
        - ``"constant"``

    Continuous-time and discrete-time parameters
    ---------------------------------------------
    The ``white_noise``, ``random_walk``, and ``gauss_markov`` processes
    are defined using continuous-time statistical parameters. These
    parameters are converted internally to the appropriate
    discrete-time representation when :meth:`gen` is called.

    The ``constant`` and ``random_constant`` processes are defined
    directly in discrete time and therefore do not require a sampling
    interval for their statistical definition.

    Process parameter conventions
    ------------------------------
    ``process_noise`` has a different interpretation depending on the
    selected process:

        White noise:
            ``process_noise`` is sqrt(PSD) of the continuous-time white
            noise.

            Units:
                [U * sqrt(s)]

            where [U] is the unit of the process output.

        Random walk:
            ``process_noise`` is sqrt(PSD) of the continuous-time
            driving white noise.

            Units:
                [U / sqrt(s)]

        Gauss-Markov:
            ``process_noise`` is sqrt(PSD) of the continuous-time
            driving white noise.

            Units:
                [U / sqrt(s)]

            ``correlation_time`` specifies the correlation time tau
            in seconds.

        Random constant:
            ``process_noise`` is the discrete-time standard deviation
            of the randomly generated constant.

            Units:
                [U]

        Constant:
            ``process_noise`` is the discrete-time constant value.

            Units:
                [U]

    Here, [U] denotes the units of the state or measurement being
    modeled.

    For the continuous-time processes, the PSD is obtained internally
    as

        PSD = process_noise ** 2.

    The corresponding discrete-time process noise is then computed
    using the sampling interval supplied to :meth:`gen`.

    Attributes:
        name (str):
            Name of the associated state or error model.

        process_enum (EnumNoiseProcess):
            Enum identifying the selected stochastic process.

        process_gen (NoiseProcess):
            Instance of the stochastic process generator associated
            with ``process_enum``.

        process_noise (float or array-like):
            Statistical parameter defining the strength or value of
            the selected stochastic process. Its meaning and units
            depend on the selected process as described above.

        correlation_time (float or array-like, optional):
            Correlation time tau in seconds for a Gauss-Markov process.
            ``None`` for processes that do not use a correlation time.

        relative_re_param (float):
            Relative covariance increase factor used when a state is
            re-parameterized, for example when a GNSS state changes
            parameterization after a pivot change or when a satellite
            re-enters the state vector.

            The default value is 1.0, meaning no covariance increase.

    Args:
        state_str (str):
            Name of the associated state or error model. Examples
            include ``"position"``, ``"clock_bias"``, ``"gyroscope_bias"``,
            or ``"accelerometer_scale_factor"``.

        process_str (str):
            Type of stochastic process to use. Supported values are:

                ``"white_noise"``
                ``"random_walk"``
                ``"random_constant"``
                ``"gauss_markov"``
                ``"constant"``

        process_noise (float or array-like, optional):
            Process statistical parameter. The interpretation depends

        correlation_time (float or array-like, optional):
            Correlation time tau in seconds for a Gauss-Markov process.
            Required when ``process_str`` is ``"gauss_markov```.
            Ignored for other process types.

        init (float or array-like, optional):
            Initial value of the stochastic process. Used by processes
            that support an initial condition, such as Random Walk and
            Gauss-Markov.

    Raises:
        ValueError:
            If ``process_str`` does not correspond to a supported
            ``EnumNoiseProcess`` value.

    """
    def __init__(self, state_str, process_str, process_noise=None, correlation_time=None, init=None):
        """
        Initialize a stochastic noise model.

        Args:
            state_str (str):
                Name of the associated state or error model.

            process_str (str):
                Type of stochastic process to use. Supported values are
                ``"white_noise"``, ``"random_walk"``,
                ``"random_constant"``, ``"gauss_markov"``, and
                ``"constant"``.

            process_noise (float or array-like, optional):
                Statistical parameter defining the process. For
                continuous-time processes this is sqrt(PSD); for
                discrete-time processes this is the standard deviation
                or value of the process, depending on the model.

            correlation_time (float or array-like, optional):
                Correlation time tau in seconds. Required for a
                Gauss-Markov process.

            init (float or array-like, optional):
                Initial value of the stochastic process.

        Raises:
            ValueError:
                If ``process_str`` does not correspond to a supported
                stochastic process.
        """
        self.name = state_str
        self.process_enum = EnumNoiseProcess.init_model(process_str)

        if self.process_enum == EnumNoiseProcess.WHITE_NOISE:
            self.process_gen = WhiteNoiseProcess(psd=process_noise)
        elif self.process_enum == EnumNoiseProcess.RANDOM_CONSTANT:
            self.process_gen = RandomConstantProcess(std=process_noise)
        elif self.process_enum == EnumNoiseProcess.RANDOM_WALK:
            self.process_gen = RandomWalkProcess(psd=process_noise, init=init)
        elif self.process_enum == EnumNoiseProcess.GAUSS_MARKOV:
            self.process_gen = GaussMarkovProcess(psd=process_noise,correlation_time=correlation_time, init=init)
        elif self.process_enum == EnumNoiseProcess.CONSTANT:
            self.process_gen = ConstantProcess(val=process_noise)
        else:
            raise ValueError(f"Expected EnumNoiseProcess enum, got {type(self.process_enum).__name__}")

        # internal variables
        self.process_noise = process_noise  # units of std (RW/WN - [units/sqrt(s)] or GM[units])
        self.correlation_time = correlation_time  # in [s]
        self.relative_re_param = 1.0

    def gen(self, n, sampling_time):
        """
        Generate a discrete-time realization of the configured process.

        For continuous-time processes (White Noise, Random Walk, and
        Gauss-Markov), the continuous-time statistical parameters stored
        in the model are discretized using ``sampling_time``.

        For Constant and Random Constant processes, the parameters are
        already defined in discrete time.

        Args:
            n (int):
                Number of discrete-time samples to generate.

            sampling_time (float):
                Sampling interval between consecutive samples, in
                seconds.

        Returns:
            numpy.ndarray:
                Generated process realization with shape
                ``(n, axis)``.

        """
        return self.process_gen.compute(n, sampling_time=sampling_time)

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
            # sigma_sq is the continuous-time sigma^2 of RW noise in [unit^2 / s]
            sigma_sq = np.multiply(self.process_noise, self.process_noise)
            process_noise = sigma_sq * time_step

        elif self.process_enum == EnumNoiseProcess.GAUSS_MARKOV:
            sigma_sq = np.multiply(self.process_noise, self.process_noise)
            alfa = np.exp(-time_step / self.correlation_time)
            process_noise = sigma_sq * (1 - alfa**2)

        else:
            raise ValueError(f"Invalid noise model: {self.process_enum}")

        return process_noise
