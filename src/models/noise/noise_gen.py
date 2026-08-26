""" This module contains the stochastic process classes """
import numpy as np

__all__ = ["RandomWalkProcess", "GaussMarkovProcess", "WhiteNoiseProcess", "RandomConstantProcess", "ConstantProcess"]

class NoiseProcess:
    """
    Base class for stochastic noise processes.

    This class defines the common interface for noise-process
    generators. Subclasses should implement the :meth:`compute`
    method to generate a realization of the corresponding stochastic
    process.


    Attributes:
        _name (str): Name of the stochastic process.
        _axis (int): Number of process axes.

    Notes
    -----
    The expected output of a noise process is a NumPy array with shape

        (n, axis)

    where ``n`` is the number of discrete-time samples and ``axis`` is
    the number of process dimensions.

    Subclasses can extend the `compute` method with additional
    parameters required by the specific stochastic process.
    """
    def __init__(self, axis=1):
        """
        Initialize a noise process.

        Args:
            axis (int, optional): Number of independent axes (or dimensions) of the process.
                The generated output has one column per axis.
        """
        self._name = "General Process"
        self._axis = axis

    def __str__(self):
        """
        Return the name of the stochastic process.

        Returns:
            str: Name of the stochastic process.
        """
        return self._name

    def compute(self, n, *args):
        """
        Generate a realization of the stochastic process.

        This base implementation returns an array of zeros. Subclasses
        should override this method to implement the corresponding
        stochastic process.

        Args:
            n (int): Number of discrete-time samples to generate.
            *args (list): Additional parameters required by a specific process.

        Returns:
            numpy.ndarray: Array of shape ``(n, axis)`` containing the generated process realization.

        Notes
        -----
        The base implementation represents a deterministic zero process.
        """
        return np.zeros((n, self._axis))


class RandomWalkProcess(NoiseProcess):
    """
    Continuous-time random walk stochastic process.

    The process is modeled as

        dx(t) / dt = w(t)

    where w(t) is zero-mean white Gaussian noise with continuous-time
    power spectral density q.

    For a discrete sampling interval dt, the process is discretized as

        x_k = x_{k-1} + n_k

    where

        n_k ~ N(0, Q_d)

    and

        Q_d = q * dt.

    Therefore, for each axis i, the random increment is Gaussian with

        n_k[i] ~ N(0, q[i] * dt)

    and standard deviation

        sigma[i] = sqrt(q[i] * dt).

    The random increments are independent between epochs and between
    axes.

    Notes
    -----
    The process variance grows linearly with time. For a constant
    sampling interval dt and initial deterministic state,

        Var[x_k] = k * q * dt

    or, equivalently, after elapsed time t,

        Var[x(t)] = q * t.

    Hence, the standard deviation of the random walk grows as

        std[x(t)] = sqrt(q * t).
    """

    def __init__(self, psd=1, init=None):
        """
        Constructor of Continuous-time random walk stochastic process.

        Args:
            psd (float or list): Continuous-time power spectral density of the random walk.
                One value is specified per axis. The units depend on the units
                of the state; for a state in [U], the units are [U^2/s].
                Therefore, the name ``psd`` stands for a continuous-time
                spectral density rather than a discrete-time variance.

            init (float or list): Initial value of the random walk for each axis. If None, all
                axes are initialized to zero. Units are [U]
        """
        # psd is continuous time power spectral density
        if isinstance(psd, float) or isinstance(psd, int):
            psd = [float(psd)]
        if isinstance(init, float) or isinstance(init, int):
            init = [float(init)]

        axis = len(psd)
        super().__init__(axis=axis)

        self._name = "Random Walk Process"
        self._psd = np.asarray(psd, dtype=float)

        if init is not None:
            if len(init) != axis:
                raise TypeError(f"Initial value for the stochastic process has dim={len(init)}, which is not "
                                f"consistent with the provided axes length={axis}")
            self._init = init
        else:
            self._init = np.zeros(axis)

    def compute(self, n, sampling_time=1.0, *args):
        """
        Generate a discrete-time realization of the random walk.

        Args:
            n (int): Number of discrete time samples to generate.
            sampling_time (float): Sampling interval between consecutive samples, in seconds.
                The default is 1 second.

        Returns:
            numpy.ndarray: Array of shape (n, axis) containing the random-walk
                realization. ``walk[k, i]`` is the value of axis ``i`` at
                discrete epoch ``k``.

        Notes
        At every epoch, an independent Gaussian increment is generated:

            n_k ~ N(0, q * sampling_time)

        and added to the previous state:

            x_k = x_{k-1} + n_k.

        Thus, for each axis i,

            E[n_k[i]] = 0
            Var[n_k[i]] = q[i] * sampling_time
            Std[n_k[i]] = sqrt(q[i] * sampling_time).

        The increments are independent between epochs and axes.

        The variance of the accumulated random walk increases linearly
        with elapsed time:

            Var[x_k] = Var[x_0] + k * q * sampling_time

        for deterministic initial conditions.
        """
        walk = np.zeros((n, self._axis))

        # set initial condition
        walk[0, :] = self._init

        # time loop (start at t = 1)
        for t in range(1, n):
            # generate gaussian noise for this epoch and all required axes
            step_n = np.random.normal(scale=np.sqrt(self._psd * sampling_time),size=self._axis)
            walk[t, :] = walk[t-1, :] + step_n

        return walk


class WhiteNoiseProcess(NoiseProcess):
    """
    Continuous-time white Gaussian noise process.

    The process is modeled as zero-mean white Gaussian noise with
    continuous-time power spectral density (PSD) ``q``:

        E[w(t)] = 0

        E[w(t) w(tau)] = q * delta(t - tau)

    For a discrete sampling interval ``dt``, the continuous-time PSD
    is converted to the standard deviation of the discrete-time
    samples according to

        sigma_d^2 = q / dt

    and therefore

        sigma_d = sqrt(q / dt).

    Each generated sample is independent of all other samples and of
    the samples generated on other axes.

    Notes
    -----
    The discrete-time samples are distributed as

        x_k ~ N(0, q / dt)

    where ``q`` is the continuous-time PSD and ``dt`` is the sampling
    interval.

    Consequently,

        E[x_k] = 0

        Var[x_k] = q / dt

        Std[x_k] = sqrt(q / dt).

    As the sampling interval becomes smaller, the variance of each
    discrete sample increases. This is required to preserve the same
    continuous-time noise spectral density.
    """

    def __init__(self, psd=1):
        """
        Constructor of Continuous-time white Gaussian noise process.

        Args:
            psd (float or list): Continuous-time power spectral density of the white noise.
                One value is specified per axis. If the process units are [U], the PSD has units U^2 * s.
        """
        if isinstance(psd, (float, int)):
            psd = [float(psd)]

        self._psd = np.asarray(psd, dtype=float)

        axis = len(self._psd)
        super().__init__(axis=axis)

        self._name = "White Noise Process"

    def compute(self, n, sampling_time=1.0, *args):
        """
        Generate a discrete-time realization of the white noise.

        Args:
            n (int): Number of discrete time samples to generate.

            sampling_time (float): Sampling interval between consecutive samples, in seconds.
                Default is 1 second.

        Returns:
            numpy.ndarray : Array of shape ``(n, axis)`` containing independent
                zero-mean Gaussian white-noise samples.

        Notes

        For each axis ``i`` and epoch ``k``:

            x[k, i] ~ N(0, psd[i] / sampling_time)

        Therefore, the standard deviation used to generate each sample
        is

            sigma[i] = sqrt(psd[i] / sampling_time).
        """

        std = np.sqrt(self._psd / sampling_time)

        return np.random.normal(
            scale=std,
            size=(n, self._axis)
        )


class RandomConstantProcess(NoiseProcess):
    """
    Generate a random constant process.

    A single random value is generated independently for each axis
    according to a zero-mean Gaussian distribution. The generated value
    is then kept constant for all samples.

    Args:
        std (float or int or list): Standard deviation of the Gaussian
            distribution used to generate the constant value for each
            axis. Units are [U].

    """

    def __init__(self, std=1):
        """
        Initialize a random constant process.

        Args:
            std (float or int or list): Standard deviation of the
                Gaussian distribution for each axis. Units are [U].
        """
        # std is the discrete-time standard deviation
        if isinstance(std, (float, int)):
            std = [float(std)]

        axis = len(std)
        super().__init__(axis=axis)

        self._name = "Random Constant Process"
        self._std = std

    def compute(self, n, *args):
        """
        Generate a realization of the random constant process.

        A single zero-mean Gaussian random value is generated for each
        axis and repeated for all ``n`` samples.

        Args:
            n (int): Number of discrete-time samples to generate.
            *args: Additional arguments, unused by this process.

        Returns:
            numpy.ndarray: Array of shape ``(n, axis)`` containing the
                random constant process. Each column contains the same
                randomly generated value for all samples.
        """
        std = np.random.normal(scale=self._std)
        return np.ones((n, self._axis)) * std

class ConstantProcess(NoiseProcess):
    """
    Generate a deterministic constant process.

    The specified value for each axis is repeated for all samples.

    Args:
        val (float or int or list): Constant value for each axis.
            Units are [U].

    """

    def __init__(self, val=1):
        """
        Initialize a constant process.

        Args:
            val (float or int or list): Constant value for each axis.
                Units are [U].
        """
        if isinstance(val, (float, int)):
            val = [float(val)]

        axis = len(val)
        super().__init__(axis=axis)

        self._name = "Constant Process"
        self._val = val

    def compute(self, n, *args):
        """
        Generate a realization of the constant process.

        Args:
            n (int): Number of discrete-time samples to generate.
            *args: Additional arguments, unused by this process.

        Returns:
            numpy.ndarray: Array of shape ``(n, axis)`` containing the
                constant process. Each column contains the corresponding
                value specified in ``val`` for all samples.
        """
        return np.ones((n, self._axis)) * self._val


class GaussMarkovProcess(NoiseProcess):
    """
    First-order Gauss-Markov stochastic process.

    The continuous-time process is defined by

        dx(t) / dt = -beta * x(t) + w(t)

    where

        beta = 1 / correlation_time

    and ``w(t)`` is zero-mean white Gaussian noise with continuous-time
    power spectral density ``psd``.

    The exact discrete-time representation for a sampling interval dt is

        x_k = phi * x_{k-1} + eta_k

    where

        phi = exp(-beta * dt)

    and

        eta_k ~ N(0, Q_d)

    with

        Q_d = (psd / (2 * beta)) * (1 - exp(-2 * beta * dt)).

    Equivalently,

        Q_d = (psd * correlation_time / 2)
              * (1 - exp(-2 * dt / correlation_time)).

    The stationary variance of the process is

        sigma_x^2 = psd / (2 * beta)
                  = psd * correlation_time / 2.

    Notes
    -----
    The autocorrelation of the stationary process decays exponentially:

        R(tau) = sigma_x^2 * exp(-|tau| / correlation_time)

    where

        sigma_x^2 = psd * correlation_time / 2.

    Therefore, ``correlation_time`` controls the temporal correlation,
    while ``psd`` controls the strength of the driving white noise.

    The generated process is Gaussian. If initialized in the stationary
    distribution, its marginal distribution is

        x(t) ~ N(0, sigma_x^2).

    A deterministic initial condition, however, results in a transient
    before the process reaches its stationary distribution.
    """

    def __init__(self,psd=1,correlation_time=1,init=None):
        """
        Constructor of First-order Gauss-Markov stochastic process.

        Args:
            psd (float or list): Continuous-time power spectral density of the driving white
                noise. One value is specified per axis.
                If the state units are [U], psd has units [U^2/s]

            correlation_time (float): Correlation time of the process in seconds. The correlation time
                determines how quickly the process loses correlation.

            init (float or list): Initial value of the process for each axis. If None, all axes
                are initialized to zero.
                Units are in [U]
        """
        if isinstance(psd, (float, int)):
            psd = [float(psd)]
        if isinstance(init, float) or isinstance(init, int):
            init = [float(init)]
        axis = len(psd)

        super().__init__(axis=axis)
        self._name = "Gauss-Markov Process"

        self._psd = np.asarray(psd, dtype=float)
        self._correlation_time = correlation_time

        if init is not None:
            if len(init) != axis:
                raise TypeError(
                    f"Initial value for the stochastic process has "
                    f"dim={len(init)}, which is not consistent with "
                    f"the provided axes length={axis}"
                )
            self._init = np.asarray(init, dtype=float)
        else:
            self._init = np.zeros(axis)

    def compute(self, n, sampling_time=1.0, *args):
        """
        Generate a discrete-time realization of the Gauss-Markov process.

        Args:
            n (int): Number of discrete time samples to generate.

            sampling_time (float): Sampling interval between consecutive samples, in seconds.

        Returns:
            numpy.ndarray : Array of shape ``(n, axis)`` containing the generated Gauss-Markov process.

        Notes

        For each axis, the discrete-time process is generated according
        to

            x_k = phi * x_{k-1} + eta_k

        where

            phi = exp(-dt / correlation_time)

        and

            eta_k ~ N(0, Q_d)

        with

            Q_d = (psd * correlation_time / 2)
                  * (1 - exp(-2 * dt / correlation_time)).

        Thus, the standard deviation of the Gaussian innovation is

            sigma_eta = sqrt(Q_d).
        """

        gm_process = np.zeros((n, self._axis))

        if n == 0:
            return gm_process

        gm_process[0, :] = self._init

        dt = sampling_time
        tau = np.asarray(self._correlation_time, dtype=float)

        phi = np.exp(-dt / tau)

        innovation_std = np.sqrt(self._psd * tau / 2.0 * (1.0 - np.exp(-2.0 * dt / tau)))

        for k in range(1, n):
            noise = np.random.normal(loc=0.0, scale=innovation_std, size=self._axis)
            gm_process[k, :] = (phi * gm_process[k - 1, :] + noise)

        return gm_process