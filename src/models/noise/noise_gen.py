import numpy as np

class NoiseProcess:
    def __init__(self, axis=1):
        self._name = "General Process"
        self._axis = axis

    def __str__(self):
        return self._name

    def compute(self, n, *args):
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

    Parameters
    ----------
    var : float or array-like of float, optional
        Continuous-time power spectral density of the random walk.
        One value is specified per axis. The units depend on the units
        of the state; for a state in meters, the units are m^2/s.
        Despite the name ``var``, this quantity is a continuous-time
        spectral density rather than a discrete-time variance.
    init : array-like of float, optional
        Initial value of the random walk for each axis. If None, all
        axes are initialized to zero.

    Notes
    -----
    The process variance grows linearly with time. For a constant
    sampling interval dt and initial deterministic state,

        Var[x_k] = k * q * dt

    or, equivalently, after elapsed time t,

        Var[x(t)] = q * t.

    Hence the standard deviation of the random walk grows as

        std[x(t)] = sqrt(q * t).
    """

    def __init__(self, psd=1, init=None):
        # psd is continuous time power spectral density
        if isinstance(psd, float) or isinstance(psd, int):
            psd = [float(psd)]

        axis = len(psd)
        super().__init__(axis=axis)

        self._name = "Random Walk Process"
        self._psd = psd

        if init is not None:
            if len(init) != axis:
                raise TypeError(f"Initial value for the stochastic process has dim={len(init)}, which is not "
                                f"consistent with the provided axes length={axis}")
            self._init = init
        else:
            self._init = np.zeros(axis)

    def compute(self, n, sampling_time=1, *args):
        """
        Generate a discrete-time realization of the random walk.

        Parameters
        ----------
        n : int
            Number of discrete time samples to generate.
        sampling_time : float, optional
            Sampling interval between consecutive samples, in seconds.
            The default is 1 second.

        Returns
        -------
        numpy.ndarray
            Array of shape (n, axis) containing the random-walk
            realization. ``walk[k, i]`` is the value of axis ``i`` at
            discrete epoch ``k``.

        Notes
        -----
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

    Parameters
    ----------
    psd : float or array-like of float, optional
        Continuous-time power spectral density of the white noise.
        One value is specified per axis.

        If the process units are U, the PSD has units U^2 / s.

        For example, if the process represents a sensor error in
        meters, ``psd`` has units m^2/s.

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
        if isinstance(psd, (float, int)):
            psd = [float(psd)]

        self._psd = np.asarray(psd, dtype=float)

        axis = len(self._psd)
        super().__init__(axis=axis)

        self._name = "White Noise Process"

    def compute(self, n, sampling_time=1, *args):
        """
        Generate a discrete-time realization of the white noise.

        Parameters
        ----------
        n : int
            Number of discrete time samples to generate.

        sampling_time : float, optional
            Sampling interval between consecutive samples, in seconds.
            Default is 1 second.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(n, axis)`` containing independent
            zero-mean Gaussian white-noise samples.

        Notes
        -----
        For each axis ``i`` and epoch ``k``:

            x[k, i] ~ N(0, psd[i] / sampling_time)

        Therefore the standard deviation used to generate each sample
        is

            sigma[i] = sqrt(psd[i] / sampling_time).
        """

        std = np.sqrt(self._psd / sampling_time)

        return np.random.normal(
            scale=std,
            size=(n, self._axis)
        )


class RandomConstantProcess(NoiseProcess):

    def __init__(self, std=1):
        """

        Args:
            std (float or int or list):
        """
        # std is discrete time standard deviation
        if isinstance(std, float) or isinstance(std, int):
            std = [float(std)]

        axis = len(std)
        super().__init__(axis=axis)
        self._name = "Random Constant Process"
        self._std = std

    def compute(self, n, *args):
        std = np.random.normal(scale=self._std)
        return np.ones((n, self._axis)) * std

class ConstantProcess(NoiseProcess):
    def __init__(self, val=1):
        """

        Args:
            val (float or int or list):
        """
        # std is discrete time standard deviation
        if isinstance(val, float) or isinstance(val, int):
            val = [float(val)]

        axis = len(val)
        super().__init__(axis=axis)
        self._name = "Constant Process"
        self._val = val

    def compute(self, n, *args):
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

    Parameters
    ----------
    dim : int, optional
        Number of samples to generate.
    psd : float or array-like of float, optional
        Continuous-time power spectral density of the driving white
        noise. One value is specified per axis. The units are
        state_unit^2 / s.
    correlation_time : float or array-like of float, optional
        Correlation time of the process in seconds. The correlation time
        determines how quickly the process loses correlation.
    axis : int, optional
        Number of independent process axes.
    initial : array-like of float, optional
        Initial value of the process for each axis. If None, all axes
        are initialized to zero.

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

    def __init__(self,psd=1,correlation_time=1,init=None,
    ):
        self._name = "Gauss-Markov Process"
        if isinstance(psd, (float, int)):
            psd = [float(psd)]
        axis = len(psd)

        super().__init__(axis=axis)

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

    def compute(self, n, sampling_time=1, *args):
        """
        Generate a discrete-time realization of the Gauss-Markov process.

        Parameters
        ----------
        n : int
            Number of discrete time samples to generate.

        sampling_time : float
            Sampling interval between consecutive samples, in seconds.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(dim, axis)`` containing the generated
            Gauss-Markov process.

        Notes
        -----
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