import numpy as np
from scipy import signal


def power_spectral_density(x, fs, nperseg=None):
    """
    Estimate the one-sided power spectral density of a signal.

    The power spectral density (PSD) is estimated using Welch's method.
    For multi-dimensional input, the PSD is independently computed along
    each column.

    Args:
        x (numpy.ndarray):
            Input signal with shape ``(n_samples,)`` or
            ``(n_samples, n_axes)``. The samples are assumed to be
            uniformly spaced in time.

        fs (float):
            Sampling frequency in Hz.

        nperseg (int, optional):
            Number of samples per segment used by Welch's method.
            If ``None``, the default value used by
            ``scipy.signal.welch`` is applied.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]:
            Estimated one-sided power spectral density and corresponding
            frequency vector.

            ``psd`` has units ``[U^2/Hz]``, where ``[U]`` is the unit
            of the input signal.

            ``f`` has units ``[Hz]``.

            For multi-axis input with shape ``(n_samples, n_axes)``,
            the PSD has shape ``(n_frequencies, n_axes)``.
    """
    f, psd = signal.welch(
        x,
        fs=fs,
        return_onesided=True,
        scaling="density",
        nperseg=nperseg,
        axis=0,
    )

    return psd, f


def _get_taus(n, ts):
    """
    Generate averaging factors and corresponding averaging times.

    The averaging factors are approximately logarithmically spaced as

        1, 2, ..., 9, 10, 20, ..., 90, 100, 200, ...

    The maximum averaging factor is selected such that at least nine
    non-overlapping bins are available for estimating the Allan variance.

    Args:
        n (int):
            Number of samples in the input time series.

        ts (float):
            Sampling period in seconds.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]:
            A tuple containing:

                - multipliers (numpy.ndarray): Integer number of samples
                  per averaging bin.

                - tau (numpy.ndarray): Corresponding averaging times in
                  seconds.
    """
    max_samples_per_bin = n // 9

    if max_samples_per_bin < 1:
        return np.array([], dtype=int), np.array([])

    multipliers = []

    scale = 1

    while scale <= max_samples_per_bin:
        for j in range(1, 10):
            multiplier = j * scale

            if multiplier > max_samples_per_bin:
                break

            multipliers.append(multiplier)

        scale *= 10

    multipliers = np.array(multipliers, dtype=int)
    tau = multipliers * ts

    return multipliers, tau


def allan_variance(x, fs):
    """
    Estimate the non-overlapping Allan variance of a multi-axis signal.

    The Allan variance is independently computed for each axis of the
    input signal using non-overlapping averaging bins.

    For an averaging factor ``m``, the averaging time is

        tau = m * Ts

    where ``Ts`` is the sampling period.

    The input signal is divided into ``M`` consecutive non-overlapping
    bins, each containing ``m`` samples. The average value of each bin is

        x_bar[k] = (1 / m) * sum(x[i])

    where the summation is performed over the samples belonging to the
    k-th bin.

    The non-overlapping Allan variance is then computed as

        sigma_A^2(tau) =
            1 / (2 * (M - 1))
            * sum(
                (x_bar[k + 1] - x_bar[k])**2
              )

    where ``M`` is the number of non-overlapping bins.

    Args:
        x (numpy.ndarray):
            Uniformly sampled input signal with shape
            ``(n_samples, n_axes)``.

        fs (float):
            Sampling frequency in Hz.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]:
            A tuple containing:

                - allan_var (numpy.ndarray): Non-overlapping Allan
                  variance with shape ``(n_tau, n_axes)``.

                - tau (numpy.ndarray): Averaging times corresponding to
                  the Allan variance estimates, in seconds.

    Raises:
        ValueError:
            If ``x`` is not a two-dimensional array or if ``fs`` is not
            positive.
    """
    x = np.asarray(x)

    if x.ndim != 2:
        raise ValueError(
            "Input signal must have shape (n_samples, n_axes)."
        )

    if fs <= 0:
        raise ValueError(
            "Sampling frequency must be positive."
        )

    n_samples, n_axes = x.shape

    ts = 1.0 / fs

    # Get averaging factors and corresponding averaging times.
    multipliers, tau = _get_taus(n_samples, ts)

    allan_var = np.zeros((len(multipliers), n_axes))

    # Compute Allan variance independently for each axis.
    for axis in range(n_axes):
        allan_var[:, axis] = _allan_variance(
            x[:, axis],
            multipliers,
        )

    return allan_var, tau


def _allan_variance(x, multipliers):
    """
    Compute the non-overlapping Allan variance of a one-dimensional signal.

    For each averaging factor ``m``, the input signal is divided into
    consecutive non-overlapping bins containing ``m`` samples. The mean
    value of each bin is computed and the Allan variance is obtained
    from the squared differences between consecutive bin averages.

    The averaging time associated with an averaging factor ``m`` is

        tau = m * Ts

    The average of the k-th bin is

        x_bar[k] =
            (1 / m)
            * sum(x[k * m : (k + 1) * m])

    Given ``M`` non-overlapping bins, the Allan variance is

        sigma_A^2(tau) =
            1 / (2 * (M - 1))
            * sum(
                (x_bar[k + 1] - x_bar[k])**2
              )

    Args:
        x (numpy.ndarray):
            One-dimensional uniformly sampled input signal.

        multipliers (numpy.ndarray):
            Integer averaging factors representing the number of samples
            contained in each non-overlapping averaging bin.

    Returns:
        numpy.ndarray:
            Non-overlapping Allan variance corresponding to each
            averaging factor in ``multipliers``. The returned array has
            shape ``(len(multipliers),)``.

    Notes:
        Any samples remaining after dividing the signal into an integer
        number of bins are discarded.

        The averaging factors are expected to be selected such that a
        sufficient number of bins is available for each Allan variance
        estimate.
    """
    x = np.asarray(x)

    n = len(x)
    allan_var = np.zeros(len(multipliers))

    for i, samples_per_bin in enumerate(multipliers):

        # Number of complete non-overlapping bins.
        n_bins = n // samples_per_bin

        if n_bins < 2:
            allan_var[i] = np.nan
            continue

        # Discard incomplete samples and divide the signal into
        # non-overlapping bins.
        bins = x[:n_bins * samples_per_bin].reshape(
            n_bins,
            samples_per_bin,
        )

        # Compute the mean value of each bin.
        bin_means = np.mean(bins, axis=1)

        # Difference between consecutive bin averages.
        diff = bin_means[1:] - bin_means[:-1]

        # Non-overlapping Allan variance.
        allan_var[i] = 0.5 * np.mean(diff**2)

    return allan_var
