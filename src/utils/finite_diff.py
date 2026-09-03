import numpy as np
import pandas
from src.errors import FiniteDifferenceError

class NumericalDifferentiator:
    """Numerically differentiate data contained in a pandas DataFrame.

        The DataFrame is expected to contain a time column and one or more data
        columns. The derivative is computed at a requested time using one of the
        available finite-difference methods.

        Supported methods are:

            - ``"forward-two-point"``
            - ``"central-two-point"``
            - ``"forward-three-point"``
            - ``"backward-three-point"``
            - ``"central-five-point"``
            - ``"optimal"``

        The ``"optimal"`` method automatically selects the finite-difference
        scheme based on the position of the requested sample:

            - First and second samples:
                forward three-point difference.
            - Last and second-to-last samples:
                backward three-point difference.
            - All other samples:
                central five-point difference.

        Notes:
            The finite-difference formulas assume uniformly spaced time samples.
            The returned derivative has the same shape as the selected data
            columns. If the input data columns contain three-dimensional vectors
            represented as separate columns, differentiation is performed
            element-wise.

        Args:
            method: Finite-difference method to use.

        Raises:
            FiniteDifferenceError: If ``method`` is not supported.
        """

    _VALID_METHODS = (
        "forward-two-point",
        "central-two-point",
        "forward-three-point",
        "backward-three-point",
        "central-five-point",
        "optimal",
    )

    def __init__(self, method: str) -> None:
        """Initialize the numerical differentiator.

        Args:
            method: Finite-difference method to use. Supported values are
                ``"forward-two-point"``, ``"central-two-point"``,
                ``"forward-three-point"``, ``"backward-three-point"``,
                ``"central-five-point"``, and ``"optimal"``.

        Raises:
            FiniteDifferenceError: If ``method`` is not supported.
        """
        if method not in self._VALID_METHODS:
            raise FiniteDifferenceError(
                "method must be one of: "
                + ", ".join(self._VALID_METHODS)
            )

        self.method = method

    def compute(
            self,
            df: pandas.DataFrame,
            time_col: int,
            data_cols: list[int] | None = None,
            time: float = -1,
    ) -> float | np.ndarray:
        """Compute the numerical derivative at a specified time.

        The derivative is calculated using the method selected when the
        differentiator was initialized.

        Args:
            df: DataFrame containing the time samples and data to
                differentiate.
            time_col: Zero-based positional index of the column containing
                the time values.
            data_cols: Zero-based positional indices of the columns to
                differentiate. If ``None``, all columns except ``time_col``
                are differentiated.
            time: Time at which the derivative is requested. The value must
                occur exactly once in the time column.

        Returns:
            A NumPy array containing the derivative of the selected data
            columns at ``time``.

        Raises:
            FiniteDifferenceError: If the input DataFrame, column indices,
                time vector, requested time, or number of samples is invalid.
            ValueError: If the time samples are not uniformly spaced.
        """
        # ------------------------------------------------------------------
        # Validate DataFrame
        # ------------------------------------------------------------------
        if not isinstance(df, pandas.DataFrame):
            raise FiniteDifferenceError(
                "df must be a pandas.DataFrame"
            )

        n_samples, n_columns = df.shape

        # ------------------------------------------------------------------
        # Validate time column
        # ------------------------------------------------------------------
        if not isinstance(time_col, (int, np.integer)):
            raise FiniteDifferenceError(
                "time_col must be an integer"
            )

        if not 0 <= time_col < n_columns:
            raise FiniteDifferenceError(
                f"time_col must be between 0 and {n_columns - 1}, "
                f"got {time_col}"
            )

        # ------------------------------------------------------------------
        # Validate number of samples
        # ------------------------------------------------------------------
        required_samples = {
            "forward-two-point": 2,
            "central-two-point": 3,
            "forward-three-point": 3,
            "backward-three-point": 3,
            "central-five-point": 5,
            "optimal": 5,
        }

        min_samples = required_samples[self.method]

        if n_samples < min_samples:
            raise FiniteDifferenceError(
                f"method '{self.method}' requires at least "
                f"{min_samples} samples, got {n_samples}"
            )

        # ------------------------------------------------------------------
        # Extract and validate time vector
        # ------------------------------------------------------------------
        time_values = np.asarray(df.iloc[:, time_col], dtype=float)

        if not np.all(np.isfinite(time_values)):
            raise FiniteDifferenceError(
                "time column contains non-finite values"
            )

        time_steps = np.diff(time_values)

        if np.any(time_steps <= 0.0):
            raise FiniteDifferenceError(
                "time values must be strictly increasing"
            )

        step = float(time_steps[0])

        # Finite-difference formulas assume a constant step.
        if not np.allclose(time_steps, step):
            raise FiniteDifferenceError(
                "time samples must be uniformly spaced"
            )

        if step <= 0.0:
            raise FiniteDifferenceError(
                "time step must be positive"
            )

        # ------------------------------------------------------------------
        # Validate data columns
        # ------------------------------------------------------------------
        if data_cols is None:
            data_cols = [
                col for col in range(n_columns)
                if col != time_col
            ]
        else:
            data_cols = list(data_cols)

            if not data_cols:
                raise FiniteDifferenceError(
                    "data_cols must contain at least one column"
                )

            if any(
                    not isinstance(col, (int, np.integer))
                    for col in data_cols
            ):
                raise FiniteDifferenceError(
                    "all data_cols must be integers"
                )

            if any(
                    col < 0 or col >= n_columns
                    for col in data_cols
            ):
                raise FiniteDifferenceError(
                    f"data_cols must contain values between "
                    f"0 and {n_columns - 1}"
                )

            if time_col in data_cols:
                raise FiniteDifferenceError(
                    "time_col cannot be included in data_cols"
                )

        # ------------------------------------------------------------------
        # Find requested time using positional indices.
        #
        # np.flatnonzero() gives us positions, which is what .iloc expects.
        # ------------------------------------------------------------------
        matches = np.flatnonzero(time_values == time)

        if len(matches) == 0:
            raise FiniteDifferenceError(
                f"No matches found for time {time}. "
                f"Time array: {time_values}"
            )

        if len(matches) > 1:
            raise FiniteDifferenceError(
                f"Multiple matches found for time {time}. "
                f"Time array: {time_values}"
            )

        index = int(matches[0])

        # ------------------------------------------------------------------
        # Select the actual method for "optimal"
        # ------------------------------------------------------------------
        method = self.method

        if method == "optimal":
            if index in (0, 1):
                method = "forward-three-point"
            elif index in (n_samples - 2, n_samples - 1):
                method = "backward-three-point"
            else:
                method = "central-five-point"

        # ------------------------------------------------------------------
        # Compute finite difference
        # ------------------------------------------------------------------
        if method == "forward-two-point":
            if index + 1 >= n_samples:
                raise FiniteDifferenceError(
                    "forward-two-point cannot be used at the last sample"
                )

            f_x = np.asarray(df.iloc[index, data_cols])
            f_x_plus_h = np.asarray(df.iloc[index + 1, data_cols])

            return forward_two_point(
                f_x,
                f_x_plus_h,
                step,
            )

        elif method == "central-two-point":
            if index == 0 or index == n_samples - 1:
                raise FiniteDifferenceError(
                    "central-two-point cannot be used at the first "
                    "or last sample"
                )

            f_x_minus_h = np.asarray(df.iloc[index - 1, data_cols])
            f_x_plus_h = np.asarray(df.iloc[index + 1, data_cols])

            return central_two_point(
                f_x_plus_h,
                f_x_minus_h,
                step,
            )

        elif method == "forward-three-point":
            if index + 2 >= n_samples:
                raise FiniteDifferenceError(
                    "forward-three-point requires two samples after "
                    "the requested sample"
                )

            f_x = np.asarray(df.iloc[index, data_cols])
            f_x_plus_h = np.asarray(df.iloc[index + 1, data_cols])
            f_x_plus_2h = np.asarray(df.iloc[index + 2, data_cols])

            return forward_three_point(
                f_x,
                f_x_plus_h,
                f_x_plus_2h,
                step,
            )

        elif method == "backward-three-point":
            if index < 2:
                raise FiniteDifferenceError(
                    "backward-three-point requires two samples before "
                    "the requested sample"
                )

            f_x = np.asarray(df.iloc[index, data_cols])
            f_x_minus_h = np.asarray(df.iloc[index - 1, data_cols])
            f_x_minus_2h = np.asarray(df.iloc[index - 2, data_cols])

            return backward_three_point(
                f_x,
                f_x_minus_h,
                f_x_minus_2h,
                step,
            )

        elif method == "central-five-point":
            if index < 2 or index + 2 >= n_samples:
                raise FiniteDifferenceError(
                    "central-five-point requires two samples before "
                    "and after the requested sample"
                )

            f_x_minus_2h = np.asarray(df.iloc[index - 2, data_cols])
            f_x_minus_h = np.asarray(df.iloc[index - 1, data_cols])
            f_x_plus_h = np.asarray(df.iloc[index + 1, data_cols])
            f_x_plus_2h = np.asarray(df.iloc[index + 2, data_cols])

            return central_five_point(
                f_x_minus_2h,
                f_x_minus_h,
                f_x_plus_h,
                f_x_plus_2h,
                step,
            )

        else:
            # This should be unreachable because the method is validated
            # during initialization.
            raise FiniteDifferenceError(
                f"Unsupported finite-difference method: {method}"
            )



def _validate_step(step) -> None:
    """Validate a finite-difference step size.

    Args:
        step: Distance between consecutive evaluation points.

    Raises:
        FiniteDifferenceError: If ``step`` is not a real number or If ``step`` is zero.
    """
    if not (isinstance(step, int) or isinstance(step, float)):
        raise FiniteDifferenceError("step must be a real number")

    if step == 0:
        raise FiniteDifferenceError("step must be non-zero")


def forward_two_point(
    f_x: float | np.ndarray,
    f_x_plus_h: float | np.ndarray,
    step: int | float,
) -> float | np.ndarray:
    """Compute the first derivative using a two-point forward difference.

    The approximation is given by:

        f'(x) ≈ (f(x + h) - f(x)) / h

    This method has a truncation error of order O(h).

    Args:
        f_x(float or np.ndarray): Function value evaluated at x.
        f_x_plus_h(float or np.ndarray): Function value evaluated at x + h.
        step: Finite-difference step size h.

    Returns:
        Approximation of the first derivative at x.

    Raises:
        FiniteDifferenceError: If ``step`` is not a real number or If ``step`` is zero.
    """
    _validate_step(step)

    return (f_x_plus_h - f_x) / step

def central_two_point(
    f_x_plus_h: float | np.ndarray,
    f_x_minus_h: float | np.ndarray,
    step: int | float,
) -> float | np.ndarray:
    """Compute the first derivative using a two-point central difference.

    The approximation is given by:

        f'(x) ≈ (f(x + h) - f(x - h)) / (2h)

    This method has a truncation error of order O(h²).

    Args:
        f_x_plus_h(float or np.ndarray): Function value evaluated at x + h.
        f_x_minus_h(float or np.ndarray): Function value evaluated at x - h.
        step(int or float): Finite-difference step size h.

    Returns:
        Approximation of the first derivative at x.

    Raises:
        FiniteDifferenceError: If ``step`` is not a real number or If ``step`` is zero.
    """
    _validate_step(step)

    return (f_x_plus_h - f_x_minus_h) / (2.0 * step)


def forward_three_point(
    f_x: float | np.ndarray,
    f_x_plus_h: float | np.ndarray,
    f_x_plus_2h: float | np.ndarray,
    step: float | int,
) -> float | np.ndarray:
    """Compute the first derivative using a three-point forward difference.

    The approximation is given by:

        f'(x) ≈ (-3f(x) + 4f(x + h) - f(x + 2h)) / (2h)

    This method has a truncation error of order O(h²).

    Args:
        f_x(float or np.ndarray): Function value evaluated at x.
        f_x_plus_h(float or np.ndarray): Function value evaluated at x + h.
        f_x_plus_2h(float or np.ndarray): Function value evaluated at x + 2h.
        step(float or int): Finite-difference step size h.

    Returns:
        Approximation of the first derivative at x.

    Raises:
        FiniteDifferenceError: If ``step`` is not a real number or If ``step`` is zero.
    """
    _validate_step(step)

    return (
        -3.0 * f_x
        + 4.0 * f_x_plus_h
        - f_x_plus_2h
    ) / (2.0 * step)


def backward_three_point(
    f_x: float | np.ndarray,
    f_x_minus_h: float | np.ndarray,
    f_x_minus_2h: float | np.ndarray,
    step: float | int,
) -> float | np.ndarray:
    """Compute the first derivative using a three-point backward difference.

    The approximation is given by:

        f'(x) ≈ (3f(x) - 4f(x - h) + f(x - 2h)) / (2h)

    This method has a truncation error of order O(h²).

    Args:
        f_x(float or np.ndarray): Function value evaluated at x.
        f_x_minus_h(float or np.ndarray): Function value evaluated at x - h.
        f_x_minus_2h(float or np.ndarray): Function value evaluated at x - 2h.
        step(float or int): Finite-difference step size h.

    Returns:
        Approximation of the first derivative at x.

    Raises:
        FiniteDifferenceError: If ``step`` is not a real number or If ``step`` is zero.
    """
    _validate_step(step)

    return (
        3.0 * f_x
        - 4.0 * f_x_minus_h
        + f_x_minus_2h
    ) / (2.0 * step)


def central_five_point(
    f_x_minus_2h: float | np.ndarray,
    f_x_minus_h: float | np.ndarray,
    f_x_plus_h: float | np.ndarray,
    f_x_plus_2h: float | np.ndarray,
    step: float | int,
) -> float | np.ndarray:
    """Compute the first derivative using a five-point central difference.

    The approximation is given by:

        f'(x) ≈ (
            f(x - 2h)
            - 8f(x - h)
            + 8f(x + h)
            - f(x + 2h)
        ) / (12h)

    This method has a truncation error of order O(h⁴).

    Args:
        f_x_minus_2h(float or np.ndarray): Function value evaluated at x - 2h.
        f_x_minus_h(float or np.ndarray): Function value evaluated at x - h.
        f_x_plus_h(float or np.ndarray): Function value evaluated at x + h.
        f_x_plus_2h(float or np.ndarray): Function value evaluated at x + 2h.
        step(float or int): Finite-difference step size h.

    Returns:
        Approximation of the first derivative at x.

    Raises:
        FiniteDifferenceError: If ``step`` is not a real number or If ``step`` is zero.
    """
    _validate_step(step)

    return (
        f_x_minus_2h
        - 8.0 * f_x_minus_h
        + 8.0 * f_x_plus_h
        - f_x_plus_2h
    ) / (12.0 * step)