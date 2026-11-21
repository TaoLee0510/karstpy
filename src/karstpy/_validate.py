"""Shared validation for count matrices and scalar parameters."""

from numbers import Real

import numpy as np
import pandas as pd


def scalar(value, name, lower=-np.inf, upper=np.inf, *, lower_open=False, upper_open=False):
    if not isinstance(value, Real) or isinstance(value, (bool, np.bool_)) or not np.isfinite(value):
        raise ValueError(f"`{name}` must be one finite number.")
    if not lower <= value <= upper or (lower_open and value == lower) or (upper_open and value == upper):
        left, right = "(" if lower_open else "[", ")" if upper_open else "]"
        raise ValueError(f"`{name}` must be in {left}{lower}, {upper}{right}.")
    return float(value)


def integer(value, name, lower=-2147483647, upper=2147483647):
    value = scalar(value, name, lower, upper)
    if value != round(value):
        raise ValueError(f"`{name}` must be an integer.")
    return int(value)


def validate_beta(beta):
    return scalar(beta, "beta", 0, 1, lower_open=True, upper_open=True)


def count_frame(counts, labels=None, *, min_columns, drop_empty=False):
    if isinstance(counts, pd.DataFrame):
        if labels is not None:
            raise ValueError("Use the DataFrame index for karyotype labels.")
        labels, columns = counts.index, counts.columns
        counts = counts.to_numpy()
    else:
        counts, columns = np.asarray(counts), None
    if counts.ndim != 2 or counts.shape[1] < min_columns:
        raise ValueError(f"`counts` needs at least {min_columns} columns, one per time point.")
    if labels is None:
        raise ValueError("`counts` must have karyotype row labels.")
    labels = list(labels)
    if len(labels) != counts.shape[0] or not all(isinstance(label, str) and label for label in labels):
        raise ValueError("`counts` must have one non-empty karyotype label per row.")
    try:
        values = counts.astype(float)
    except (ValueError, TypeError) as error:
        raise ValueError("`counts` must contain finite non-negative integers.") from error
    rounded = np.rint(values)
    if not np.all(np.isfinite(values) & (values >= 0) & (abs(values - rounded) <= 1e-8)):
        raise ValueError("`counts` must contain finite non-negative integers.")
    if np.any(rounded > np.iinfo(np.int32).max):
        raise ValueError("`counts` exceeds the integer range.")
    if columns is None:
        columns = [f"t{i}" for i in range(values.shape[1])]
    if drop_empty:
        keep = np.any(rounded != 0, axis=1)
        if not keep.any():
            raise ValueError("At least one karyotype must have non-zero counts.")
        rounded = rounded[keep]
        labels = np.asarray(labels)[keep]
    return pd.DataFrame(rounded.astype(np.int64), index=labels, columns=columns)
