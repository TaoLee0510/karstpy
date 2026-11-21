"""Prepare karyotype counts and their observation times."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ._validate import count_frame, scalar, validate_beta
from .karyotypes import format_karyotypes, parse_karyotypes


@dataclass
class KarstData:
    """Aligned counts, karyotype labels, and elapsed times."""

    counts: pd.DataFrame
    labels: list[str]
    dt: float
    beta: float
    metadata: dict
    times: np.ndarray | None = None

    def __repr__(self):
        return f"KarstData({len(self.labels)} karyotypes, {self.counts.shape[1]} time points, dt={self.dt}, beta={self.beta})"


def resolve_times(n_time, dt, times):
    if times is not None:
        if dt is not None:
            raise ValueError("Give `dt` or `times`, not both.")
        values = np.asarray(times)
        if values.shape != (n_time,) or values.dtype.kind not in "iuf" or not np.all(np.isfinite(values)):
            raise ValueError(f"`times` must give {n_time} finite times, one per time point.")
        values = values.astype(float)
        if np.any(np.diff(values) <= 0):
            raise ValueError("`times` must be strictly increasing.")
        return values - values[0]
    if n_time == 2:
        interval = 1.0 if dt is None else scalar(dt, "dt", np.finfo(float).eps)
        return np.array([0.0, interval])
    intervals = np.asarray(dt)
    if (intervals.shape != (n_time - 1,) or intervals.dtype.kind not in "iuf"
            or not np.all(np.isfinite(intervals) & (intervals > 0))):
        raise ValueError(f"With {n_time} time points give `times` or {n_time - 1} positive intervals.")
    return np.r_[0.0, np.cumsum(intervals, dtype=float)]


def prepare_karst_data(counts, dt=None, beta=0.00005, metadata=None, times=None, *, labels=None):
    """Prepare a count matrix with one labelled row per karyotype and one column per time point."""
    frame = count_frame(counts, labels, min_columns=2, drop_empty=True)
    n_time = frame.shape[1]
    elapsed = resolve_times(n_time, dt, times)
    beta = validate_beta(beta)
    if metadata is not None and not isinstance(metadata, dict):
        raise ValueError("`metadata` must be a dictionary.")
    canonical = format_karyotypes(parse_karyotypes(frame.index, unique=False))
    if len(set(canonical)) != len(canonical):
        for label in dict.fromkeys(canonical):
            matches = [old for old, new in zip(frame.index, canonical) if new == label]
            if len(matches) > 1:
                raise ValueError(f"Karyotype labels {' and '.join(map(repr, matches))} are the same karyotype.")
    frame.index = canonical
    frame.columns = [f"t{i}" for i in range(n_time)]
    return KarstData(frame, canonical, elapsed[-1], beta, dict(metadata or {}), elapsed if n_time > 2 else None)


def karst_data_times(data):
    """Return each time point measured from the first observation."""
    if not isinstance(data, KarstData):
        raise ValueError("`data` must be a KarstData object.")
    return np.array([0.0, data.dt]) if data.times is None else data.times.copy()
