"""Diagnostics for sparse karyotype counts."""

import numpy as np
import pandas as pd

from ._validate import count_frame, scalar
from .data import KarstData
from .karyotypes import parse_karyotypes


def _nearest_median(karyotypes, cells):
    nearest = np.zeros(len(cells))
    singleton = np.flatnonzero(cells == 1)
    for start in range(0, len(singleton), 200):
        rows = singleton[start:start + 200]
        distance = np.zeros((len(rows), len(cells)))
        for chromosome in karyotypes.T:
            distance += abs(chromosome[rows, None] - chromosome[None, :])
        distance[np.arange(len(rows)), rows] = np.inf
        nearest[rows] = distance.min(axis=1)
    order = np.argsort(nearest, kind="stable")
    cumulative = np.cumsum(cells[order])
    middle = np.searchsorted(cumulative, [(cumulative[-1] - 1) // 2, cumulative[-1] // 2], side="right")
    return nearest[order[middle]].mean()


def check_sparse_joint_counts(counts, major_count=3, *, labels=None):
    """Summarise repeated, shared and nullisomic karyotypes at each time point."""
    if isinstance(counts, KarstData):
        counts = counts.counts
    frame = count_frame(counts, labels, min_columns=1)
    major_count = scalar(major_count, "major_count", 1)
    karyotypes = parse_karyotypes(frame.index).astype(np.int64)
    nullisomic = (karyotypes == 0).any(axis=1)
    values = frame.to_numpy()
    observed_at_multiple_times = (values > 0).sum(axis=1) > 1
    rows = []
    for name, cells in zip(frame.columns, values.T):
        seen = cells > 0
        total = cells.sum()
        rows.append({
            "sample": name,
            "cells": total,
            "karyotypes": seen.sum(),
            "repeated_fraction": cells[cells >= major_count].sum() / total if total else np.nan,
            "nearest_distance": _nearest_median(karyotypes[seen], cells[seen]) if total > 1 else np.nan,
            "nullisomy_fraction": cells[nullisomic].sum() / total if total else np.nan,
            "shared": (seen & observed_at_multiple_times).sum() if values.shape[1] > 1 else np.nan,
        })
    return pd.DataFrame(rows)
