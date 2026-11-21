"""Convert karyotypes between dot-separated labels and copy numbers."""

import numpy as np

from . import _core


def parse_karyotypes(labels, *, unique=True):
    """Parse chromosome copy numbers as an integer matrix, one row per label."""
    if isinstance(labels, str):
        labels = [labels]
    return _core.parse_karyotypes(list(labels), unique)


def format_karyotypes(karyotypes):
    """Format a copy-number matrix as a list of dot-separated labels."""
    matrix = np.asarray(karyotypes)
    if matrix.ndim == 1:
        matrix = matrix[None, :]
    if matrix.ndim != 2 or matrix.dtype.kind not in "ifu":
        raise ValueError("Karyotypes must be a numeric vector or matrix.")
    limits = np.iinfo(np.int32)
    if not np.all(np.isfinite(matrix) & (matrix == np.rint(matrix)) & (matrix >= limits.min) & (matrix <= limits.max)):
        raise ValueError("Copy numbers must be finite integers in the supported range.")
    return _core.format_karyotypes(matrix.astype(np.int32, copy=False))
