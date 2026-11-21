"""Isolated streams matching R 4.4.2's default seeded draws."""

import numpy as np

from ._core import RandomStream
from ._validate import integer


def random_stream(seed):
    """Create a local stream without modifying another stream's state."""
    return RandomStream(integer(seed, "seed", 0))


def normal_array(stream, shape):
    """Fill an array in the same order as R's matrix constructor."""
    size = int(np.prod(shape))
    return np.fromiter((stream.normal() for _ in range(size)), dtype=float, count=size).reshape(shape, order="F")
