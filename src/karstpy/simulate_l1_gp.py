"""Lazy sine-wave fitness landscapes."""

from dataclasses import dataclass, field
import math

import numpy as np

from ._validate import integer, scalar
from .karyotypes import format_karyotypes, parse_karyotypes
from .rng import random_stream


@dataclass
class GRFLandscape:
    n_chr: int
    n_centroids: int
    lambda_: float
    centroids: np.ndarray
    founder: np.ndarray
    founder_label: str
    founder_fitness: float
    founder_raw: float
    scale: float
    seed: int
    centroid_min: int
    centroid_max: int
    min_cn: int
    max_cn: int
    family: str = "paper_sine_grf"
    covariance: None = None
    formula: str = "founder_fitness + scale * (sum_i sin(||k-r_i||_2/lambda)/(pi*sqrt(M)) - raw(founder))"
    _cache: tuple | None = field(default=None, repr=False, compare=False)

    def _signature(self):
        return self.lambda_, self.scale, self.founder_fitness, self.founder_raw, self.centroids.shape, self.centroids.tobytes()


def _raw_grf(karyotypes, centroids, wavelength):
    values = np.zeros(len(karyotypes))
    for centroid in centroids:
        radius = ((karyotypes.astype(float) - centroid) ** 2).sum(axis=1) ** 0.5
        values += np.fromiter((math.sin(x / wavelength) for x in radius), float, count=len(values))
    return values / (math.pi * math.sqrt(len(centroids)))


def simulate_grf_landscape(
    n_chr=22, n_centroids=30, lambda_=0.8, founder=None, founder_cn=2,
    founder_fitness=1, scale=1, seed=1, centroid_min=-10, centroid_max=20,
    min_cn=1, max_cn=8, cache=True,
):
    """Generate a deterministic oracle without enumerating the full karyotype lattice."""
    n_chr = integer(n_chr, "n_chr", 1)
    n_centroids = integer(n_centroids, "n_centroids", 1)
    lambda_ = scalar(lambda_, "lambda_", 0, lower_open=True)
    founder_cn = integer(founder_cn, "founder_cn")
    founder_fitness = scalar(founder_fitness, "founder_fitness")
    scale = scalar(scale, "scale")
    seed = integer(seed, "seed", 0)
    centroid_min = integer(centroid_min, "centroid_min")
    centroid_max = integer(centroid_max, "centroid_max", centroid_min)
    min_cn = integer(min_cn, "min_cn")
    max_cn = integer(max_cn, "max_cn", min_cn)
    if founder is None:
        founder = np.full(n_chr, founder_cn)
    else:
        founder = np.asarray(founder)
        if founder.shape != (n_chr,) or founder.dtype.kind not in "iuf" or not np.all(founder == np.rint(founder)):
            raise ValueError(f"`founder` must be an integer vector of length {n_chr}.")
    if not np.all((founder >= min_cn) & (founder <= max_cn)):
        raise ValueError("`founder` must lie inside [min_cn, max_cn].")
    if not isinstance(cache, (bool, np.bool_)):
        raise ValueError("`cache` must be Boolean.")
    founder = founder.astype(np.int32)
    stream = random_stream(seed)
    indices = stream.sample(centroid_max - centroid_min + 1, n_centroids * n_chr, True)
    centroids = np.asarray(indices).reshape(n_centroids, n_chr, order="F") + centroid_min
    landscape = GRFLandscape(
        n_chr, n_centroids, lambda_, centroids, founder, format_karyotypes(founder)[0],
        founder_fitness, _raw_grf(founder[None, :], centroids, lambda_)[0],
        scale, seed, centroid_min, centroid_max, min_cn, max_cn,
    )
    if cache:
        landscape._cache = landscape._signature(), {}
    return landscape


def simulate_l1_gp_landscape(
    n_chr=22, min_cn=1, max_cn=8, diploid_cn=2, diploid_fitness=1,
    lower=None, upper=None, ell=None, seed=1, use_full_range=None,
    include_table=None, n_centroids=30, lambda_=None, founder=None,
    founder_cn=None, founder_fitness=None, scale=1, centroid_min=-10,
    centroid_max=20, cache=True,
):
    """Construct a lazy GRF landscape using the older entry point."""
    return simulate_grf_landscape(
        n_chr=n_chr, n_centroids=n_centroids, lambda_=lambda_ if lambda_ is not None else (0.8 if ell is None else ell),
        founder=founder, founder_cn=diploid_cn if founder_cn is None else founder_cn,
        founder_fitness=diploid_fitness if founder_fitness is None else founder_fitness,
        scale=scale, seed=seed, centroid_min=centroid_min, centroid_max=centroid_max,
        min_cn=min_cn, max_cn=max_cn, cache=cache,
    )


def predict_landscape_fitness(landscape, karyotypes, use_cache=True):
    """Evaluate a GRF landscape on labels or integer copy-number vectors."""
    if not isinstance(landscape, GRFLandscape):
        raise ValueError("`landscape` must be a GRFLandscape object.")
    if not isinstance(use_cache, (bool, np.bool_)):
        raise ValueError("`use_cache` must be Boolean.")
    values = np.asarray(karyotypes)
    if values.ndim == 0 and isinstance(values.item(), str):
        values = values.reshape(1)
    if values.ndim == 1 and values.size and all(isinstance(item, str) for item in values):
        values = parse_karyotypes(values.tolist(), unique=False)
    else:
        values = np.atleast_2d(values)
    if values.ndim != 2 or values.shape[1] != landscape.n_chr:
        raise ValueError(f"`karyotypes` must have {landscape.n_chr} columns.")
    labels = format_karyotypes(values)
    cache = landscape._cache
    stored = cache[1] if use_cache and cache is not None and cache[0] == landscape._signature() else {}
    missing = [i for i, label in enumerate(labels) if label not in stored]
    if missing:
        raw = _raw_grf(values[missing], landscape.centroids, landscape.lambda_)
        result = landscape.founder_fitness + landscape.scale * (raw - landscape.founder_raw)
        stored.update((labels[i], x) for i, x in zip(missing, result))
    return np.array([stored[label] for label in labels])
