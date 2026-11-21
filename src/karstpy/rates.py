"""Sensitivity of fitted landscapes to externally supplied rates."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ._validate import integer, scalar
from .sparse_joint import fit_sparse_joint


@dataclass
class KarstRateFit:
    rates: np.ndarray
    fits: list
    parameter: str
    point: int
    central: int

    def __repr__(self):
        return f"KarstRateFit({self.parameter}={self.rates.tolist()}, point={self.rates[self.point]:g})"

    def predict(self, karyotypes):
        """Return point predictions and intervals spanning the given rates."""
        predictions = [fit.predict(karyotypes) for fit in self.fits]
        means = np.column_stack([p.relative_fitness for p in predictions])
        deviations = np.column_stack([p.relative_sd for p in predictions])
        critical = 1.9599639845400536
        low = (means - critical * deviations).min(axis=1)
        high = (means + critical * deviations).max(axis=1)
        per_division = means / [fit.diagnostics["birth_rate"] for fit in self.fits]
        centre = means[:, self.point]
        result = pd.DataFrame({
            "karyotype": predictions[0].karyotype,
            "relative_fitness": centre, "relative_low": low, "relative_high": high,
            "relative_sd": np.maximum(high - centre, centre - low) / critical,
            "per_division": per_division[:, self.point],
            "per_division_low": per_division.min(axis=1), "per_division_high": per_division.max(axis=1),
        })
        agreement = np.nan
        if len(result) > 2:
            correlation = pd.DataFrame(means).corr(method="spearman").to_numpy()
            pair = correlation[np.triu_indices_from(correlation, k=1)]
            if np.isfinite(pair).any():
                agreement = np.nanmin(pair)
        result.attrs["rank_agreement"] = agreement
        return result


def fit_sparse_joint_rates(counts, *, birth_rate_range=None, beta_range=None, n_rates=3, point="lowest", **kwargs):
    """Fit one landscape per log-spaced rate within a supplied rate range."""
    if (birth_rate_range is None) == (beta_range is None):
        raise ValueError("Give one of `birth_rate_range` and `beta_range`.")
    parameter = "birth_rate" if beta_range is None else "beta"
    interval = np.asarray(birth_rate_range if beta_range is None else beta_range)
    if interval.shape != (2,):
        raise ValueError(f"`{parameter}_range` must give two rates.")
    lower = scalar(interval[0], f"{parameter}_range lower", 0, lower_open=True)
    upper = scalar(interval[1], f"{parameter}_range upper", lower, lower_open=True)
    if parameter == "beta" and upper >= 1:
        raise ValueError("`beta_range` must lie in (0, 1).")
    if parameter in kwargs:
        raise ValueError(f"Give `{parameter}_range` without a separate `{parameter}`.")
    if point not in ("lowest", "centre"):
        raise ValueError("`point` must be 'lowest' or 'centre'.")
    rates = np.geomspace(lower, upper, integer(n_rates, "n_rates", 2))
    fits = [fit_sparse_joint(counts, **{parameter: rate}, **kwargs) for rate in rates]
    central = (len(fits) - 1) // 2
    return KarstRateFit(rates, fits, parameter, 0 if point == "lowest" else central, central)
