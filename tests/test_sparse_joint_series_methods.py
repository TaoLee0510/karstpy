import numpy as np
import pytest

from karstpy import as_alfak, fit_sparse_joint


@pytest.mark.parametrize("method", ["anchored", "pairwise", "state_space"])
def test_two_time_points_use_the_same_fit(fit, counts, method):
    another = fit_sparse_joint(counts, dt=3, history=9, series=method)
    np.testing.assert_array_equal(fit.summary.relative_fitness, another.summary.relative_fitness)


@pytest.mark.parametrize("method,objective", [("pairwise", 116.93968496574617), ("state_space", 74.33207179004155)])
def test_series_method_matches_r_on_example(series_counts, series_times, method, objective):
    fit = fit_sparse_joint(series_counts, times=series_times, history=9, series=method)
    assert fit.diagnostics["series"] == method
    assert fit.diagnostics["posterior_objective"] == pytest.approx(objective, rel=1e-6)
    assert fit.summary[["expected_t0", "expected_t3"]].notna().any().all()
    assert as_alfak(fit, minobs=10, nboot=5, xval=False).landscape_posterior_samples.shape[1] == 5
