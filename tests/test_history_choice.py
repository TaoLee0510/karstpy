import numpy as np

from karstpy import fit_sparse_joint


def test_history_selects_finite_posterior_objective(counts):
    fit = fit_sparse_joint(counts, dt=3, history=[9, 3, 15])
    profile = fit.diagnostics["history_profile"]
    valid = profile[np.isfinite(profile.objective)]
    best = valid.loc[valid.posterior_objective.idxmin()]
    assert fit.diagnostics["history"] == best.history
    assert fit.diagnostics["posterior_objective"] == best.posterior_objective


def test_history_is_irrelevant_for_free_initial_state(counts):
    a = fit_sparse_joint(counts, dt=3, history=[9, 3], initial_state="free")
    b = fit_sparse_joint(counts, dt=3, history=3, initial_state="free")
    np.testing.assert_array_equal(a.summary.relative_fitness, b.summary.relative_fitness)
