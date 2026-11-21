from dataclasses import replace

import numpy as np

from karstpy.sparse_joint import rebuild_posterior


def test_rebuild_uses_exact_stored_mode(fit):
    solver = rebuild_posterior(fit)
    assert solver.stored
    assert solver.value == fit.diagnostics["posterior_objective"]
    np.testing.assert_array_equal(solver.mode["g"], fit.posterior["g"])
    np.testing.assert_allclose(solver.expected, fit.summary[["expected_t0", "expected_t1"]])


def test_rebuild_does_not_move_an_inexact_stored_mode(fit):
    altered = replace(fit, posterior=dict(fit.posterior, g=fit.posterior["g"] + 0.0001))
    solver = rebuild_posterior(altered)
    np.testing.assert_array_equal(solver.mode["g"], altered.posterior["g"])
    assert np.max(abs(solver.objective.gradient(solver.point))) > 0.001
