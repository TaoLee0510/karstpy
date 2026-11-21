import numpy as np
import pandas as pd
import pytest

from karstpy import fit_sparse_joint
from karstpy._problem import build_problem
from karstpy.options import get_option, option_context


def test_scope_restores_the_fit_options():
    assert not get_option("choose_accepted")
    with option_context(choose_accepted=True):
        assert get_option("choose_accepted")
        with option_context(choose_accepted=False):
            assert not get_option("choose_accepted")
        assert get_option("choose_accepted")
    assert not get_option("choose_accepted")


def test_marginal_gradient_matches_finite_difference(fit):
    problem = build_problem(fit.data, fit.settings)
    objective = problem.objective(fit.settings)
    x = objective.par
    value = objective.fn(x)
    best = objective.last_par_best.copy()
    gradient = objective.gr(x)
    shift = np.eye(len(x)) * 2e-4
    finite = [(objective.fn(x + d) - objective.fn(x - d)) / 4e-4 for d in shift]
    np.testing.assert_allclose(gradient, finite, rtol=1e-4, atol=2e-5)
    assert objective.value_best <= value
    np.testing.assert_array_equal(best[objective.fixed], x)


def test_fitted_scale_checks_and_spectrum_are_reported(fit):
    check = fit.diagnostics["optimisation"]
    assert check["optimizer"] == "L-BFGS-B"
    assert fit.diagnostics["convergence"] == 0
    assert check["min_eigen"] >= 0.002
    assert check["projected_gradient"] < 0.05


@pytest.mark.parametrize(("method", "value", "gradient"), [
    ("anchored", 133.5477295081503, [-0.0856283784208, 0.00481515464913, -0.0501339115740]),
    ("state_space", 140.52278001796137, [-0.0791664705490, 0.00417234749261, -0.0459085293263, 0.0482792763600]),
])
def test_laplace_with_nearly_equal_growth_rates(method, value, gradient):
    # R 4.4.2 / TMB 1.9.19 with only the log-phi1 stability fix.
    # In the legacy formula this fixture amplifies roundoff in the third derivatives.
    counts = pd.DataFrame(
        [[94, 124, 112, 128], [12, 14, 18, 14], [28, 32, 36, 33], [15, 27, 21, 28],
         [29, 31, 31, 38], [16, 20, 14, 16], [26, 32, 48, 62], [0, 0, 0, 1]],
        index=["2.2.2", "1.2.2", "3.2.2", "2.1.2", "2.3.2", "2.2.1", "2.2.3", "1.3.2"],
    )
    fit = fit_sparse_joint(counts, times=[0, 0.65, 1.35, 2], beta=0.001, max_cn=4,
                           history=2, series=method, intermediate_weight=1)
    objective = build_problem(fit.data, fit.settings).objective(fit.settings)
    for offset in (0, 0.01, -0.015):
        fixed = objective.par + offset
        actual_value = objective.fn(fixed)
        actual_gradient = objective.gr(fixed)
    assert actual_value == pytest.approx(value, abs=1e-8, rel=0)
    np.testing.assert_allclose(actual_gradient, gradient, rtol=0, atol=1e-8)


def test_integral_derivatives_when_growth_gap_crosses_series_cutoff(fit):
    problem = build_problem(fit.data, fit.settings)
    data = problem.tmb_data
    # Isolate one directed edge and reuse the tape while changing its rate gap.
    one_edge = {name: data[name][:1] for name in ("edge_from", "edge_to", "edge_log_prob")}
    objective = problem.objective(fit.settings, {**data, **one_edge, "balance": 0})
    source, target = data["edge_from"][0], data["edge_to"][0]
    net = 2 * np.exp(data["log_self"]) - 1
    source_index, target_index = objective.slices["g"][[source, target]]
    step = np.zeros_like(objective.initial)
    step[source_index] = 1e-4
    weights = (
        np.isin(objective.row, [source_index, target_index])
        & np.isin(objective.column, [source_index, target_index])
    ).astype(float)
    assert weights.sum() == 4

    for gap in (-2.0, -0.5, 0.0, 0.5, 2.0):
        point = objective.initial.copy()
        point[target_index] = 0.0
        point[source_index] = np.log((net[target] + gap / data["dt"]) / net[source])
        plus, minus = point + step, point - step
        width = 2 * step[source_index]
        numeric_grad = (objective.value(plus) - objective.value(minus)) / width
        numeric_hess = (objective.gradient(plus) - objective.gradient(minus)) / width
        numeric_third = weights @ (objective.tape.hessian_values(plus) - objective.tape.hessian_values(minus)) / width

        assert objective.gradient(point)[source_index] == pytest.approx(numeric_grad, rel=2e-7, abs=1e-6)
        np.testing.assert_allclose(objective.hessian(point, random=False)[:, source_index].toarray().ravel(),
                                   numeric_hess, rtol=2e-7, atol=1e-6)
        third = objective.tape.hessian_reverse(point, weights)
        assert third[source_index] == pytest.approx(numeric_third, rel=2e-7, abs=1e-6)
