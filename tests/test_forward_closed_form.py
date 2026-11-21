import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm

from karstpy import _core
from karstpy.forward import closed_form_forward, exact_forward_reference, log_phi1, log_phi2, sparse_division_kernel


def test_exact_kernel_retains_full_daughter_loss():
    beta = 0.012
    k = sparse_division_kernel([[2, 2], [1, 2], [2, 1], [3, 2]], beta)
    assert list(zip(k.source, k.target)) == [(0, 1), (3, 0), (0, 3), (1, 0), (0, 2), (2, 0)]
    np.testing.assert_allclose(np.exp(k.log_prob[0]), _core.pij(2, 1, beta) * _core.pij(2, 2, beta))
    np.testing.assert_allclose(np.exp(k.log_self[0]), _core.pij(2, 2, beta)**2)


def test_divided_differences_match_integrals():
    dt, a, b, c = 2.2, 0.3, 0.1, -0.4
    np.testing.assert_allclose(np.exp(log_phi1(a, b, dt)), quad(lambda t: np.exp(a * (dt - t) + b * t), 0, dt)[0])
    exact = quad(lambda t: quad(lambda s: np.exp(a * (dt - t) + b * (t - s) + c * s), 0, t)[0], 0, dt)[0]
    np.testing.assert_allclose(np.exp(log_phi2(a, b, c, dt)), exact)
    np.testing.assert_allclose(np.exp(log_phi2(a, a, a, dt)), np.exp(a * dt) * dt**2 / 2)


def test_uniformisation_and_dyson_against_dense_exponential():
    k = sparse_division_kernel([[2, 2], [2, 3], [3, 3]], 5e-5)
    birth = np.array([1, 0.9, 1.1])
    n0 = np.array([100, 15, 4])
    net = birth * (2 * np.exp(k.log_self) - 1)
    a = np.diag(net)
    a[k.target, k.source] = 2 * birth[k.source] * np.exp(k.log_prob)
    expected = np.log(expm(a * 4) @ n0)
    reference = exact_forward_reference(np.log(n0), np.log(birth), k, 4)
    first = closed_form_forward(np.log(n0), np.log(birth), k, 4)
    second = closed_form_forward(np.log(n0), np.log(birth), k, 4, order=2)
    np.testing.assert_allclose(reference, expected, rtol=1e-12)
    assert np.max(abs(second - expected)) < np.max(abs(first - expected))
