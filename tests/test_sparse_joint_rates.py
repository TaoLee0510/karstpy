import numpy as np
import pytest

from karstpy import fit_sparse_joint_rates


@pytest.mark.parametrize("parameter,interval", [("birth_rate", (0.6, 1.2)), ("beta", (5e-5, 8e-4))])
def test_rate_range_covers_all_intervals(counts, parameter, interval):
    fit = fit_sparse_joint_rates(counts, dt=3, history=9, **{parameter + "_range": interval})
    assert fit.parameter == parameter
    np.testing.assert_allclose(fit.rates, np.geomspace(*interval, 3))
    labels = fit.fits[0].summary.karyotype
    out = fit.predict(labels)
    for variant in fit.fits:
        prediction = variant.predict(labels)
        np.testing.assert_array_less(out.relative_low - 1e-12, prediction.relative_fitness - 1.9599639845400536 * prediction.relative_sd)
        np.testing.assert_array_less(prediction.relative_fitness + 1.9599639845400536 * prediction.relative_sd, out.relative_high + 1e-12)


@pytest.mark.parametrize("kwargs", [{}, {"beta_range": (0, 1e-4)}, {"beta_range": (1e-3, 1e-4)},
                                    {"beta_range": (1e-4, 1e-3), "beta": 1e-4}])
def test_invalid_rate_ranges(counts, kwargs):
    with pytest.raises(ValueError):
        fit_sparse_joint_rates(counts, **kwargs)
