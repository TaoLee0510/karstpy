import numpy as np
import pandas as pd
import pytest

from karstpy import KarstDomainCoverageError, check_sparse_joint_counts, fit_sparse_joint
from karstpy._problem import one_step_neighbours


def test_domain_has_majors_and_all_neighbours(fit):
    majors = fit.summary.loc[fit.summary.major, "karyotype"].tolist()
    expected = set(majors) | set(one_step_neighbours(fit.trend["major"], 1, 4))
    assert set(fit.summary.karyotype) == expected
    assert fit.summary.karyotype.iloc[:len(majors)].tolist() == majors


def test_predict_inside_outside_and_with_repeated_labels(fit):
    inside = fit.summary.karyotype.iloc[0]
    result = fit.predict([inside, "4.4", "02.2"])
    np.testing.assert_array_equal(result.in_domain, [True, False, True])
    assert result.relative_fitness.iloc[0] == result.relative_fitness.iloc[2]
    assert np.isfinite(result.relative_sd).all()
    assert (result.relative_sd >= result.relative_sd_model).all()


def test_typed_error_includes_coverage():
    counts = pd.DataFrame([[30, 0], [0, 1]], index=["2.2", "6.6"])
    with pytest.raises(KarstDomainCoverageError) as caught:
        fit_sparse_joint(counts, history=9)
    assert caught.value.coverage.cells_in_domain.tolist() == [30, 0]


def test_count_check_can_use_a_single_time_point():
    counts = pd.DataFrame([3, 1, 1], index=["2.2", "2.3", "0.2"], columns=["t0"])
    row = check_sparse_joint_counts(counts).iloc[0]
    assert row.repeated_fraction == 0.6
    assert row.nullisomy_fraction == 0.2
    assert row.nearest_distance == 0.0
