import numpy as np
import pandas as pd
import pytest

from karstpy import KarstSamplesLeftOutWarning, fit_sparse_joint, karst_data_times, option_context, prepare_karst_data


def test_prepare_preserves_elapsed_times(series_counts):
    data = prepare_karst_data(series_counts, times=[8, 9, 13, 15])
    np.testing.assert_array_equal(karst_data_times(data), [0, 1, 5, 7])
    with pytest.raises(ValueError, match="not both"):
        prepare_karst_data(series_counts, dt=[1, 4, 2], times=[8, 9, 13, 15])


def test_series_template_agrees_for_two_time_points(counts, fit):
    with option_context(sparse_joint_template="series"):
        other = fit_sparse_joint(counts, dt=3, history=9)
    np.testing.assert_allclose(other.summary.relative_fitness, fit.summary.relative_fitness, atol=1e-7)


def test_intermediate_time_point_without_coverage_is_left_out():
    counts = pd.DataFrame([[100, 0, 89], [2, 0, 4], [0, 80, 0]], index=["2.2", "2.3", "6.6"])
    with pytest.warns(KarstSamplesLeftOutWarning):
        fit = fit_sparse_joint(counts, times=[0, 2, 4], history=4, max_cn=7)
    assert "count_t2" in fit.summary and "count_t1" not in fit.summary
    assert fit.input is not None
    assert fit.diagnostics["dropped_samples"].timepoint.tolist() == ["t1"]
    assert "6.6" in fit.dropped
