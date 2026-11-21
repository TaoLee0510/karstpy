import pytest

import karstpy


def test_current_entry_points():
    for name in ("fit_karst", "fit_sparse_joint", "fit_sparse_joint_rates", "as_alfak", "alfak_compat"):
        assert callable(getattr(karstpy, name))
        assert name in karstpy.__all__


@pytest.mark.parametrize("name", ["fit_local_posterior", "fit_extrapolation_layer", "fit_graph_posterior", "plot_karst"])
def test_removed_pipeline_exports(name):
    assert not hasattr(karstpy, name)
