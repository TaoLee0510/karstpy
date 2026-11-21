import pandas as pd
import pytest

from karstpy import fit_karst, fit_sparse_joint, prepare_karst_data


def test_fit_wrapper_uses_prepared_beta(counts):
    data = prepare_karst_data(counts, dt=3, beta=0.002)
    assert fit_karst(data, sparse_joint_control={"history": 9}).data.beta == data.beta
    with pytest.raises(ValueError, match="either directly"):
        fit_karst(data, beta=0.002, sparse_joint_control={"beta": 0.002})


@pytest.mark.parametrize("beta", [0, 1, -0.1, float("nan")])
def test_beta_must_be_strictly_positive_and_below_one(counts, beta):
    with pytest.raises(ValueError, match="beta"):
        prepare_karst_data(counts, beta=beta)


def test_canonical_duplicate_rows_are_rejected():
    counts = pd.DataFrame([[1, 2], [5, 6]], index=["2.02.3", "2.2.3"])
    with pytest.raises(ValueError, match="same karyotype"):
        prepare_karst_data(counts)


@pytest.mark.parametrize("method", ["anchored", "pairwise", "state_space"])
def test_free_balance_level_is_refused_for_series(series_counts, series_times, method):
    with pytest.raises(ValueError, match="two time points"):
        fit_sparse_joint(series_counts, times=series_times, series=method, balance_level="free")
