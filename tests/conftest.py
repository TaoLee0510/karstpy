import numpy as np
import pandas as pd
import pytest

from karstpy import fit_sparse_joint


@pytest.fixture
def series_counts():
    return pd.DataFrame(
        [[180, 240, 190, 175], [36, 31, 51, 59], [2, 0, 1, 0], [6, 3, 4, 2]],
        index=["2.2", "2.3", "1.2", "3.2"],
    )


@pytest.fixture
def counts(series_counts):
    return series_counts.iloc[:, :2]


@pytest.fixture
def series_times():
    return np.array([0, 3, 6, 9])


@pytest.fixture
def fit(counts):
    return fit_sparse_joint(counts, dt=3, history=9)
