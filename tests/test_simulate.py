import numpy as np
import pandas as pd
import pytest

from karstpy import _core, simulate_grf_landscape, simulate_sparse_counts


def test_seeded_stream_uses_r_mersenne_twister():
    # R 4.4.2 set.seed(1); writeBin(runif(3), ...)
    stream = _core.RandomStream(1)
    np.testing.assert_array_equal([stream.uniform() for _ in range(3)],
                                  [0.2655086631421, 0.37212389963679016, 0.5728533633518964])


@pytest.mark.parametrize("dispersion", [0, 35])
def test_simulation_is_reproducible(dispersion):
    land = simulate_grf_landscape(n_chr=3)
    a = simulate_sparse_counts(land, beta=0.008, seed=19, n0=100, n1=200, overdispersion=dispersion)
    b = simulate_sparse_counts(land, beta=0.008, seed=19, n0=100, n1=200, overdispersion=dispersion)
    pd.testing.assert_frame_equal(a["counts"], b["counts"])
    np.testing.assert_array_equal(a["counts"].sum(axis=0), [100, 200])
    assert len(a["truth_observed"]) == len(a["counts"])
    assert (a["counts"].sum(axis=1) > 0).all()
