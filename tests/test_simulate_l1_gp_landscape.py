from dataclasses import replace

import numpy as np
import pandas as pd

from karstpy import predict_landscape_fitness, simulate_grf_landscape, simulate_l1_gp_landscape


def test_lazy_landscape_anchors_founder():
    landscape = simulate_grf_landscape(n_chr=22, n_centroids=12, founder_fitness=1.5, seed=112)
    assert landscape.centroids.shape == (12, 22)
    assert predict_landscape_fitness(landscape, landscape.founder_label)[0] == 1.5
    labels = pd.Series([landscape.founder_label], dtype=object)
    np.testing.assert_array_equal(predict_landscape_fitness(landscape, labels), [1.5])


def test_degenerate_centroid_range_and_legacy_interface_agree():
    legacy = simulate_l1_gp_landscape(n_chr=4, ell=1.2, centroid_min=2, centroid_max=2)
    oracle = simulate_grf_landscape(n_chr=4, lambda_=1.2, centroid_min=2, centroid_max=2)
    assert np.all(legacy.centroids == 2)
    np.testing.assert_array_equal(predict_landscape_fitness(legacy, ["1.2.2.2"]),
                                  predict_landscape_fitness(oracle, ["1.2.2.2"]))


def test_changed_oracle_does_not_inherit_obsolete_cache():
    landscape = simulate_grf_landscape(n_chr=2, seed=42)
    original = predict_landscape_fitness(landscape, ["2.3"])
    modified = replace(landscape, founder_fitness=landscape.founder_fitness + 2)
    np.testing.assert_allclose(predict_landscape_fitness(modified, ["2.3"]), original + 2)
