import pickle
from pathlib import Path

import numpy as np
import pytest

from karstpy import alfak_compat, as_alfak, write_alfak_outputs
from karstpy.alfak import alfak_neighbours


def test_alfak_neighbor_order():
    assert alfak_neighbours(["2.2", "2.3"]) == ["1.2", "3.2", "2.1", "1.3", "3.3", "2.4"]


def test_alfak_outputs_preserve_means_medians_sds_and_frequencies(fit):
    output = as_alfak(fit, minobs=4, nboot=9, n0=None, nb=None, landscape_data_output=True)
    calls = output.landscape.k.to_list()
    predicted = fit.predict(calls)
    np.testing.assert_allclose(output.landscape["mean"], predicted.fitness_mean)
    np.testing.assert_allclose(output.landscape["median"], predicted.fitness_mean)
    np.testing.assert_allclose(output.landscape.sd, predicted.relative_sd)
    np.testing.assert_allclose(output.bootstrap_res["final_frequencies"].sum(axis=1), 1)
    np.testing.assert_allclose(output.landscape_data["mean"].predict(calls).ravel(), output.landscape["mean"])
    assert output.rebuild["xval_folds"]["rebuilt"] > 0


def test_correct_efflux_preserves_requested_growth_level(fit):
    output = as_alfak(fit, minobs=4, nboot=7, n0=100, nb=400, correct_efflux=True, first_time=2, xval=False)
    central_fitness = output.bootstrap_res["final_fitness"][-1]
    frequency = output.bootstrap_res["final_frequencies"][-1]
    assert frequency @ central_fitness == pytest.approx(np.log(4) / 3)


def test_alfa_k_files_can_be_read_and_fit_can_be_reloaded(fit, tmp_path):
    rdata = pytest.importorskip("rdata", minversion="1.0")
    output = as_alfak(fit, minobs=4, nboot=7, landscape_data_output=True)
    paths = write_alfak_outputs(output, tmp_path)
    assert len(paths) == 6
    loaded = rdata.read_rds(tmp_path / "landscape.Rds")
    np.testing.assert_allclose(loaded["mean"], output.landscape["mean"])
    with (tmp_path / "karst_fit.pkl").open("rb") as file:
        restored = pickle.load(file)
    np.testing.assert_allclose(restored.summary.relative_fitness, fit.summary.relative_fitness)


def test_alfak_wrapper_validates_counts_columns(counts, tmp_path):
    with pytest.raises(ValueError, match="numeric column"):
        alfak_compat({"x": counts.rename(columns=lambda x: f"t{x}"), "dt": 3}, Path(tmp_path), nboot=5)
