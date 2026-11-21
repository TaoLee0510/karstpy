"""Time a fixed-scale evaluation on a 5,000-state, 22-chromosome domain."""

from itertools import islice
from time import perf_counter
import resource
import sys

import numpy as np
import pandas as pd

from karstpy import prepare_karst_data
from karstpy._problem import build_problem, domain_from_labels, one_step_neighbours
from karstpy.karyotypes import format_karyotypes


def peak_mib():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        peak *= 1024
    return peak / 1024**2


def main():
    rng = np.random.default_rng(42)
    major = rng.integers(1, 4, size=(150, 22))
    labels = format_karyotypes(major)
    major_labels = set(labels)
    neighbors = (label for label in one_step_neighbours(major, 1, 4) if label not in major_labels)
    all_labels = labels + list(islice(neighbors, 5000 - len(labels)))
    assert len(all_labels) == 5000
    counts = pd.DataFrame(rng.poisson(50, size=(len(labels), 2)), index=labels)
    data = prepare_karst_data(counts, dt=2, beta=0.001)
    domain = domain_from_labels(data.counts, all_labels, labels)
    settings = dict(
        max_cn=4, min_cn=1, birth_rate=1, shell_power=0.85, birth_rate_log_sd=0.05, shell_sd=0.5,
        initial_frequency=(np.log(1e-5), 3), initial_state="balance", history=2, balance_sd=1, inner_control={},
    )
    start = perf_counter()
    problem = build_problem(data, settings, domain)
    prepared = perf_counter()
    objective = problem.objective(settings)
    taped, setup_peak = perf_counter(), peak_mib()
    value = objective.fn(objective.par)
    valued, value_peak = perf_counter(), peak_mib()
    gradient = objective.gr(objective.par)
    finished, gradient_peak = perf_counter(), peak_mib()
    assert np.isfinite(value) and np.isfinite(gradient).all()
    print(f"{len(all_labels)} states, {objective.random.size} latent variables, {problem.kernel.source.size} edges")
    print(f"Model data {prepared - start:.3f} s")
    print(f"Native tapes and sparse pattern {taped - prepared:.3f} s, peak RSS {setup_peak:.1f} MiB")
    print(f"Laplace objective {valued - taped:.3f} s, peak RSS {value_peak:.1f} MiB")
    print(f"Laplace gradient {finished - valued:.3f} s, peak RSS {gradient_peak:.1f} MiB")
    print(f"Tape and factor pattern {taped - start:.3f} s")
    print(f"Laplace objective and gradient {finished - taped:.3f} s")
    print(f"Peak RSS {gradient_peak:.1f} MiB")


if __name__ == "__main__":
    main()
