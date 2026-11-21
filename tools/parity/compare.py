"""Compare the locally built Python package with generated R fixtures."""

import argparse
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from karstpy import KarstConvergenceWarning, fit_sparse_joint, simulate_grf_landscape, simulate_sparse_counts
from karstpy.tmb import LaplaceObjective


def relative_error(value, expected):
    expected = np.asarray(expected)
    return np.max(np.abs(value - expected)) / max(1, np.max(np.abs(expected)))


def compare_native(path):
    data = json.loads(path.read_text())
    parameters = {k: np.asarray(v) for k, v in data["parameters"].items()}
    random = {name: np.ones_like(parameters[name], dtype=bool) for name in data["random"]}
    if data["model"] == "sparse_joint_pairwise":
        random["zeta"] = np.asarray(data["member"], dtype=bool)
    objective = LaplaceObjective(data["model"], data["data"], parameters, random, {"log_scale": data["free"]})
    records = []
    for i, case in enumerate(data["joint"], start=1):
        x = np.asarray(case["point"])
        hessian = objective.hessian(x, random=False).toarray()
        ref = np.asarray(case["he"])
        records.append(dict(
            model=data["model"], point=i, value_error=relative_error(objective.value(x), case["fn"]),
            gradient_error=relative_error(objective.gradient(x), case["gr"]), hessian_error=relative_error(hessian, ref),
            numeric_sparsity_identical=np.array_equal(abs(hessian) > 1e-10, abs(ref) > 1e-10),
        ))
    for record, case in zip(records, data["laplace"]):
        point = np.asarray(case["fixed"])
        value = objective.fn(point)
        value_mode = objective.last_par
        gradient = objective.gr(point)
        ref = np.asarray(case["best"])
        # TMB removes inactive parameters. The standalone tape retains constants.
        kept = np.sort(np.r_[objective.random, objective.fixed])
        record.update(laplace_error=relative_error(value, case["value"]),
                      marginal_gradient_error=relative_error(gradient, case["grad"]),
                      value_mode_error=relative_error(value_mode[kept], case["mode"]),
                      gradient_mode_error=relative_error(objective.last_par[kept], case["gradient_mode"]),
                      best_mode_error=relative_error(objective.last_par_best[kept], ref))
    return records


def compare_fit(root, case):
    name = case["name"]
    counts = pd.read_csv(root / f"{name}.counts.csv", index_col=0)
    counts.index = counts.index.astype(str)
    expected = pd.read_csv(root / f"{name}.summary.csv", dtype={"karyotype": str}, index_col=0).reset_index(drop=True)
    record = {key: case[key] for key in ("name", "case", "seed", "rate_group")}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", KarstConvergenceWarning)
            fit = fit_sparse_joint(counts, times=case["times"], beta=case["beta"], birth_rate=case["birth_rate"],
                                   max_cn=4, series=case["series"], history=case["history"])
    except (ValueError, RuntimeError) as error:
        return dict(record, error=str(error))
    observed = fit.summary
    domain = observed.karyotype.tolist() == expected.karyotype.tolist() and np.array_equal(observed.major, expected.major)
    merged = observed.merge(expected, on="karyotype", suffixes=("_python", "_r"))
    left = fit.diagnostics.get("dropped_samples", pd.DataFrame(columns=["sample"]))
    return dict(record, domain_equal=domain,
                history_equal=fit.diagnostics["history"] == case["selected_history"],
                coverage_equal=np.array_equal(left["sample"], np.atleast_1d(case["left_out"])) or len(left) == 0 and not case["left_out"],
                convergence_equal=fit.diagnostics["convergence"] == case["convergence"],
                r_convergence=case["convergence"], python_convergence=fit.diagnostics["convergence"],
                r_history=case["selected_history"], python_history=fit.diagnostics["history"],
                objective_error=relative_error(fit.diagnostics["objective"], case["objective"]),
                posterior_objective_error=relative_error(fit.diagnostics["posterior_objective"], case["posterior_objective"]),
                fitness_error=(merged.relative_fitness_python - merged.relative_fitness_r).abs().max(),
                relative_sd_error=(merged.relative_sd_python - merged.relative_sd_r).abs().max())


def compare_simulations(root):
    results = []
    for seed in range(1, 11):
        landscape = simulate_grf_landscape(n_chr=3, min_cn=1, max_cn=4, seed=seed)
        simulated = simulate_sparse_counts(landscape, beta=0.01, dt=2, seed=seed, n0=90, n1=110,
                                          initial_population=250, overdispersion=0 if seed % 2 else 70)
        expected = pd.read_csv(root / f"sim-{seed}.csv", index_col=0)
        expected.index = expected.index.astype(str)
        truth = np.fromfile(root / f"sim-{seed}.bin", dtype="<f8")
        results.append(dict(seed=seed, exact_counts=simulated["counts"].equals(expected),
                            exact_truth=np.array_equal(truth, simulated["truth_observed"])))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixtures", type=Path)
    parser.add_argument("--native-only", action="store_true", help="check the four native templates only")
    parser.add_argument("--require-native-parity", action="store_true", help="fail when a native comparison exceeds tolerance")
    args = parser.parse_args()
    root = args.fixtures
    native = pd.DataFrame(record for path in sorted(root.glob("native-*.json")) for record in compare_native(path))
    if len(native) != 12 or native.model.nunique() != 4:
        raise SystemExit("Expected 12 native points across four likelihood templates.")
    native.to_csv(root / "native-comparison.csv", index=False)
    maxima = native.select_dtypes(include="number").drop(columns="point").max()
    print("Native maxima")
    print(maxima.to_string())
    if args.require_native_parity:
        joint = ["value_error", "gradient_error", "hessian_error"]
        marginal = ["laplace_error", "marginal_gradient_error", "value_mode_error", "gradient_mode_error", "best_mode_error"]
        if (not native.numeric_sparsity_identical.all()
                or not (native[joint] <= 1e-10).all().all()
                or not (native[marginal] <= 1e-8).all().all()):
            raise SystemExit("Native parity exceeded the tolerance. See native-comparison.csv.")
        print("Native tolerance checks passed.")
    if args.native_only:
        return
    simulations = pd.DataFrame(compare_simulations(root))
    simulations.to_csv(root / "simulation-comparison.csv", index=False)
    cases = json.loads((root / "cases.json").read_text())
    fits = pd.DataFrame(compare_fit(root, case) for case in cases)
    fits.to_csv(root / "fit-comparison.csv", index=False)
    print("Simulation matches")
    print(simulations[["exact_counts", "exact_truth"]].sum().to_string())
    columns = ["domain_equal", "history_equal", "coverage_equal", "convergence_equal"]
    print("Fit matches")
    print(fits.groupby("case")[columns].sum().to_string())
    print("Fit discrepancies")
    errors = ["objective_error", "posterior_objective_error", "fitness_error", "relative_sd_error"]
    print(fits.groupby("case")[errors].max().to_string())
    if "error" in fits:
        print(fits.loc[fits.error.notna(), ["name", "error"]].to_string(index=False))


if __name__ == "__main__":
    main()
