# karstpy

karstpy reconstructs aneuploid karyotype fitness landscapes from cell counts
at two or more time points. It ports the sparse joint count model from
[KARST 1.1.3](https://github.com/TaoLee0510/karst/tree/6ae2a24). The four C++
likelihoods come from that source version, with a stable evaluation of the
shared forward integral for nearly equal growth rates.

## Install

Use Python 3.10 or later. Version 1.1.3 provides wheels for CPython 3.10–3.14
on Linux x86_64 and aarch64 and on macOS 15 or later on Intel and Apple
silicon. Once the release is available on PyPI, install it with:

```sh
python -m pip install 'karstpy==1.1.3'
python -m pip install 'karstpy[alfak]==1.1.3'  # add RDS export dependencies
```

Installation from a source checkout uses a C and C++17 compiler. pip
installs the build backend, CMake if needed, and pybind11. Runtime
dependencies are NumPy, SciPy, and pandas. R is used for reference checks
and can read the optional ALFA-K export.

```sh
python -m pip install .
python -m pip install '.[alfak]'  # add RDS export dependencies
```

## Fit and predict

Use a pandas DataFrame with karyotype labels as the row index and one column
per time point, in chronological order. Each label gives the copy numbers of
the same chromosomes in the same order. Counts must be non-negative integers.
You can also pass a matrix with `labels=`. Supply observation times with
`times=`, or the intervals between consecutive observations with `dt=`.

```python
import pandas as pd
from karstpy import check_sparse_joint_counts, fit_karst

counts = pd.DataFrame(
    [[180, 240], [36, 31], [2, 0], [6, 3]],
    index=["2.2", "2.3", "1.2", "3.2"],
)
print(check_sparse_joint_counts(counts))

fit = fit_karst(counts, dt=3, beta=5e-5, sparse_joint_control={"birth_rate": 1})
print(fit.predict(["2.2", "2.3", "4.4"])[["karyotype", "relative_fitness", "relative_sd"]])
```

`relative_fitness` is net growth relative to the count-weighted mean of the
major karyotypes, per unit of time supplied. `relative_sd` is the calibrated
standard deviation. `relative_sd_model` is the standard deviation from the
fitted model. The baseline birth rate and missegregation rate are supplied
parameters. Use `fit_sparse_joint_rates` to show sensitivity across a supplied
range of either rate.

With multiple time points, `fit_sparse_joint(..., series=...)` supports
`"anchored"`, `"pairwise"`, and `"state_space"`. The fitted domain, counts,
expectations and precision are in `fit.summary`. Use `fit.diagnostics` for
domain coverage, convergence, and the history profile when several candidate
histories were supplied. The [simulation notebook](examples/karst-simulated-landscape.ipynb)
walks through a seeded example.

## ALFA-K output

```python
from karstpy import as_alfak, write_alfak_outputs

output = as_alfak(fit, minobs=10, nboot=45, seed=1, landscape_data_output=True)
write_alfak_outputs(output, "alfak-results")
```

The RDS files contain the landscape table, posterior draws, bootstrap layout,
cross-validation results and optional landscape predictor. The Python fit is
saved as `karst_fit.pkl`. Use `alfak_compat` to fit directly from ALFA-K's
input layout and write its files in one call. RDS export uses rdata on Python
3.11 and later and an XDR writer for these output types on Python 3.10.

## Numerical validation

karstpy uses TMBad and Eigen for automatic differentiation and
sparse Cholesky, and SciPy L-BFGS-B for the bounded scale optimisation. Repeated
optimiser evaluations and individual posterior draws can differ from the R
implementation because the optimizer and Cholesky permutation differ. See
[the local parity workflow](tools/parity/README.md) for the measured
differences and the commands that reproduce the comparison.

```sh
python -m pip install '.[test,alfak]'
python -m pytest
```

## License

GNU General Public License version 2 only. See [LICENSE](LICENSE).
Upstream KARST source retains its [original MIT notice](LICENSE-karst).
