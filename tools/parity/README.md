# Local R comparison

The full fit and simulation reference is KARST 1.1.3 at commit `6ae2a24`,
built with TMB 1.9.19 and R 4.4.2. Use an isolated library with the reference
package installed. The generator also needs R's `jsonlite`. The fixture
folders below are ignored by Git and excluded from the source distribution.

```sh
Rscript --vanilla tools/parity/generate.R _reference/Rlib _reference/build/parity
python tools/parity/compare.py _reference/build/parity
python tools/parity/measure_domain.py
```

The fit comparison uses ten seeds for each of the two-time-point, short
interval, anchored, pairwise, and state-space cases. It also compares all
three rate points for ten externally specified birth-rate ranges. The
comparison writes CSVs for the native checks, simulations, and fitted
results. Objective, derivative, and mode errors use the maximum absolute
distance divided by the larger of one and the maximum absolute reference
value. Fitness and standard deviation differences are absolute.

## Fixed-parameter Laplace checks

The isolated derivative reference applies the stabilization introduced in
Python dev2 and keeps R's CppAD backend. The checked-in patch reproduces that
reference independently of the Python operator. Start with the phase-0 R
build source in `_reference/build/karst-1.1.3` and use a separate build tree
and library.

```sh
mkdir -p _reference/build/karst-laplace-stable _reference/Rlib/laplace-stable
rsync -a --exclude='*.o' --exclude='*.so' --exclude='.DS_Store' \
  _reference/build/karst-1.1.3/ _reference/build/karst-laplace-stable/
patch -d _reference/build/karst-laplace-stable -p1 < tools/parity/laplace-stable.patch
R CMD INSTALL --preclean --library=_reference/Rlib/laplace-stable _reference/build/karst-laplace-stable
Rscript --vanilla tools/parity/generate.R _reference/Rlib/laplace-stable \
  _reference/build/parity-stable --native-only
python tools/parity/compare.py _reference/build/parity-stable --native-only --require-native-parity
```

The strict check covers 12 points across the four likelihoods. Joint
objective, gradient, and Hessian errors must be at most `1e-10`, with the
same numerical Hessian nonzero pattern. The Laplace value, gradient, and
latent modes must agree to `1e-8`. The native fixture JSON retains 17 digits
for double-precision round trips.

The original R integral expression loses precision when the growth rates
nearly coincide. A `2e-14` change in one latent log birth rate shifted the
original R marginal gradient by up to `7.26e-4` for anchored and `2.69e-3`
for state-space at fixed outer parameters. The stabilized evaluation uses a
series in the squared rate difference near zero, which avoids cancellation
in the third derivatives needed by the Laplace gradient.

## Results on macOS arm64

| Check | Result |
|---|---|
| Domains, major karyotypes, histories, left-out time points, and convergence labels | All 80 fits agreed |
| Maximum absolute relative fitness difference | `2.67e-6` |
| Maximum scaled objective difference at the reported posterior mode | `4.24e-8` |
| Same-seed simulations | Counts and truth matched bit for bit in all 10 cases |
| Joint objective, gradient, and Hessian against stabilized R | Maximum scaled errors `5.16e-14`, `2.59e-15`, and `4.50e-16`. All numerical nonzero patterns agreed |
| Marginal Laplace value and gradient against stabilized R | Maximum scaled errors `5.74e-14` and `4.35e-11`. The gradient errors for anchored and state-space were at most `1.56e-11` and `1.55e-11`, respectively |
| Marginal Laplace value and gradient against original R | Maximum scaled errors `5.40e-9` and `1.39e-3` |

The bounded scale optimizer is L-BFGS-B. R's PORT source credits point to
Netlib, where access to the source is described but GPL-2-compatible
redistribution rights for these routines are unclear. The stage-2 iteration
sequence condition applies to a PORT implementation. See the
[Netlib PORT distribution notes](https://www.netlib.org/port/readme).

Eigen's Cholesky permutation can differ from R Matrix's CHOLMOD permutation.
This changes individual draws at the same seed. The export checks verify
the reported means and uncertainties, normalized frequencies, named matrix
dimensions, retry metadata, and the R prediction method on the saved
predictor.

```sh
Rscript --vanilla tools/parity/check_alfak.R _reference/Rlib path/to/alfak-results
```

## 5,000-state benchmark

The benchmark evaluates one fixed scale on a 22-chromosome domain with
10,047 latent variables and 9,700 directed edges. These are medians of three
interleaved runs in separate processes on the same macOS arm64 host with
CPython 3.13, using the dev1 and dev2 wheels and the dev3 Release build.

| Version | Setup | Laplace value and gradient | Peak RSS | Hessian tape operators |
|---|---:|---:|---:|---:|
| dev1 | 4.674 s | 4.430 s | 4,099 MiB | 25,732,319 |
| dev2 | 6.387 s | 5.832 s | 5,072 MiB | 35,255,444 |
| dev3 | 3.866 s | 4.100 s | 3,555 MiB | 21,039,656 |

Dev2 recorded both the polynomial and exponential branches of each forward
integral and used conditional expressions to select the evaluation. Taking
the gradient and sparse Hessian duplicated their operations across the
derivative tapes. The Hessian operator counts were measured separately
after tape optimization. All three builds had 1,201,598 structural Hessian
entries, and the benchmark's peak RSS was already reached during setup in
each run.

Dev3 records the smooth log-integral correction as a unary TMBad operator.
Each sweep evaluates the relevant branch using a series near equal growth
rates or the closed form farther away, with analytic derivatives through
third order for the Laplace gradient. The benchmark reports setup and the
two Laplace calls separately.
