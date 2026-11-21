# karstpy 1.1.3

- First stable Python release of the KARST 1.1.3 sparse joint models,
  simulation, and ALFA-K export.
- Provides wheels for CPython 3.10–3.14 on Linux x86_64 and aarch64 and
  macOS 15 or later on Intel and Apple silicon.
- Uses the original R package for fit and simulation comparisons and an
  isolated R build with the stable integral for fixed-point Laplace checks.

# karstpy 1.1.3.dev3

- Record the stabilized forward integral as a compact TMBad operator with
  derivatives through third order. This reduces sparse Hessian tape size,
  setup time, and memory on large domains.
- Report separate setup, Laplace value, and Laplace gradient timings in the
  5,000-state benchmark.

# karstpy 1.1.3.dev2

- Stabilized the shared division integral used by the Laplace approximation.
  The new evaluation avoids cancellation in the higher derivatives when
  adjacent karyotypes have nearly equal growth rates.
- Exported reference parameters at round-trip precision and added regression
  tests for the anchored and state-space Laplace gradients.

# karstpy 1.1.3.dev1

- Implemented sparse joint model fitting, prediction, and rate sensitivity.
- Added two-time-point, anchored, pairwise, and state-space likelihoods using
  the original KARST C++ templates.
- Added seeded landscape and count simulation using R 4.4.2 random streams.
- Added ALFA-K-compatible summaries, posterior draws, and RDS files.
- Added a runnable notebook, Python tests, local R comparisons, and a Linux
  and macOS wheel build workflow.

# karstpy 1.1.3.dev0

- Established the Python source and wheel build and the isolated R 1.1.3
  reference environment.
