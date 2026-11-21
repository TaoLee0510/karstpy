# R math routines

The source files come from R 4.4.2, principally from `src/nmath`.
They use the GNU General Public License, version 2 or later. The original
notices appear in the source files. The parent project includes the GPL-2
text in `LICENSE`.

`Rmath.h` is configured from R's `Rmath.h0.in`. The local `Rconfig.h`
selects the subset of features needed on Linux and macOS.
`cpp/rng.cpp` supplies the standalone uniform generator using the default
Mersenne Twister and rejection sampling from R's `src/main/RNG.c`.
