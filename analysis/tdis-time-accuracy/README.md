# TDIS time accuracy analysis

Notes and scripts for reproducing the analysis of simulation time accuracy in TDIS
(`src/Timing/tdis.f90`, `src/Timing/ats.f90`). Related: #3204 (first step length), #3202
(ATS retry at end of simulation).

| File | Contents |
|---|---|
| `analysis.md` | The analysis: error sources, current behavior, options, measurement plan |
| `first_step_errors.py` | First-step error tables and the 20,000-case random survey |
| `step_errors.py` | Middle-step and last-step error within one geometric period |
| `plot_errors.py`, `figures/` | Figures comparing the first-step methods (current, #3204, `expm1`/`log1p`, hybrid) |
| `first_step.f90`, `first_step_cmp.py` | Same first-step forms in Fortran, compared with Python |
| `tdis_val_repro.py` | Writes small models that show the zero-`PERLEN` validation bug |

The scripts use only the Python standard library, except `tdis_val_repro.py`, which needs
FloPy, and `plot_errors.py`, which needs numpy and matplotlib. Run from the repository root.

## First-step and per-step errors

```shell
python3 -I analysis/tdis-time-accuracy/first_step_errors.py
python3 -I analysis/tdis-time-accuracy/step_errors.py
```

The figures need numpy and matplotlib (both in the pixi environment) and take a few seconds:

```shell
pixi run python -I analysis/tdis-time-accuracy/plot_errors.py
```

## Fortran comparison

Compiles the first-step forms (integer exponent, real exponent, `expm1`/`log1p` through
`bind(c)`) and checks they are bit-identical to the Python evaluation. Inside the pixi
environment, so the compiler and linker are on the path:

```shell
pixi run bash -c '
  gfortran -O1 -o /tmp/first_step analysis/tdis-time-accuracy/first_step.f90 &&
  /tmp/first_step | python3 -I analysis/tdis-time-accuracy/first_step_cmp.py'
```

Try other optimization levels and compilers. To keep the compiler from folding the literal
inputs, pass `vol` and feed the inputs on stdin:

```shell
printf '1.0 10 1.0001\n1.0 10 1.000001\n1.0 3 1.000000001\n100.0 50 0.999\n10.0 100 1.005\n' |
  /tmp/first_step vol | python3 -I analysis/tdis-time-accuracy/first_step_cmp.py
```

Checked so far: gfortran 15.3 on arm64 macOS at -O0, -O1, and -O3, with both kinds of input,
all bit-identical. ifx has not been tested.

## Zero-PERLEN validation bug

```shell
pixi run python -I analysis/tdis-time-accuracy/tdis_val_repro.py <workspace> bin/mf6
```

This writes four simulations into `<workspace>`. Run `mf6` in each, with a time limit in case
a model hangs:

```shell
for d in <workspace>/*/; do (cd "$d" && perl -e 'alarm 20; exec @ARGV' "$OLDPWD/bin/mf6"); done
```

Observed with a `develop` build from November 2026: `neg_perlen_control` fails with a PERLEN
error, `neg_perlen_after_zero` and `neg_tsmult_after_zero` end in normal termination with
invalid times, and `zero_nstp_after_zero` crashes with SIGILL.
