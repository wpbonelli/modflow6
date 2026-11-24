# TDIS time accuracy analysis

Notes and scripts for reproducing the analysis of simulation time accuracy in TDIS
(`src/Timing/tdis.f90`, `src/Timing/ats.f90`). Related: #3204 (first step length), #3202
(ATS retry at end of simulation).

| File | Contents |
|---|---|
| `analysis.md` | The analysis: error sources, current behavior, options, measured drift |
| `first_step_errors.py` | First-step error tables and the 20,000-case random survey |
| `step_errors.py` | Middle-step and last-step error within one geometric period |
| `plot_errors.py`, `figures/` | Figures comparing the first-step methods (current, #3204, `expm1`/`log1p`, hybrid) |
| `first_step.f90`, `first_step_cmp.py` | Same first-step forms in Fortran, compared with Python |
| `tdis_drift_instrumentation.patch` | Throwaway TDIS logging for the drift measurement (apply to `develop`) |
| `build_cases.py`, `run_drift.sh` | Synthetic drift models, and a runner for each first-step method |
| `analyze_drift.py`, `plot_drift.py`, `drift_summary.csv` | Drift against exact references, figures, and the measured results |
| `tdis_val_repro.py` | Writes small models that show the zero-`PERLEN` validation bug |

The scripts use only the Python standard library, except `tdis_val_repro.py` and
`build_cases.py`, which need FloPy, and `plot_errors.py` and `plot_drift.py`, which need numpy
and matplotlib. Run from the repository root.

## First-step and per-step errors

```shell
python3 -I analysis/tdis-time-accuracy/first_step_errors.py
python3 -I analysis/tdis-time-accuracy/step_errors.py
```

The figures need numpy and matplotlib (both in the pixi environment) and take a few seconds:

```shell
pixi run python -I analysis/tdis-time-accuracy/plot_errors.py
```

## Drift measurement

Build an instrumented `mf6` in a separate worktree. The ATS cases also need the #3202 fix
(`fix-tdis-endofsim-retry`); without it, a retry on the final step ends the simulation early.

```shell
git worktree add -b tdis-drift ../mf6-drift develop
cd ../mf6-drift
git am ../modflow6/analysis/tdis-time-accuracy/tdis_drift_instrumentation.patch
pixi run --manifest-path ../modflow6/pixi.toml meson setup _bd
pixi run --manifest-path ../modflow6/pixi.toml meson compile -C _bd ./src/mf6:executable
```

Build the synthetic models, copy in the real ones, and run each under all three first-step
methods. Set `MF6_TDIS_FORCE_RETRY` only for the ATS cases. The variables are described in
`analysis.md`.

```shell
D=<workspace>; T=analysis/tdis-time-accuracy; MF6=$PWD/../mf6-drift/_bd/src/mf6
pixi run python -I $T/build_cases.py $D/cases
mkdir -p $D/real && for m in \
    modflow6-testmodels/mf6/{test020_NevilleTonkinTransient,test028_sfr,test045_lake2tr_dev} \
    modflow6-testmodels/mf6/{test051_uzfp3_lakmvr_v2_dev,test061_csub_holly,test205_gwtbuy-henrytidal} \
    modflow6-largetestmodels/{test1200_gwtbuy-goswami,test1201_gwtbuy-elderRa60} \
    modflow6-examples/examples/{ex-gwe-prt,ex-gwf-sagehen,ex-gwt-keating}; do
  cp -R ../$m $D/real/; done
R="$(ls -d $D/real/test* $D/real/ex-gwf-sagehen) $D/real/ex-gwe-prt/gwe $D/real/ex-gwt-keating/mf6gwt"
$T/run_drift.sh $MF6 $(ls -d $D/cases/* | grep -v ats_) $R
MF6_TDIS_FORCE_RETRY=10 $T/run_drift.sh $MF6 $D/cases/ats_*
python3 -I $T/analyze_drift.py $D/out $D/cases/* $R
pixi run python -I $T/plot_drift.py $D/out $T/figures
```

The GWE and GWT examples run their transport simulation from the flow output already in the
example directory. Analysis takes about 20 seconds for all 31 models.

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
