# TDIS validation stops at the first zero-length period

## Summary

`check_tdis_timing` (`src/Timing/tdis.f90:600`) validates each stress period's `PERLEN`,
`NSTP`, and `TSMULT`. When it reaches a period with `PERLEN = 0`, it writes a warning and
**returns** instead of moving on to the next period (line 639). Later periods are not checked
at all, so invalid input there is accepted. Zero-length periods are common, usually as an
initial steady-state period, so this affects many real models.

Found while auditing TDIS time accuracy (`analysis.md`). It is a validation
bug, not an accuracy problem.

## Code

```fortran
      ! -- Warn if perlen is zero
      if (perlen(kper) == DZERO) then
        write (iout, fmtpwarn) kper
        return
      end if
```

The other checks in the loop (`NSTP <= 0`, `TSMULT <= 0`, `PERLEN < 0`, step too small) store an
error and also `return`. The routine's docstring says this is so the caller can attach the file
name before stopping. For errors that is acceptable (only the first error gets reported). For
the zero-length warning it is wrong: no error has been stored, so the caller carries on and the
simulation runs.

The `return` was probably meant to skip the step-length loop for that period, since a zero
`dt` would fail the `tstart == tend` check. `cycle` does that.

The checks are also ordered so that `PERLEN = 0` is tested before `TSMULT`. A zero-length period
with `TSMULT <= 0` is never checked, and with `TSMULT = -1` and even `NSTP` the first step
formula is `0 * 2 / (1 - 1) = 0/0` (NaN).

## Reproduction

Two-period GWF models run with `bin/mf6` (built from `develop` plus the current branch, gfortran,
arm64). Period 1 is steady state, period 2 transient.

| Period data | Result |
|---|---|
| `(1.0, 1, 1.0), (-10.0, 1, 1.0)` (control) | Error: "PERLEN cannot be less than 0.0 for stress period 2" |
| `(0.0, 1, 1.0), (-10.0, 1, 1.0)` | Runs, step length -10, total time -10, "Normal termination" |
| `(0.0, 1, 1.0), (10.0, 3, -2.0)` | Runs, steps 3.33, -6.67, 13.33 (time goes backward mid-period), "Normal termination" |
| `(0.0, 1, 1.0), (10.0, 0, 1.0)` | SIGILL (floating-point trap from `perlen / float(0)`), backtrace, no error message |

The repro script is `tdis_val_repro.py` on the `tdis-time-accuracy-analysis` branch, in
`analysis/tdis-time-accuracy/`. The tests below replace it.

## Fix

In `check_tdis_timing`:

1. Replace the `return` after the zero-length warning with `cycle`, so later periods are still
   checked. `tstart` stays where it is, which is right for a zero-length period.
2. Move the zero-length check after the `NSTP`, `TSMULT`, and negative `PERLEN` checks, so a
   zero-length period still has those validated.
3. Optional: replace the error `return`s with `cycle` too, so every bad period is reported in one
   run instead of one per run. `store_error` already accumulates errors, and
   `tdis_source_timing` checks `count_errors()` after the call, so the docstring's reason for
   returning doesn't require stopping at the first one. Skip the step-length loop for any period
   that failed a check, to avoid follow-on errors. Update the docstring.

No change in behavior for valid input. Models that currently run with invalid input after a
zero-length period will now fail with an error. That is the intended effect, but it needs a
release note.

## Tests

Add to `autotest/test_gwf_errors.py`, which already builds a simple model with configurable
`tdiskwargs` and checks for expected error messages:

- `(0.0, 1, 1.0), (-10.0, 1, 1.0)`: expect "PERLEN cannot be less than 0.0 for stress period 2".
- `(0.0, 1, 1.0), (10.0, 1, -2.0)`: expect the `TSMULT` error for period 2.
- `(0.0, 1, 1.0), (10.0, 0, 1.0)`: expect the `NSTP` error for period 2 (currently a crash).
- `(0.0, 2, -1.0)`: expect the `TSMULT` error for period 1 (ordering fix).
- If item 3 is done: two bad periods, expect both errors.
- A valid model with a zero-length first period still runs (covered by existing tests, e.g.
  `test_gwf_pertim.py`; confirm one is in the targeted run).

Run `test_gwf_errors.py` and a few existing tests with zero-length periods (from the grep:
`test_gwf_pertim.py`, `test_gwf_maw04.py`, `test_gwf_newton01.py`).

## Related, not in scope

- The step-length check compares each step with the time at the start of the period, not the
  time the step actually starts at: `tstart` is only advanced once per period (line 677). A step
  that is negligible next to its real start time later in the period isn't caught. Making it
  stricter could reject models that run today, so it should be a separate change.
- The step-length check duplicates the step formula from `tdis_set_delt` without the last-step
  correction. See #3204 and `analysis.md`.
- The zero-length warning says "PERLEN must not be zero for transient periods", but TDIS doesn't
  know which periods are transient, so a transient zero-length period still only gets a warning.
