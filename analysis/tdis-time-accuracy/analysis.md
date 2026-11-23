# Simulation time accuracy in TDIS

## Summary

MF6 advances simulated time by repeated floating-point addition (`totim = totim + delt`). The
rounding error in `totim` therefore grows with the number of additions, not only with the size
of `t`. Some of this error is inherent to floating point, but the growth with step count is not.
TDIS already partly corrects for it at stress period boundaries, though not for ATS periods and
not across periods.

Computing times directly from the inputs (closed form where possible, compensated summation
elsewhere) would bound the error at a few ulps of `t` for every step. The change is small in
code but affects every model, because `totim` feeds output, time series, observations, and
exchanges. This note describes the current behavior and the options, and proposes measuring the
actual drift before deciding.

`totim` and `delt` can't both be accurate and exactly consistent (`totim(k) = totim(k-1) +
delt(k)`) in floating point. The current scheme mostly favors accurate `delt` and lets `totim`
drift; the options below need to choose deliberately (see [Decoupling `totim` and
`delt`](#decoupling-totim-and-delt)).

This came up while working on PRT release time deduplication. PRP compares `totimc` with
release times parsed from input and needs a tolerance to absorb the difference. A fixed absolute
default can't suit every time scale, and a more accurate `totim` would let that tolerance shrink
to a small relative value.

Related issue: [#3204](https://github.com/MODFLOW-ORG/modflow6/issues/3204) (real exponent for
the TDIS first time step length). It fixes one of the error sources below and can land
independently (see [Relationship to #3204](#relationship-to-3204)).

## What is unavoidable

Doubles near `t` are spaced about `eps * |t|` apart (eps ≈ 2.2e-16). Any time stored as a
double has relative error of about eps no matter how it was computed, so absolute error
necessarily grows with `t`. The only alternative is exact representation, such as integer
ticks or rationals. That doesn't fit MF6: `perlen`, `tsmult`, and ATS produce arbitrary real
step lengths, so there is no common tick unit.

A consequence: any step length computed as a difference of two times, `t1 - t0`, has absolute
error of about `ulp(t1)`, so its relative error is about `ulp(t1) / (t1 - t0)`. Small steps late
in a long simulation lose many digits this way. This matters for the last step of each period
(below) and for option 2.

## What is avoidable

- Error that grows with the **number of operations**. Each `totim + delt` can add up to half an
  ulp of `totim`, so after N additions the worst case is about `N * eps * t / 2` (typically
  closer to `sqrt(N) * eps * t`). Geometric stepping adds more: `delt = tsmult * delt`
  compounds a relative error in `delt` itself at every step.
- Error in **how the first step is computed** (#3204 and the `tsmult` near 1 case below).

## Current behavior (`src/Timing/tdis.f90`)

- **Per step.** `tdis_set_timestep` sets `totimc = totim`, then `totim = totimsav + delt` and
  `pertim = pertimsav + delt` (lines 181 to 186).
- **Fixed stepping (`tdis_set_delt`).**
  - The first step uses a closed-form `delt` (lines 252 to 259), but the closed form is not
    accurate. It has two error sources: the integer exponent in `tsmult**nstp` (#3204), and
    cancellation in `1 - tsmult**nstp` when `tsmult` is near 1, which #3204 doesn't fix. See
    [First time step length](#first-time-step-length) for the mechanisms, measurements, and
    the fix.
    - The constant-step branch uses `float(nstp)` (line 258), a single-precision conversion.
      It is exact for any realistic `NSTP`; #3204 changes it to `real(nstp, DP)` for clarity.
  - Middle steps compound with `delt = tsmult * delt` (line 264). With an accurate first step,
    the measured middle-step error stays within about 5 to 7 ulp of `delt` over 100 to 200
    steps, and `totim` within 2 to 4 ulp.
  - The last step sets `delt = topertim + perlen - totim` "to avoid accumulation errors"
    (line 262). This brings `totim` back to `topertim + perlen` within one rounding at the end
    of each period. **The cost is in `delt`:** the last step absorbs all of the period's
    accumulated `totim` error, and it is computed as a difference of times, so its relative
    error is about `ulp(perlen) / delt_last`. For `tsmult < 1` the last step is the smallest, so
    the error is amplified:

    | PERLEN | NSTP | TSMULT | last `delt` error (current) | with accurate first step |
    |---|---|---|---|---|
    | 10 | 100 | 1.005 | 4205 ulp (9e-13 rel) | 19 ulp |
    | 10 | 100 | 0.95 | 9653 ulp (1.3e-12 rel) | 9653 ulp |
    | 1000 | 200 | 0.97 | 3505 ulp (7e-13 rel) | 3505 ulp |

    With `tsmult > 1`, most of the last-step error comes from the first step (#3204). With
    `tsmult < 1`, it is dominated by cancellation in the subtraction and #3204 doesn't change
    it. `delt` feeds storage terms, so this is a (small) error in the physics of that step, not
    only in reported time.
  - `topertim = totim` at the start of the period (line 251) is itself carried over from the
    previous period. So rounding still accumulates across periods, at about one rounding per
    period.
- **`check_tdis_timing` (lines 600 to 682)** repeats the step calculation to validate step
  lengths: integer exponent (line 662), `float(nstp)` (line 659) and compounding, but no
  last-step correction. Any change to how steps are computed must update both copies. A shared
  helper would keep them in sync. The check is also weaker than it looks: `tend = tstart + dt`
  (line 666) compares each step with the time at the start of the *period*, because `tstart`
  is only advanced once per period (line 677). A step that is negligible next to the time it
  actually starts at, later in the period, is not caught.
- **ATS periods.**
  - `ats_set_delt` clips the last step so that `pertim` lands on `perlen`
    (`ats.f90:632`), and `pertim` is snapped to `perlen` at the end of the period
    (`tdis.f90:198`).
  - `totim` gets no matching correction, so its error accumulates over all ATS steps in all
    periods.
  - `topertim` is set to zero at `kstp == 1` (line 168) and is only updated in `tdis_set_delt`,
    which isn't called for adaptive periods. So `topertim` isn't the period start time during
    ATS periods. It is only used inside TDIS, so this looks harmless, but it should be
    confirmed.
  - **Retry path (`tdis_delt_reset`) differs from the normal path.** After a failed step, ATS
    divides `delt` by `dtfailadj` (`ats.f90:669`) and `tdis_delt_reset` recomputes `totim` and
    `pertim` from the saved values. Then:
    - `ats_set_endofperiod` ends the period when `abs(pertim - perlen) < dtmin`
      (`ats.f90:699`). A retried step can land short of `perlen` but within `dtmin` of it: for
      example, a remaining `1.5*dtmin` retried with `dtfailadj = 1.2` gives `delt = 1.25*dtmin`
      and stops `0.25*dtmin` short. Unlike `tdis_set_timestep`, `tdis_delt_reset` doesn't snap
      `pertim` to `perlen`, and neither path snaps `totim`. So the period ends with `totim` up to
      `dtmin` behind, and the next period starts there. This error is **proportional to
      `dtmin`, not to roundoff**, and it carries forward to all later periods.
    - In the last period, `totim = totalsimtime` (line 239) hides this at the end of the
      simulation. That may be why the snap exists only on the retry path (see open questions).
      The snap doesn't adjust `delt`, so `totim /= totimsav + delt` on that step.
    - `endofsimulation` is set to `.true.` but never reset to `.false.` anywhere (only at
      allocation, line 499). If the first attempt at the final step reaches the end of the
      simulation and fails, and the retry with a shorter `delt` no longer reaches it,
      `endofperiod` is reset but `endofsimulation` stays `.true.`. The time loop in
      `mf6core.f90:47` would then stop after that step, before `totalsimtime`. Confirmed and
      fixed in [#3202](https://github.com/MODFLOW-ORG/modflow6/pull/3202) (open):
      `test_gwf_ats03.py` stopped at t = 0.8 of a 100-day period. #3202 also snaps `pertim` to
      `perlen` in `tdis_delt_reset`. It doesn't snap `totim` for periods before the last, so the
      `dtmin`-sized `totim` shortfall above remains, and `pertim` and `totim - period start` now
      disagree by up to `dtmin` on that step instead of agreeing.
- **End of simulation.**
  - `totim = totalsimtime` is applied only in `tdis_delt_reset` (line 239), the path for
    retrying a failed step.
  - The normal path, `tdis_set_timestep`, sets `endofsimulation` but doesn't snap `totim`.
    I haven't checked whether that's intentional.
- **`totalsimtime`.** It is a plain running sum of `perlen` (line 588), so its error grows with
  the number of periods. It is also summed differently from how `totim` reaches the end
  (`topertim + perlen - totim` added to `totim`), so the two can differ by an ulp even without
  ATS.

Net effect: in non-ATS runs the `totim` error is bounded by steps per period within a period,
plus the number of periods across the run, and the last step of each period carries a `delt`
error of about `ulp(perlen) / delt_last`. In ATS runs the roundoff error is bounded by the total
number of ATS steps, plus up to `dtmin` per period on the retry path. These are upper bounds.
Apart from the tables in this note, which cover single periods computed outside MF6, the actual
numbers haven't been measured.

## First time step length

For a geometric period, `tdis_set_delt` (line 254) and `check_tdis_timing` (line 660) compute

```fortran
delt = perlen(kper) * (DONE - tsmult(kper)) / (DONE - tsmult(kper)**nstp(kper))
```

Write `r = tsmult`, `n = nstp`, and `h = r - 1`. For `r` in [0.5, 2], which covers any
realistic `TSMULT`, `h = r - 1` (and `1 - r`) is computed exactly, so the numerator adds no
error beyond the final multiply. All the error is in the denominator, from two separate sources.

### Source 1: the integer exponent (#3204)

gfortran evaluates `r**n` with integer `n` by repeated squaring (`__powidf2` in libgcc). That
takes about `log2(n)` squarings plus some multiplies, each rounding with up to eps/2 relative
error. Squaring also doubles any relative error already there, so an error made early is
multiplied by roughly the remaining power. The total grows roughly in proportion to `n`. #3204
measures up to 66 ulp in the first step (`perlen=10, nstp=100, tsmult=1.005`). How `**` is
evaluated is up to the compiler, so gfortran and ifx/ifort can produce different step lengths
from the same input.

`r**real(n, DP)` calls `pow`. `pow` computes `log(r)` in extra precision internally,
multiplies by `n`, exponentiates, and rounds once, so the result is within 1 ulp. That is #3204's
fix.

### Source 2: cancellation for `tsmult` near 1 (not covered by #3204)

Doubles have constant *relative* precision, so the absolute gap between neighbouring doubles
grows with magnitude:

| Near | Gap between neighbouring doubles |
|---|---|
| 1 | 2.2e-16 (2^-52) |
| 1e-3 | 2.2e-19 (2^-62) |
| 1e-9 | 2.1e-25 (2^-82) |

For `tsmult = 1.0001, nstp = 10`, the quantity we need is `r^n - 1 ≈ 1.0005e-3`. Computed
directly, it could be stored with a gap of about 2e-19, so about 16 significant digits. But the
formula first computes `r^n ≈ 1.0010005`, which is stored with a gap of 2.2e-16. Even correctly
rounded, it can be off by about 1.1e-16. The subtraction `1 - r^n` is then exact (both operands
are within a factor of 2), so it adds no new error. But it can't recover digits that were never
stored, so the 1.1e-16 error survives. Relative to a result of 1e-3 that is about 1e-13, roughly
500 times the precision the result could have had.

In bits: a double near 1 has 52 bits after the binary point. In `1.0010005`, the first 10 of
them are zeros that only pad out the leading "1.00", leaving about 42 bits for the 0.0010005
part. In general, the relative error of `1 - r^n` is about `eps / (2 |1 - r^n|)`, or
`eps / (2 n |h|)` when `n*h` is small. No way of computing the power fixes this, because the
loss happens when a number near 1 is rounded. A real exponent removes only Source 1.

### Fix: `expm1`/`log1p`, with a real exponent outside their range

Use `r^n - 1 = exp(n ln r) - 1` and compute it as `expm1(n * log1p(h))`, so that the number near
1 is never formed:

1. `log1p(h)` returns `ln(1 + h)` with about 1 ulp of relative error, even though the result is
   tiny. Its input `h` is exact. (`log(r)` would do as well here, because `r = 1 + h` exactly and
   good math libraries compute `log` near 1 accurately. The two gave identical results in every
   case tested. `log1p` is the conventional partner of `expm1`.)
2. Multiplying by `n` rounds once, so `y = n * log1p(h)` is within about 1.5 ulp.
3. `expm1(y)` returns `e^y - 1` with about 1 ulp of relative error. It never forms `e^y` near 1
   and subtracts.
4. A relative error δ in `y` becomes about `κ·δ` in `expm1(y)`, where
   `κ = y·e^y / (e^y - 1)`. For small `y`, κ ≈ 1, so step 2's error passes through unamplified.

The denominator is then within a few ulp, the numerator is exact, and the multiply and divide
add one rounding each. The result is a few ulp regardless of `n` and `h`. Both error sources are
gone: there is no chain of multiplications (Source 1), and the small quantity is computed as a
small number rather than as a difference of numbers near 1 (Source 2).

**Limit.** For large positive `y`, κ ≈ `y`. The rounding of `y` (about `eps * y / 2`) is
magnified by the exponential into roughly `y/2` ulp. `pow` avoids this by keeping the logarithm
in extra precision; composing double-precision `log1p` and `expm1` can't. For large negative `y`
(`tsmult < 1`), κ goes to 0 and there is no problem.

**Hybrid.** For `|y| < 1`, use `expm1`/`log1p`. There κ is at most about 1.6, and this is the
range where cancellation matters. Otherwise use `1 - r**real(n, DP)`. There `r^n` is beyond `e`
or below `1/e`, so the subtraction can't cancel badly and `pow` is within 1 ulp:

```fortran
h = tsmult(kper) - DONE
y = real(nstp(kper), DP) * log1p(h)
if (abs(y) < DONE) then
  delt = perlen(kper) * h / expm1(y)
else
  delt = perlen(kper) * (DONE - tsmult(kper)) / &
         (DONE - tsmult(kper)**real(nstp(kper), DP))
end if
```

Fortran has no `expm1`/`log1p` intrinsics. They are available through `bind(c)` interfaces to
the C99 math library, which gfortran, ifx, and the Windows runtime all provide, or as a local
implementation. Like `pow`, they are within about 1 ulp, and results can differ by 1 ulp between
math libraries.

The switch into the geometric branch is an exact `tsmult /= DONE` test, so any `tsmult` a few
ulps from 1 takes the cancelling path at full strength.

### Measurements

All errors are first-step `delt` against exact rational arithmetic for the same double inputs.

Near 1, where cancellation dominates:

| PERLEN | NSTP | TSMULT | `y` | integer `**` | real `**` | `expm1`/`log1p` |
|---|---|---|---|---|---|---|
| 1 | 10 | 1.01 | 0.10 | 8.8 ulp | 5.2 ulp | 0.2 ulp |
| 1 | 10 | 1.001 | 0.010 | 109 ulp | 49.5 ulp | 0.5 ulp |
| 1 | 5 | 1.001 | 0.005 | 150 ulp | 150 ulp | 0.4 ulp |
| 1 | 10 | 1.0001 | 0.001 | 1373 ulp | 226 ulp | 0.2 ulp |
| 1 | 10 | 1.000001 | 1e-5 | 2.4e5 ulp | 7.6e4 ulp | 0.3 ulp |
| 1 | 3 | 1 + 1e-9 | 3e-9 | 6.0e6 ulp | 6.0e6 ulp | 0.3 ulp |
| 100 | 50 | 0.999 | -0.050 | 37 ulp | 5 ulp | 0.0 ulp |

Large `|y|`, where `expm1`/`log1p` alone loses accuracy (`perlen = 1`):

| NSTP | TSMULT | `y` | real `**` | `expm1`/`log1p` |
|---|---|---|---|---|
| 100 | 1.005 | 0.50 | 0.9 ulp | 1.1 ulp |
| 100 | 1.05 | 4.9 | 0.6 ulp | 4.6 ulp |
| 100 | 1.2 | 18 | 0.2 ulp | 6.2 ulp |
| 200 | 1.2 | 36 | 0.7 ulp | 10.7 ulp |
| 1000 | 1.5 | 405 | 0.3 ulp | 44 ulp |
| 100 | 0.95 | -5.1 | 0.6 ulp | 0.6 ulp |
| 500 | 0.9 | -53 | 0.0 ulp | 0.0 ulp |

The cases from #3204:

| PERLEN | NSTP | TSMULT | real `**` | `expm1`/`log1p` | hybrid |
|---|---|---|---|---|---|
| 10 | 100 | 1.005 | 0.7 ulp | 0.3 ulp | 0.3 ulp |
| 1095.75 | 100 | 1.05 | 0.3 ulp | 4.3 ulp | 0.3 ulp |
| 5000 | 100 | 1.02 | 0.3 ulp | 0.7 ulp | 0.3 ulp |
| 365.25 | 6 | 1.3 | 0.9 ulp | 1.9 ulp | 0.9 ulp |

Over 20,000 random cases (`PERLEN` in {1, 10, 365.25, 1e4}, `|tsmult - 1|` from 1e-8 to 0.5
log-uniform, `NSTP` from 2 to 1000), the worst errors were:

| Method | Worst error |
|---|---|
| hybrid | 2.7 ulp |
| `expm1`/`log1p` alone | 367 ulp |
| real `**` alone | 2.1e7 ulp |

### How these were measured, and whether Fortran agrees

The tables come from Python scripts that evaluate the same expressions in double precision
(integer power by repeated squaring, like `__powidf2`; real power and `expm1`/`log1p` from the
platform math library) and compare them with exact rational arithmetic (`fractions.Fraction`).
The scripts are on the `tdis-time-accuracy-analysis` branch in `analysis/tdis-time-accuracy/`:
`first_step_errors.py` for these tables, `step_errors.py` for the middle-step and last-step
figures, and `first_step.f90` with `first_step_cmp.py` for the Fortran comparison. The integer
column also reproduces #3204's numbers exactly (65.7, 23.7, 12.7, and 1.1 ulp).

The same expressions compiled with gfortran 15.3 on arm64 macOS gave **bit-identical** results
for all three forms, at -O0, -O1, and -O3, with inputs both as literals and read at run time. It
also reproduces #3204's 65.7 ulp. Expected differences elsewhere:

- **ifx/ifort** (not tested) may evaluate the integer power differently. That changes the
  integer-`**` column but not the cancellation, which no power method avoids.
- **Other math libraries** (glibc, the Windows runtime) can differ by up to 1 ulp in `pow`,
  `expm1`, and `log1p`, which is small next to the cancellation.
- **FMA contraction** (see the local note `macos-gfortran-fma-findings.md`, not on this branch) doesn't apply: the formula has no
  `a*b + c` pattern. The middle-step and last-step tables above were computed in Python only;
  they use separate multiply and add statements in different routines, and the macOS build has
  no FMA at -O1, but they haven't been checked in Fortran.

## Consumers that make their own time assumptions

Places outside TDIS compute or compare times in ways the options have to account for:

- **`totimc + delt` used as the end of the step** instead of `totim`: `TimeSelect.f90:226`
  (used by PRT release scheduling, PRT output control, and tracking methods) and
  `prt.f90:1125` (particle tracking `tmax`).
  `totimc + delt` is generally not bit-equal to `totim`, and on the retry path's
  end-of-simulation snap it can differ by more. The options should name one canonical
  end-of-step time and switch these sites to it.
- **`TimeSelect.f90` boundary cases** use an absolute `-epsilon(DZERO)` lower bound on the first
  step and `huge` on the last step.
- **BMI `update_until`** (`mf6bmi.f90:131`) treats times within a relative `1e-10` as equal,
  and its comment says this is "so that roundoff accumulated in the simulation time does not
  add a time step". This is evidence that drift is noticeable in practice, and a ceiling the
  measurement should be compared against.
- **Time-series interpolation** (`TimeSeries.f90`) matches breakpoints with `is_close`, whose
  default is relative `DSAME = 100 * eps` (about 50 to 100 ulp of `t`). Drift beyond that can
  change which time-series record a step is matched to. This is the most concrete threshold for
  the measurement step.

## Options

0. **Accurate first step (#3204, plus the `tsmult` near 1 case).** Use the hybrid from [First
   time step length](#first-time-step-length): `expm1`/`log1p` for
   `|nstp * log1p(tsmult - 1)| < 1`, and #3204's real exponent otherwise. Apply it in both
   `tdis_set_delt` and `check_tdis_timing`, through one shared function. This is independent of
   the rest and can land first.
1. **Precomputed period start times.** Compute `perstart(k)` once when TDIS reads its period
   data, using compensated (Kahan or Neumaier) summation of `perlen`, and set `totalsimtime`
   the same way. Each period then starts from an accurate value instead of a carried one. This
   is cheap and local to TDIS.
2. **Closed-form step times for fixed stepping.** Use `totim = perstart(kper) + k*perlen/nstp`
   for `tsmult = 1`, and `perstart(kper) + delt1*(tsmult^k - 1)/(tsmult - 1)` otherwise.
   Snap the last step to `perstart(kper+1)` exactly.
   - The power `tsmult^k` has the same problems as the first step: it needs a real exponent
     and the same hybrid for `tsmult^k - 1` (`expm1(k*log1p(tsmult - 1))` while that argument
     is small), not `tsmult**k`.
   - **Don't derive `delt` as the difference of consecutive times.** That makes every step's
     relative error about `ulp(t) / delt`, which is far worse than compounding late in a long
     simulation (for `t = 1e6` and `delt = 1`, about 5e5 ulp of `delt`). It extends the
     last-step problem above to every step. Compute `delt` from its own closed form instead
     (`delt1 * tsmult^(k-1)`, or `perlen/nstp`), and accept that `totimc + delt` and `totim`
     can differ at the ulp level. See below.
3. **Compensated accumulation for ATS.** Carry a compensation term alongside `totim` (and
   `pertim`). On the last step of a period, snap `totim` to `perstart(kper+1)`. Step retries
   (`tdis_delt_reset`) restore both the value and the compensation from saved copies. The
   snap should apply on both the normal and retry paths, which also fixes the `dtmin`-sized
   loss above.
4. **Snap `totim` to `totalsimtime` on the normal end-of-simulation path.** This is consistent
   with the retry path, but it depends on the open question above. With option 1 or 3 in place,
   this becomes a special case of snapping to `perstart(kper+1)`.
5. **Fix the retry path** independently of accuracy. #3202 resets `endofsimulation` in
   `tdis_delt_reset` and snaps `pertim` at the end of a period. What remains is snapping `totim`
   at the end of every period on both paths (needs `perstart` from option 1, or at least
   `topertim` set correctly for ATS periods).

With options 0 to 3, every `totim` is within a few ulps of the value computed exactly from the
inputs, regardless of step or period count, and every `delt` is within a few ulps of its
nominal value.

### Decoupling `totim` and `delt`

The current scheme keeps `totim(k) = totim(k-1) + delt(k)` up to one rounding per step and
lets the error collect in `totim` (and in the last `delt` of each period). Options 1 to 3 make
`totim` and `delt` each accurate on their own, which means that identity holds only to within
about `ulp(t)`. That is the right trade: `delt` matters for the physics (storage terms) and
`totim` for scheduling and output, and neither should pay for the other. But it requires
that:

- end-of-step comparisons use `totim`, not `totimc + delt` (see consumers above);
- code that reconstructs elapsed time from summed `delt` (budgets, PRT tracking time within a
  step) tolerates ulp-level differences from `totim`.

Not recommended:

- **Integer ticks.** See above; there is no common unit.
- **Quad precision for `totim`.** Compiler support varies, it means converting everywhere
  `totim` is read, and the same benefit is available in double precision with compensated
  summation.

## Impact

`totim`, `totimc`, `pertim`, and `delt` are read widely:

- binary output headers (heads, budgets, concentrations)
- OC and print timing
- observations
- time-series and time-array-series interpolation
- PRT release scheduling and tracking
- the API (`get_current_time`, `update_until`, etc.)
- model-to-model exchanges, which share TDIS, so they stay consistent with each other

Expected effects:

- **Ulp-level changes in reported times and in `delt`.** `delt` feeds storage terms, so heads
  and concentrations could also change at roundoff level. The last step of each geometric
  period would change by more (thousands of ulp for `tsmult < 1`), though still around 1e-12
  relative.
- **Comparisons.** Most autotest comparisons use tolerances, but exact comparisons of times
  (e.g. headers in binary file comparisons, snapshot tests, `==` on times in Python tests) may
  need updating. Any regenerated snapshots should be listed and justified.
- **Time-based branching.** Code that branches on time near boundaries (time-series breakpoints,
  OC/observation selection, PRT release windows `(totimc, totimc+delt]`) could change which
  side of a boundary a value falls on. Accurate times should make these boundaries behave
  *more* predictably, but individual results could still change.
- **ATS retry fixes (option 5)** change behavior beyond roundoff where they apply: period start
  times shift by up to `dtmin`. (#3202 already changes results for simulations that stopped
  early, which now run to the end.)
- **Users.** The changes are below any meaningful precision for most models. Models with many
  steps, long simulation times, or `tsmult != 1` would see the largest differences. A release
  note is needed.

## Preliminary step: measure drift (not yet done)

Before committing to the change, measure how much error the current scheme actually accumulates:

1. **Instrument TDIS in a throwaway build.** For each step, log `totim` and `delt` alongside:
   - a reference computed with compensated summation (or in `real(REAL128)` where the compiler
     supports it), and
   - the closed-form value for fixed-stepping periods.

   Report the difference in ulps and as a relative error. Report `delt` error separately for
   the last step of each period. No other behavior should change.
2. **Run a representative set of models,** chosen to stress each mechanism:
   - many steps per period with `tsmult = 1`
   - `tsmult > 1` with large `NSTP` (first-step error, #3204)
   - `tsmult < 1` (last-step cancellation)
   - `tsmult` within about 1e-3 of 1 (first-step cancellation)
   - many periods
   - ATS with many steps and retries, including a retry on the last step of a period and of
     the simulation
   - a long simulation in seconds (T around 1e9 or more)

   Draw these from autotest and the example models. Build with both gfortran and ifx, since the
   first step is compiler-dependent.
3. **Summarize.** For each model, report the maximum and final drift in ulps and the step at
   which drift crosses common tolerances (`DSAME * t` as used by time-series matching, the BMI
   `1e-10 * t`, and the PRT release-time tolerance). Also note any model where drift is large
   enough to flip a time comparison.

If measured drift is consistently within a few ulps, the case for options 1 to 3 is weak and the
current scheme can be documented as-is. If drift reaches hundreds or thousands of ulps in
realistic models, that supports doing them. Option 0 and option 5 are worth doing either way.

## Relationship to #3204

[#3204](https://github.com/MODFLOW-ORG/modflow6/issues/3204) is option 0 without the `tsmult`
near 1 case:

- It removes the integer-exponent error in the first step and the compiler dependence. That is
  the main source of `totim` error within geometric periods with `tsmult > 1`, and of most of
  the last-step `delt` error there.
- It doesn't address accumulation (compounding, carry-over across periods, ATS `totim`,
  `totalsimtime`), cancellation for `tsmult` near 1, or the last-step cancellation for
  `tsmult < 1`.
- Its autotest survey is a first data point for the measurement step: of 105 simulations with
  `NSTP > 1` and `TSMULT != 1`, 15 of 17 distinct period settings get a different first step, by
  at most 1.2e-14 relative (roughly 50 to 100 ulp). Period end times are unchanged.
- It trades exact cross-compiler agreement in the algorithm for platform `pow`, which can differ
  by 1 ulp between math libraries. Option 2 inherits the same trade-off; `expm1`/`log1p` has it
  too.
- Landing it first makes later drift measurements cleaner, since first-step error would no
  longer be mixed in. The hybrid keeps #3204's real exponent for `|y| >= 1` and only changes
  periods with `tsmult` near 1. It's best done together with #3204, so outputs change only once.
- The cancellation case hasn't been raised on #3204 yet.

## Relationship to the PRT change

The PRT release-time tolerance change (local note `prt-release-time-tolerance-plan.md`, not on this branch) doesn't depend on
this work. It keeps an explicit floating-point floor and derives the default tolerance from
simulation length. If TDIS time accuracy improves, the PRT default could later shrink to a small
relative tolerance, and possibly the option could be retired. PRT also uses `totimc + delt` as
the end of the step (`prt.f90:1125`), which the decoupling above affects.

## Open questions

- Is skipping the `totim = totalsimtime` snap on the normal end-of-simulation path intentional?
  Or was the retry-path snap added to cover the `dtmin`-sized shortfall described above?
- Is `topertim` intentionally left at zero for ATS periods?
- Do any downstream tools (FloPy, the API consumers, the PRT/MODPATH comparison tests) depend
  on the current exact `totim` values?

## Unrelated issue found during the audit

- `check_tdis_timing` returns as soon as it finds a period with `perlen == 0` (line 639,
  after the warning). Later periods are then not validated at all, so a negative `PERLEN`,
  `NSTP <= 0`, or `TSMULT <= 0` in a later period is not reported. Confirmed: negative `PERLEN`
  and negative `TSMULT` run to normal termination, and `NSTP = 0` crashes. See
  `zero-perlen-validation.md`.
