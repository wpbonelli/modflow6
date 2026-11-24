# XT3D across discontinuities: #3104, #3208, #415

Working notes and plan. Not intended for merge as-is.

## Summary

Three open issues touch the same corner of XT3D: how a connection's flux is
computed when a sharp discontinuity (K contrast or HFB) is nearby.

| Issue | Symptom | Cause | Status |
|---|---|---|---|
| [#3104](https://github.com/MODFLOW-ORG/modflow6/issues/3104) | Discharge wrong by 16–107% across a K contrast on refined DISU, isotropic K | Secondary gradients estimated naively across the contrast | Provost prototyping a "heterogeneity correction" |
| [#3208](https://github.com/MODFLOW-ORG/modflow6/issues/3208) | Head overshoot / circulation next to an HFB with rotated anisotropy | Same as #3104; the HFB is the discontinuity | Open |
| [#415](https://github.com/MODFLOW-ORG/modflow6/issues/415) | HFB multiplier option (`hydchr <= 0`) wrong under XT3D | Wrong formula for the primary-face correction | Open since 2020, still present |

#3104 and #3208 are the same defect showing up in different ways. #415 is a
separate formula bug in the same routine, but its second point (what a
multiplier means for XT3D's multi-term flux) has to be settled for any
#3208 fix.

Not related: #2826 (dispersion ellipse orientation, fixed in 6.8.1 by #2849)
was an `acos`/`atan2` sign error in `calcdispellipse`, upstream of XT3D.

## Mechanism (#3104 and #3208)

XT3D computes the flux across a face (the *primary* connection) from

- the head difference across that face (`ahat` / `chat01`), plus
- an along-face gradient estimated from the cell's other (*secondary*)
  connections (`bhat` / `chati0`, `chat1j`), weighted by the off-axis
  terms `sigma(2)`, `sigma(3)` of `vn . K . rmat`.

`abhats` (`src/Model/ModelUtilities/Xt3dAlgorithm.f90`) estimates each
secondary gradient as `(h_nbr - h_cell) / (dl0 + dln)`, i.e. the plain head
difference over the full connection length. A K contrast or barrier on that
secondary connection is ignored, so a sharp head jump is read as a steep
smooth gradient and leaks into the primary flux.

Two ingredients are needed together:

| | Off-axis `sigma` terms from | Discontinuity |
|---|---|---|
| #3104 | grid geometry (sub-vertical connections from 1/4/1 refinement) | 100x K contrast |
| #3208 | rotated anisotropy (`angle1` = 30/45 deg) on a regular grid | HFB |

The HFB correction (`xt3d_fhfb`, `xt3d_flowjahfb` in
`src/Model/ModelUtilities/Xt3dInterface.f90`) only scales the barrier face's
own coefficients. Secondary connections that cross the barrier are not
touched, and `xt3d_load` builds the K tensor from NPF arrays only, so the
barrier is invisible to them.

How the error shows up differs:

- #3104: the bad gradient lands on the primary face, so net discharge is
  wrong.
- #3208: net flow through the barrier is correct (~0); the error is
  circulation parallel to the barrier on both sides, seen as overshoot.

## Mechanism (#415)

Standard HFB, `hydchr <= 0`: `cond = |hydchr| * cond`
(`gwf-hfb.f90`, in `hfb_fc` and `condsat_modify`).

XT3D: `term = -condhfb`, then each coefficient is reduced by `term * c`,
leaving `(1 - |hydchr|) * c`. Same branch in `xt3d_fhfb` and
`xt3d_flowjahfb`.

## Evidence

Scripts in this directory. Run from the repo root with a built `bin/mf6`:

```
pixi run python xt3d-discontinuity/exp_hfb.py $PWD/bin/mf6 <scratch dir>
pixi run python xt3d-discontinuity/exp3104.py $PWD/bin/mf6 <scratch dir>
```

Results below from `6.9.0.dev0` (develop at 68e79f9ed9).

### #3208 (20x20, K=100, K22=1, CHD 1 | 0, barrier between cols 10/11)

`net` = net flow across the barrier line; `circ` = gross inflow at the
right-hand CHD (should be ~0; measures spurious circulation).

| Case | hmax | hmin | net | circ |
|---|---|---|---|---|
| aniso 45, HFB 1e-12, XT3D | 1.284 | -0.284 | 3.5e-12 | 14.6 |
| same, XT3D off | 1.000 | 0.000 | 2.0e-11 | 3.5e-10 |
| isotropic, HFB, XT3D | 1.000 | 0.000 | 2.0e-11 | 2.4e-11 |
| angle 0, HFB, XT3D | 1.000 | 0.000 | 2.0e-11 | 1.1e-09 |
| aniso 45, no HFB, XT3D | 1.014 | -0.014 | 18.4 | 0.73 |

Overshoot grows smoothly with barrier strength (`hydchr` 1e3 to 1e-12:
hmax 1.014 to 1.284).

Barrier replaced by a column of low-K cells (no HFB):

| Column width | hmax (XT3D) | hmax (std) |
|---|---|---|
| 1.0 | 1.158 | 1.000 |
| 0.1 | 1.327 | 1.000 |
| 0.01 | 1.362 | 1.000 |

Isotropic low-K column with XT3D: no overshoot. So the HFB code is not the
cause; any sharp discontinuity plus off-axis terms triggers it.

The aniso 45 / no-HFB overshoot (1.014) is the usual non-monotonicity of
MPFA-type schemes with 100:1 rotated anisotropy, and is a separate matter.

### #3104 (six-cell 1/4/1 DISU)

| K middle layer | Std error | XT3D error |
|---|---|---|
| 0.01 | 0.00% | +107.32% |
| 0.1 | 0.00% | +76.09% |
| 1.0 | 0.00% | 0.00% |

### #415 (isotropic regular grid, so XT3D should equal std)

| hydchr | Intended | Std net | XT3D net |
|---|---|---|---|
| -0.1 | x0.1 | 71.4 | 104.7 |
| -0.5 | x0.5 | 100.0 | 100.0 (coincidence) |
| -1.0 | x1 | 105.3 | ~0 |
| -1.5 | x1.5 | 107.1 | does not converge |
| 0 | blocked | ~0 | 105.3 |

## Plan

### 1. Fix #415 now (small, independent)

- In `xt3d_fhfb` and `xt3d_flowjahfb`, change the `condhfb <= 0` branch to
  `term = 1 + condhfb` (i.e. `1 - |m|`), so the corrected coefficients are
  `|m| * c`.
- Check the Newton path in `xt3d_fhfb` (qsat accumulation) is consistent.
  Its comment assumes the correction cannot reverse flow, which holds for
  `0 <= |m| <= 1` but not for `|m| > 1`; decide whether multipliers above 1
  are allowed under XT3D.
- Test: isotropic regular grid, multipliers {0, -0.1, -0.5, -1, -1.5},
  assert XT3D heads and FLOW-JA-FACE match the standard formulation.
- Release note under fixes.

### 2. Coordinate #3104 / #3208 with the heterogeneity correction

Owner: Provost's prototype from #3104. Questions to put to him and Huite:

- Does the correction work purely from NPF K, or can it accept an extra
  resistance per connection? For #3208 it must include HFB resistance
  (`1/hydchr`) on secondary connections, not only the primary face.
- Huite's proposals on #3208 (correct the slope for the barrier jump;
  weight links by barrier factor so they vanish as resistance -> inf) look
  like the same correction specialised to barriers. Agree on one approach.
- Limiting behaviour: as resistance -> inf, connections across the barrier
  should drop out; as resistance -> 0, recover current XT3D.

### 3. Decide the multiplier option under XT3D (#415 point 2)

A barrier-aware secondary gradient needs a resistance. A multiplier gives
none directly. Options:

- Convert to an equivalent resistance from the unmodified standard
  conductance of that connection, then treat like a positive `hydchr`.
- Disallow `hydchr <= 0` with XT3D and error at read time.

### 4. Tests (with step 2)

- #3208 model with HFB: heads within [0, 1] for the strong-barrier case;
  net barrier flow ~0; circulation ~0.
- Same with a thin low-K column instead of HFB: as width -> 0 the two should
  converge (Huite's consistency check).
- #3104 six-cell model: discharge within a tolerance of analytical.
- Regression: existing XT3D tests, and Provost et al. (2025) dipping
  aquifer cases.

## Open questions

- Does the heterogeneity correction change results for models with no
  discontinuity? If so, does it need an option flag?
- Parallel / GWF-GWF exchanges: secondary connections across a model
  boundary, and HFBs on exchanges, need the same treatment.
- Cost: barrier lookup per secondary connection in `xt3d_load`.
