# Ask the discretization for geometry

**Proposal:** model packages should get cell and connection geometry by *asking the
discretization object* (type-bound accessors on `DisBaseType`, overridden per grid type)
instead of reading the raw `dis%top`, `dis%bot` and `dis%area` arrays and re-deriving
volumes, thicknesses and flow areas themselves.

The idea came up while reviewing the fully 3D DISU proposal
([#3207](https://github.com/MODFLOW-ORG/modflow6/issues/3207); see
`full3d-disu-feasibility.md` on branch `full3d-disu-feasibility`). This doc argues it is
worth doing on its own merits, regardless of what happens with that proposal.

Counts and line references below are from `upstream/develop` at the time of writing.

## Executive summary

`DisBaseType` exposes `top`, `bot`, and `area` as public arrays, and packages use them to
re-derive cell volumes, thicknesses, and connection flow areas themselves. There are ~250
such references across 37 source files. The proposal is to route derived geometry through
type-bound accessors on the discretization (with per-grid-type overrides), so each
geometric concept has one implementation owned by the grid. The pattern already exists
and works: SWF calls `dis%get_flow_width`, which DISV1D overrides, so SWF packages never
need to know the grid type. This would extend it to the quantities that matter most
numerically.

### Motivation

- **Duplication.** `get_cell_volume` exists but has only 2 callers. 17 sites in GWT
  MST/IST, GWE EST, and FMI write out `area*(top-bot)` inline. The saturated flow area of
  a connection (vertical → `hwva`; horizontal → `hwva` × saturated thickness, with
  staggered/convertible special cases) is reimplemented in NPF, XT3D, DSP, CND,
  `GwfGwfConnection`, the GWF-GWF exchange, and the central-difference weights.
  `calcdispcoef` in DSP and CND is 127 lines each, and only 3 differ. Every copy is a
  chance for a fix to miss the others.
- **Encapsulation.** What volume, thickness, and flow area mean depends on the grid type
  (prisms, DIS2D/DISV1D, the synthetic interface-model DISU). With raw-array access every
  package implicitly hard-codes those conventions, so adding or changing a grid type means
  auditing all of them.
- **Less orientation branching.** About 18 files branch on `ihc`. Many of those branches
  only choose between "area" and "width × thickness", and would collapse into a single
  accessor.
- **Testability.** Accessors can be unit-tested per grid type (staggered connections,
  dry/pinched cells, `idomain` reduction) rather than only through integration tests.
- **Parallel/interface models.** The interface DISU currently fills `top`/`bot`/`area`
  with placeholders and overwrites them from virtual model data, because consumers read
  the raw arrays. Accessors make the required contract explicit and narrow what has to be
  transferred. HFB also bypasses the dis object and fetches `TOP`/`BOT` via memory paths.
- **Room to evolve.** Non-prism or polyhedral cells (cf. #3207), time-varying thickness,
  alternative saturation functions, and capability flags (like the existing
  `supports_layers`) become local to the discretization instead of spreading `if`
  branches across packages.

### Sketch

Default implementations in `DisBaseType` reproduce today's prism behavior exactly:

- cells: `get_cell_volume(n)`, `get_sat_volume(n, h)` / `get_dvdh(n, h)`,
  `get_cell_thickness(n)`, `get_sat_thickness(n, h)`, `get_cell_center(n, x, y, z)`
- connections: `get_flow_area(ipos, satn, satm)` alongside the existing `get_flow_width`
  and `connection_normal/vector`
- capabilities: e.g. `has_vertical_extent`, `supports_unconfined`

Real elevation comparisons (boundary elevations against the cell bottom, highest active
cell, LAK/MAW/SFR connection geometry) can keep using elevations. The goal is to stop
packages *re-deriving* geometry, not to ban `top`/`bot`.

### Plan

Each step should be behavior-preserving (identical results to round-off) and mergeable on
its own:

1. Replace the 17 inline cell-volume formulas with an accessor.
2. Add `get_flow_area`, and merge the duplicated DSP/CND coefficient code.
3. Route NPF/XT3D/exchange/advection flow-area logic through it, collapsing the `ihc`
   branching. This is the riskiest step and needs review against the staggered, Newton,
   and convertible-cell models.
4. HFB: use the dis object instead of memory-path pointers.
5. Migrate the remaining thickness/center derivations (STO, CSUB, advanced packages).
6. Add capability flags.
7. Narrow the virtual-model data the interface model needs.

Steps 1–2 are worthwhile even if we go no further.

### Concerns

- **Performance:** a virtual call per cell or connection in `fc`/`cq` loops. Mitigate by
  caching static quantities in `dis_ar`, adding array-valued variants for hot loops, and
  benchmarking large models before and after.
- **Churn:** keep PRs small and component-scoped.
- **Intentional differences** between existing copies must be preserved explicitly
  through arguments rather than merged away.

The rest of this document gives the evidence and detail behind each point.

## The current situation

`DisBaseType` exposes geometry as public pointer arrays:

```fortran
real(DP), dimension(:), pointer, contiguous :: top  !< (size:nodes) cell top elevation
real(DP), dimension(:), pointer, contiguous :: bot  !< (size:nodes) cell bottom elevation
real(DP), dimension(:), pointer, contiguous :: area !< (size:nodes) cell area, in plan view
```

The model and package code uses them directly:

- **252 references** to `dis%top`, `dis%bot` or `dis%area` in **37 source files**: GWF 15,
  SWF 5, GWT/GWE/TSP 9, PRT 3, exchanges 3, plus a few others.
- **HFB bypasses the discretization object entirely.** It fetches `TOP`/`BOT` from the
  memory manager by path (`gwf-hfb.f90:127-128`).
- **The parallel virtual model ships `TOP`, `BOT` and `AREA` across processes**
  (`VirtualModel.f90:120-122`), because the interface model's synthetic DISU has to
  reproduce the raw arrays that downstream code reads.

The base class already has some accessors, but they're used inconsistently:

| Accessor | Exists | Used by | Bypassed by |
|---|---|---|---|
| `get_cell_volume(n, x)` | yes | 2 call sites (`gwf.f90:675`, `tsp-adv.f90:203`) | 17 inline `area*(top-bot)` in GWT MST (8), GWE EST (6), GWT IST (2), FMI (1) |
| `get_area(n)` | yes | RCH, EVT, CSUB, SWF packages | direct `dis%area(n)` elsewhere |
| `get_flow_width(n, m, ...)` | yes; **overridden by DISV1D** | SWF DFW/STO/EVP/PCP, SWF-GWF exchange | n/a (this one works as intended) |
| `connection_normal/vector` | yes | NPF, XT3D, DSP/CND, exchanges | n/a |
| `supports_layers`, `is_3d/2d/1d` | yes | capability checks | n/a |

So the pattern already exists, and SWF shows it working: `get_flow_width` lets DISV1D
supply its own per-cell widths without SWF packages knowing the grid type. What's
missing is applying it consistently, especially to the two things that matter most
numerically: **cell volume/saturation** and **connection flow area**.

## Motivation (independent of 3D grids)

### 1. The same physics is written out many times

- **Cell volume.** `get_cell_volume` exists, but 17 call sites in transport and energy
  re-derive `area(n) * (top(n) - bot(n))` inline.
- **Saturated flow area of a connection.** "Vertical → `hwva`; horizontal → `hwva` ×
  saturated thickness, with a special case for staggered connections and for convertible
  or unsaturated cells" is re-implemented in NPF (via `GwfConductanceUtils`), XT3D
  (`xt3d_areas`), GWT DSP, GWE CND, GWF-GWF spdis edge properties (`GwfGwfConnection.f90`),
  the GWF-GWF exchange conductance, and the central-difference advection weights.
- **Thickness `top - bot`.** Computed inline in at least 25 files (12 in
  `GwfConductanceUtils.f90` alone, 8 in `gwt-mst.f90`, 7 in `gwf-csub.f90`, …).
- **DSP and CND are near copies.** `calcdispcoef` in `gwt-dsp.f90` and `gwe-cnd.f90` is
  127 lines in both, and only 3 lines differ.

Every copy is a chance for the copies to drift apart: a fix to how a dry or staggered cell
contributes flow area in NPF does not reach DSP/CND/XT3D unless someone remembers. One
implementation per concept means one place to get it right, one place to review, and one
place to document.

### 2. The discretization should own its geometric conventions

What "volume", "thickness" and "flow area" mean depends on the grid type:

- DIS/DISV/DISU: prisms with plan-view `area`, `top` and `bot`.
- DIS2D/DISV1D (SWF): `top`/`area` have different or partial meaning, and widths can
  differ per cell, which is why `get_flow_width` is overridden.
- Interface models (parallel / XT3D exchanges): a synthetic DISU whose geometry is
  assembled from neighboring models.

When packages read raw arrays, each package implicitly takes a position on these
conventions, and adding or changing a grid type means auditing every package. When
packages call `dis%get_cell_volume(n)` or `dis%get_flow_area(ipos, ...)`, the grid type
answers for itself, and packages don't need to know which grid they're on. That's the
same encapsulation the SWF code already relies on.

### 3. Grid-specific logic stops spreading through package code

Branches on connection orientation (`ihc == 0`, `ihc == 2`, …) appear in about 18 files,
including HFB (8), NPF (6), the GWF-GWF exchange (6), XT3D (5) and BUY (5). Some are
legitimate physics (buoyancy really cares about vertical), but many exist only to choose
between "area" and "width × thickness". Behind a flow-area accessor those branches collapse
into one implementation in the discretization. Any new connection type then gets added in
one place, not in every consumer.

### 4. Geometry becomes unit-testable

Geometry math embedded inside `calcdispcoef` or `sto_fc` can only be tested end to end
through integration tests. Accessors on the discretization can be tested directly with
the Fortran unit-test suite (`autotest/Test*.f90`): build a small DIS/DISV/DISU, check
volumes, saturated volumes, flow areas and normals against hand calculations, including
the edge cases (staggered connections, dry cells, pinched cells, `idomain` reductions)
that are currently covered only indirectly.

### 5. Parallel and interface models get simpler and safer

The interface model's DISU is built in `GridConnection%getDiscretization` by filling
`top`/`bot`/`area` with placeholders (`-huge`) and later overwriting them from virtual
model data. It has to look exactly like a real DISU at the raw-array level, because that's
what consumers read. If consumers call accessors instead, the contract becomes explicit:
the interface DISU must answer these specific queries correctly. The virtual model can
then send exactly the data those answers need. Any shortfall shows up as one unimplemented
accessor instead of a wrong number somewhere downstream.

### 6. Room to evolve the geometry model without touching every package

Changes that today would mean editing dozens of packages become local to the
discretization:

- **Different cell shapes:** non-prismatic or polyhedral cells (the #3207 case), or
  sloping/deformed layers.
- **Time-varying geometry:** for example, cell thickness updated by compaction (CSUB
  already tracks it separately).
- **Alternative saturation functions:** the quadratic smoothing in
  `get_cell_volume`/`sQuadraticSaturation` lives in one place instead of being reapplied
  inconsistently.
- **Alternative flow-area definitions:** for example, cross-sectional geometry in SWF,
  which is already partly behind `get_flow_width` and `swf-cxs`.

None of these is on the roadmap as a commitment. The point is that the current structure
makes each of them expensive, and the proposed one makes them cheap.

### 7. Capability checks happen once, in one place

`supports_layers` and `is_3d/2d/1d` already show the pattern. With geometry behind
accessors, a grid type that can't answer a question (for example, saturated volume for a
confined-only grid) can say so with one capability flag that packages check once at
allocate-and-read time. The alternative, which the #3207 fork illustrates, is about 60
grid-specific `if` branches spread over 14 files.

## Sketch of the API

All accessors take **reduced** node numbers / `ia/ja` positions, consistent with the
existing methods. Default implementations in `DisBaseType` reproduce today's prism
behavior exactly, so the migration changes no results.

**Cell queries**

| Accessor | Prism default | Replaces |
|---|---|---|
| `get_cell_volume(n)` (no head argument) | `area * (top - bot)` | 17 inline formulas |
| `get_cell_volume(n, x)` (exists) / `get_sat_volume(n, h)` | `area * (top - bot) * sat(h)` | — |
| `get_dvdh(n, h)` | `area * dsat/dh * (top - bot)` | ad hoc derivatives in STO/Newton |
| `get_cell_thickness(n)`, `get_sat_thickness(n, h)` | `top - bot`, `(top - bot) * sat` | inline `top - bot` |
| `get_cell_center(n, x, y, z)` | `(xc, yc, (top + bot) / 2)` | scattered center calcs |
| `get_area(n)` (exists) | plan area | direct `dis%area(n)` |
| `get_top(n)`, `get_bot(n)` | `top(n)`, `bot(n)` | elevation comparisons (see below) |

**Connection queries**

| Accessor | Prism default | Replaces |
|---|---|---|
| `get_flow_area(ipos, satn, satm)` | `hwva` (vertical) or `hwva * thksat` (horizontal / staggered) | NPF, XT3D, DSP, CND, GwfGwfConnection, exchange copies |
| `get_flow_width(...)` (exists) | `hwva` | — |
| `connection_normal/vector` (exist) | unchanged | — |

**Capabilities**

`supports_layers` (exists), plus for example `has_vertical_extent` and
`supports_unconfined`.

Not every raw reference should go. Some uses really are elevation comparisons: validating
a boundary elevation against the cell bottom, RCH/EVT finding the highest active cell,
LAK/MAW/SFR connection geometry, PRT release-point checks. Those can keep reading
elevations through `get_top`/`get_bot` (or the raw arrays, at first). The target is to
remove *re-derived geometry* (volumes, thicknesses, flow areas) from packages, not to ban
elevations.

## Performance

The obvious objection: a type-bound (virtual) call per cell or per connection inside
`fc`/`cq` loops costs more than an array read, and gfortran won't inline across a
polymorphic dispatch without LTO and devirtualization. Mitigations:

- **Cache what is static.** Full cell volume, thickness and fully saturated flow area don't
  change during a run. Compute them once in `dis_ar` (or lazily) into arrays owned by the
  discretization. The accessor then returns an array element, or the caller grabs the
  cached array once per loop.
- **Dispatch once per loop, not per element.** Provide array-valued variants where a hot
  loop needs them, for example `call dis%get_sat_volumes(h, vol)`.
- **Measure.** Benchmark a few large models (e.g., from `modflow6-examples` or the
  large-test repo) before and after each step. If an accessor shows up in profiles, cache
  it. Most per-cell geometry work is negligible next to the linear solve.

## Migration plan

Each step is behavior-preserving and independently mergeable. Verify each with the
affected unit and integration tests, expecting results identical to round-off.

1. **Cell volume.** Replace the 17 inline `area*(top-bot)` with a no-head
   `get_cell_volume(n)` (or a cached volume array). Mechanical and low risk.
2. **Flow area.** Add `get_flow_area`, then migrate DSP and CND to it and to each other:
   fold the duplicated `calcdispcoef` into a shared routine parameterized by the
   coefficient arrays.
3. **NPF/XT3D/exchanges.** Route the saturated-flow-area logic in NPF, `xt3d_areas`,
   `GwfGwfConnection`, the GWF-GWF exchange and central-difference weights through the same
   accessor. This is the step where the `ihc` branching collapses. It's also the riskiest
   step, so it needs careful review against the staggered/Newton/convertible test models.
4. **HFB.** Stop pulling `TOP`/`BOT` via memory paths and use the discretization object.
5. **Thickness/center accessors.** Migrate the remaining derived-geometry call sites
   (STO, CSUB, LAK/MAW/SFR where they compute thickness rather than compare elevations).
6. **Capabilities.** Add the capability flags. New grid types (for example a polyhedral
   discretization) then override accessors and flags, not packages.
7. **Interface model.** Narrow the virtual-model data to what the accessors need.

Steps 1–2 are useful even if nothing else happens.

## Risks and costs

- **Churn.** Many files are touched, so there will be conflicts with open branches. Doing
  it per step and per component keeps each PR reviewable.
- **Hidden semantics.** Some inline formulas may differ subtly on purpose (for example,
  which cell's saturation is used for a vertical connection between convertible cells).
  Unifying them must preserve those differences explicitly via arguments, not erase them.
  The identical-results requirement catches accidental changes.
- **Over-abstraction.** The API should cover what packages actually compute. It should not
  become a general geometry library, and elevations that packages genuinely need can stay
  as elevations.
- **Performance.** Covered above: cache and measure.

## Relationship to #3207

With this in place, the fully 3D proposal becomes a new discretization type that answers
these queries for polyhedra (volume from input or vertices, flow area = face area, normal
from geometry), declares `supports_unconfined = .false.`, and needs no changes to the
transport and storage packages. The fork's placeholder trick (`top=1`, `bot=0`,
`area=volume`), its 60 `is_full_3d()` branches, and its silent dependence on thickness
being 1 all disappear. But as argued above, the refactor stands on its own whether or not
#3207 goes ahead.
