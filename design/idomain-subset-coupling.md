# Design: per-model active domains for GWF-coupled models

- Status: PRT prototype validated; GWT prototype functional but has a known
  numerical caveat at domain-subset boundaries (see "Prototype status:
  GWT" below) -- not ready to propose upstream as-is.
- Related: [discussion #2420](https://github.com/MODFLOW-ORG/modflow6/discussions/2420),
  [issue #1129](https://github.com/MODFLOW-ORG/modflow6/issues/1129)
- Author: wpbonelli (drafted with Claude)
- Date: 2026-08-18

## Problem

MT3DMS had `ICBUND`, which let a user turn off transport in parts of the grid
without changing the flow solution there. MODFLOW 6 has no equivalent: GWT,
GWE, and PRT models coupled to a GWF model are currently required to have an
**identical** `IDOMAIN` array to that GWF model. A user who wants to exclude
part of the flow domain from transport/energy/particle-tracking today has to
zero out all transport parameters, sources, and initial concentrations in
that region instead — which (per `aprovost-usgs` in the discussion) doesn't
meaningfully reduce runtime, since the transport calculation still executes
there.

Issue #1129 asks for GWT's active domain to be allowed to be *smaller* than
GWF's, and notes GWE and PRT would likely want the same. This doc proposes a
concrete mechanism, scoped to keep the first implementation tractable.

## Current constraint, in code

Each model type owns its own DIS object and its own independent `idomain`
array (not a shared pointer):

- `src/Model/Discretization/Dis.f90:31` — structured `idomain(ncol,nrow,nlay)`
- `src/Model/Discretization/Disu.f90:47` — `idomain(nodes)`
- `src/Model/Discretization/Disv.f90` — analogous, DISV
- `src/Model/ModelUtilities/DiscretizationBase.f90:59-60` — `nodereduced` /
  `nodeuser`, the per-model map between "user" node numbers (dense, includes
  inactive cells) and "reduced" node numbers (compressed, active cells only).
  Comment: "0 is idomain = 0". This mapping is derived independently by each
  model from its own `idomain`.

The GWF-GWT, GWF-GWE, and GWF-PRT exchanges each hard-require the two
models' `idomain` arrays to be identical, in an identical block of code:

- `src/Exchange/exg-gwfgwt.f90:223-260`
- `src/Exchange/exg-gwfgwe.f90:195-260`
- `src/Exchange/exg-gwfprt.f90:191-252`

```fortran
! -- Check to make sure sizes are identical
if (gwtmodel%dis%nodes /= gwfmodel%dis%nodes .or. &
    gwtmodel%dis%nodesuser /= gwfmodel%dis%nodesuser) then
  ! ... store_error(..., terminate=.TRUE.)
end if

! -- Make sure idomains are identical
select type (gwfdis => gwfmodel%dis)
type is (DisType)
  select type (gwtdis => gwtmodel%dis)
  type is (DisType)
    if (.not. all(gwfdis%idomain == gwtdis%idomain)) then
      ! ... "Ensure discretization packages, including IDOMAIN, are identical."
    end if
  end select
! ... same pattern for DisvType, DisuType
end select
```

**Why it's this strict today:** immediately after this check, the exchange
wires up direct pointers from GWF arrays into the transport/particle
model's FMI object:

```fortran
gwtmodel%fmi%gwfhead => gwfmodel%x
```

`tsp-fmi.f90` (shared by GWT/GWE) and `prt-fmi.f90` then index those pointers
directly using the *transport model's own* reduced node numbers
(`tsp-fmi.f90:439-488`, `prt-fmi.f90:118-142,254`), e.g.:

```fortran
do n = 1, this%dis%nodes
  if (this%gwfsat(n) > DZERO) then           ! gwfsat indexed by GWT's own node number
    this%ibdgwfsat0(n) = 1
  ...
```

and, for connection-indexed flow data, using the *transport model's own*
`con%ia`/`con%ja` positions to index straight into `gwfflowja`
(`tsp-fmi.f90:462-494`):

```fortran
do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
  m = this%dis%con%ja(ipos)
  flownm = this%gwfflowja(ipos)   ! ipos is a GWT connection position, used to index GWF's flowja
```

This only produces correct results if GWT's reduced node numbering *and*
connection ordering are bit-for-bit identical to GWF's — which is only
guaranteed when `idomain` is identical. That's the real reason for today's
restriction; it isn't fundamental to the physics.

**Existing precedent for relaxing whole-grid identity:** GWT-GWT and GWE-GWE
exchanges already avoid a shared-grid requirement. `exg-gwfgwt.f90:370-408`
links a GWT-GWT exchange to its parent GWF-GWF exchange purely by comparing
explicit node lists (`nodem1`/`nodem2`), not by requiring a shared DIS
object or matching whole-grid arrays. `DiscretizationBase.f90` already
exposes the primitives needed to translate between two models' reduced
numbering:

- `get_nodeuser(noder)` — reduced node → user node (`DiscretizationBase.f90:272-282`)
- `get_nodenumber(nodeu, icheck)` — user node → reduced node
  (`DiscretizationBase.f90:69-80,284-315`)

## Proposed scope

Support the case **child domain ⊆ parent (GWF) domain**: a coupled
GWT/GWE/PRT model may mark inactive (`idomain = 0`) any cell that is active
in GWF, but may not mark active a cell that GWF has marked inactive. This
is exactly the ICBUND use case in the original discussion, and exactly what
issue #1129 asks for ("smaller than the GWF active domain").

Explicitly **out of scope** for this iteration:

- A transport/energy/particle domain larger than or disjoint from GWF's
  active domain. That requires either extrapolating flow data into cells
  GWF never solved, or genuinely different grids — a much bigger lift with
  no concrete ask behind it yet.
- Independent grids/geometries between coupled models (different NCOL/NROW/
  NLAY, different DISU node count, etc). We still require the same DIS
  subtype and the same user-node addressing (`nodesuser` and grid shape
  match); only the *active/inactive flags* may differ.

## Design

### 1. Relax the exchange check

Replace the equality test with a subset test, in all three exchange types
(`exg-gwfgwt.f90`, `exg-gwfgwe.f90`, `exg-gwfprt.f90`):

```fortran
! child idomain must be a subset of parent idomain
if (any(gwtdis%idomain /= 0 .and. gwfdis%idomain == 0)) then
  ! store_error: "<model> active domain is not contained within GWF active domain"
end if
```

Grid-shape / `nodesuser` / DIS-subtype checks stay as-is — those are still
required so that user-node addressing lines up between the two models.

### 2. Build a node map at `exg_ar` time

One pass, O(nodesuser), computed once when the exchange is set up:

```fortran
do n = 1, childmodel%dis%nodes         ! child's reduced node numbering
  nodeu = childmodel%dis%get_nodeuser(n)
  node_map(n) = gwfmodel%dis%get_nodenumber(nodeu, 0)   ! GWF's reduced node number, or 0/less if inactive in GWF (shouldn't happen given the subset check)
end do
```

Store `node_map` on the FMI object (`TspFmiType`, `PrtFmiType`). Every place
that currently does `this%gwfhead(n)`, `this%gwfsat(n)`, `this%gwfspdis(:, n)`
becomes `this%gwfhead(this%node_map(n))`, etc. This covers
`tsp-fmi.f90:239,440-488` and `prt-fmi.f90:118-142,254`.

### 3. Build a connection map for `gwfflowja`

Harder, because `flowja` is indexed by *position in the ja array*, not by
node, and a child model with fewer active cells will have a differently
sized/ordered `ja`. Also computed once, at `exg_ar` time:

```fortran
do n = 1, childmodel%dis%nodes
  gwf_n = node_map(n)
  do ipos = childmodel%dis%con%ia(n) + 1, childmodel%dis%con%ia(n + 1) - 1
    m = childmodel%dis%con%ja(ipos)
    gwf_m = node_map(m)
    ! find position of (gwf_n, gwf_m) in gwfmodel%dis%con%ia/ja
    flowja_map(ipos) = <search gwfmodel%dis%con%ja(gwfmodel%dis%con%ia(gwf_n)+1 : gwfmodel%dis%con%ia(gwf_n+1)-1) for gwf_m>
  end do
end do
```

Each inner search is over a short adjacency list (grid connectivity, not
graph-wide), so this is cheap and one-time. Runtime hot loops then do
`flownm = this%gwfflowja(this%flowja_map(ipos))` instead of
`this%gwfflowja(ipos)` — no change to per-timestep cost.

### 4. Everything downstream is unaffected

A child-model cell with `idomain = 0` already gets no packages, no matrix
row, no budget/output entry — existing machinery doesn't care *why* a cell
is inactive. So once (2) and (3) are in place, no further changes should be
needed in the transport/energy/particle solve itself. This is also why the
feature genuinely saves runtime (unlike zeroing parameters): excluded cells
simply don't participate in the matrix or particle tracking at all.

## Staging plan

1. **PRT first.** No matrix solve — the mapping only needs to feed
   velocity/flow lookups for particle movement. Lowest-risk place to land
   and validate the node-map / connection-map infrastructure.
2. **GWT next.** Needs the connection map for the dry/rewet logic in
   `tsp-fmi.f90:462-494` (weighted concentration from surrounding flows) and
   wherever else `flowja` feeds advection/dispersion.
3. **GWE last.** `aprovost-usgs` flagged that GWE keeps some dry cells
   "active" for conduction, unlike GWF. That's an orthogonal dry-cell
   concern, not an IDOMAIN-subset concern, but needs confirmation that it
   doesn't interact badly with a shrunk GWE domain before shipping.

## Testing

- Unit-level: node map and flowja map construction against small synthetic
  DIS/DISV/DISU pairs with hand-checkable idomain subsets.
- Regression: a GWF domain with a GWT (or PRT) domain that excludes a
  known sub-region (e.g., a corner), compared against a reference run where
  that corner's transport parameters/sources are simply zeroed out (today's
  workaround). Mass balance and concentrations elsewhere in the domain
  should match; the excluded region should be absent from GWT output
  entirely rather than present-with-zero-values.
- Confirm runtime improves for a case with a large excluded region, as
  motivation for doing this over the zero-out workaround.

## Prototype status (2026-08-18)

A first working prototype of the PRT case is implemented on this branch:

- `src/Exchange/exg-gwfprt.f90`: `exg_ar` now requires only that GWF and
  PRT share `nodesuser` (same grid shape), then verifies the subset
  property (no PRT-active cell may be GWF-inactive) via a new
  `check_and_map_domains` subroutine. That subroutine also builds, only
  when the two domains actually differ:
  - `gwf2loc` (GWF reduced node → PRT reduced node, 0 if inactive in PRT)
  - `loc2gwf` (PRT reduced node → GWF reduced node)
  - `loc2gwfja` (PRT reduced connection position → GWF reduced connection
    position), built via a neighbor search over each PRT connection's
    corresponding GWF adjacency list
  All three are built purely from `get_nodenumber`/`get_nodeuser`, which
  turned out to already be generic across DIS/DISV/DISU on the base
  `DisBaseType` — no `select type` needed for the mapping itself (only
  the `nodesuser` check needed touching).
- `src/Model/ParticleTracking/prt-fmi.f90` (`PrtFmiType`): gained the three
  map arrays (public, null unless built) and now routes every per-node GWF
  array access (`gwfhead`, `gwfsat`) and the storage-flow accumulation
  (`gwfstrgss`/`gwfstrgsy`) through `loc2gwf` when it's associated, and
  translates GWF boundary-package `nodelist` entries through `gwf2loc`
  before using them to index PRT's own per-node arrays.
- `src/Solution/ParticleTracker/Method/MethodDis.f90` and `MethodDisv.f90`:
  the `gwfflowja` face-flow lookup now goes through `loc2gwfja` when set.
- When the two domains are identical (the common case today), none of the
  above maps are built and every access site takes the original direct-index
  branch — zero behavior change for existing simulations.

**Validation** (`design/prototype/smoke_test.py` — an ad hoc smoke test,
not yet a proper autotest): three runs on the 10x10x1
`FlopyReadmeCase` grid (existing PRT test fixture) sharing one GWF flow
field:
1. PRT idomain == GWF idomain (all active): runs to normal termination
   (exercises the unmodified path).
2. PRT idomain excludes one GWF-active cell mid-grid (subset case): runs to
   normal termination. Of the 9 released particles, the 6 whose paths never
   approach the excluded cell match case 1 bit-for-bit (validates the
   node/connection maps against real, non-trivial index shifts — removing
   one cell from the middle of the reduction order changes reduced node
   numbers for most of the grid). The 3 particles whose paths do approach
   the excluded cell terminate there with `istatus = TERM_NO_EXITS`, i.e.
   PRT correctly treats the excluded-but-GWF-active cell as having no exit
   face — exactly the ICBUND-like behavior discussion #2420 asked for.
3. PRT active in a cell GWF marks inactive (violation of the subset rule):
   rejected at `exg_ar` with the new, more specific error message.

**Existing test to reconcile next:** `autotest/test_prt_exg.py` already has
an `"idmu"` case built for exactly scenario 2 above (GWF fully active, PRT
excludes one cell) and currently marks it `xfail`. That should flip to a
real assertion once this lands; its `"idmn"` case (PRT active where GWF is
inactive) should remain a rejection, now via the new error message. Running
that test requires the full `TestFramework`/MODPATH 7 comparison
scaffolding (`bin/mf6`, `bin/downloaded/mp7`, etc.) which wasn't set up in
this pass — the ad hoc smoke test above substitutes for now.

**Not yet done:**
- GWE still hard-requires identical IDOMAIN (unchanged); see GWT below for
  why GWE needs its own careful look (dry-cell handling) before attempting
  the same approach.
- `idomain < 0` (vertical pass-through) interaction with the subset check
  is untested.
- The `loc2gwfja` build does an O(connections-per-cell) neighbor search per
  connection; fine for a prototype, could be tightened later if profiling
  ever shows it matters (it's a one-time `exg_ar` cost, not per-timestep).
- No unit tests for `check_and_map_domains` in isolation.

## Prototype status: GWT (2026-08-18)

Also implemented on this branch, reusing the PRT design's node-map/
connection-map machinery, adapted for GWT:

- `src/Exchange/exg-gwfgwt.f90`: same relaxed check + `check_and_map_domains`
  as PRT (`gwf2loc`/`loc2gwf`/`loc2gwfja` on `gwtmodel%fmi`), plus a hard
  rejection if GWF has BUY or VSC active alongside a subset GWT domain
  (those packages read GWT's concentration/temperature using GWF's own
  node numbering, in `gwf-buy.f90`/`gwf-vsc.f90`, which isn't mapped —
  scoped out, see below).

**Different mechanism than PRT, by necessity.** GWT has ~25 call sites
across `gwt-mst.f90`, `gwt-dsp.f90`, `gwt-ist.f90`, `gwt-src.f90`, and
`tsp-adv.f90` that read `this%fmi%gwfsat(n)` / `gwfhead(n)` / `gwfspdis(:,n)`
/ `gwfstrgss(n)` / `gwfstrgsy(n)` directly (PRT only had a handful). Patching
every call site the way PRT's prototype does would be large and easy to
under-cover. Instead, `TspFmiType` (`src/Model/TransportModel/tsp-fmi.f90`)
gained `gwfhead_raw`/`gwfsat_raw`/`gwfspdis_raw`/`gwfflowja_raw`/
`gwfstrgss_raw`/`gwfstrgsy_raw` (aliased into GWF's arrays exactly as
before) plus a new `translate_gwf_arrays` method that runs once per
`fmi_ad` and repopulates the *public* `gwfhead`/`gwfsat`/`gwfspdis`/
`gwfflowja`/`gwfstrgss`/`gwfstrgsy` -- now GWT-owned, GWT-sized arrays --
via `loc2gwf`/`loc2gwfja`. Every existing consumer keeps reading those same
field names, unmodified, and gets correctly-translated values. This only
activates when the domains differ (`associated(this%loc2gwf)`); identical
domains still take the original direct-alias path with zero added cost.

This "translate once" approach does **not** work for `gwfpackages(ip)%
nodelist` (used ~8 places across `tsp-ssm.f90` and one place in
`tsp-mvt.f90` to look up which GWT/GWF cell a boundary-package entry
applies to): that array is `mem_reassignptr`'d directly into the live GWF
boundary package's own memory, so translating it in place would corrupt
GWF's own copy. Those ~9 call sites were each patched individually instead,
translating the fetched value through `gwf2loc` right after reading it and
before the pre-existing `if (n <= 0) cycle` guards -- which already existed
in every case, so 0 → "not part of GWT's domain" falls out for free.

**A real numerical caveat was found, not just a scoping gap.** Ad hoc
testing (`design/prototype/smoke_test_gwt.py`, a 10x10x1 grid with a
diagonal CHD-driven flow field, mirroring the PRT smoke test) showed:
baseline (identical domains) and subset (GWT excludes one GWF-active cell)
both run to normal termination and report ~0% global mass-balance
discrepancy, but the subset run's concentration at the cell bordering the
exclusion **overshoots to ~7.3** against a source concentration of 1.0 --
clearly unphysical for a conservative tracer with no dispersion.

Root cause: unlike PRT (a particle simply terminates with `TERM_NO_EXITS`
when its cell has no exit face -- an outcome PRT's algorithm already
handles correctly, since it doesn't require water-balance closure), GWT's
finite-volume advection scheme assumes each active cell's water balance
closes over *its own* set of connections. When a GWT-active cell borders a
cell that's active in GWF but excluded from GWT, the real, generally
nonzero GWF flux across that face is currently just dropped (GWT's own
`dis%con` simply has no connection there, by construction, since it's
built from GWT's own IDOMAIN) rather than folded in as a boundary term.
That breaks the per-cell balance GWT's scheme depends on and produces
local artifacts. The aggregate mass balance still looked fine in this test
because the discrepancy is a *local misallocation*, not a global leak.

This is precisely the kind of model-specific "ramification" `aprovost-usgs`
flagged in issue #1129 ("PRT is a simply different kind of animal that
isn't solving a matrix equation") -- now concretely demonstrated rather
than speculative.

**What was done about it, given the scope of a real fix:** rather than
attempt the proper fix (below) in this pass, `exg-gwfgwt.f90` now counts
these "dropped boundary connections" (GWT-active cell ↔ GWF-active/
GWT-excluded cell pairs) after building the maps, and if any exist, writes
a `WARNING` to the listing file naming the exchange and the count, so a
user isn't silently misled. It's a warning, not a hard error: a user
excluding a genuinely low/no-flow region (the common ICBUND use case from
discussion #2420) may see no practical effect, and there's no way to know
the actual flux magnitude at `exg_ar` time (before the first solve).

**What the real fix looks like:** for each GWT-active cell with one or more
dropped connections, compute the net GWF flux across them at each `fmi_ad`
(available via `gwfflowja_raw`, already resolved by the map-building code)
and inject it into GWT's transport equation as an explicit boundary term --
conceptually a weak sink/source at that cell, using the cell's own
concentration on outflow (the same `omega` convention already used
throughout `tsp-ssm.f90`'s `ssm_term`). That likely means: a new per-cell
"dropped outflow" array built alongside `loc2gwfja`, consumed either in
`tsp-adv.f90` alongside the existing face-flow terms, or as a synthetic
SSM-like term. This needs real design attention (in particular, getting
the sign/direction and inflow-concentration convention right) before it's
attempted -- it's the actual blocking item for taking this past prototype
stage for GWT.

**Not yet done:**
- The boundary-flux fix described above.
- GWE (would need the same `TspFmiType` mechanism since it shares
  `tsp-fmi.f90`, but its own exchange, `exg-gwfgwe.f90`, is untouched; it
  would also need to resolve the boundary-flux issue, plus the dry-cell
  wrinkle `aprovost-usgs` flagged).
- `gwfconn2gwtconn` (linking to GWT-GWT/parallel-decomposition connections)
  is untouched; like PRT's analogous TODO, multi-model domain decomposition
  combined with a subset domain is unexplored.
- No pytest-based regression test for GWT (only the ad hoc smoke test).

## Open questions for maintainers

- Should `idomain < 0` (vertical pass-through) cells be treated as "active"
  or "inactive" for the subset check? Current `DiscretizationBase.f90:59`
  comment treats `-1` (pass-through) distinctly from `0` (inactive) in the
  reduced-numbering scheme; need to decide how a child model's `-1` vs `0`
  interacts with a parent's `-1`/`0`/`1`.
- Any interaction with the existing GWT-GWT / GWE-GWE multi-model coupling
  path (`exg-gwfgwt.f90:370-408`) when the "parent" GWF-GWF exchange's two
  GWF models don't have matching idomain to each other? (Believed out of
  scope — that's a separate GWF-GWF constraint — but worth confirming.)
- Error UX: should mismatch (child active where parent inactive) be a hard
  error immediately, or should MF6 auto-correct by forcing those cells
  inactive with a warning? Proposal above assumes hard error for now,
  matching current behavior's strictness.
- The big one, from the GWT prototype: should MF6 *block* a subset GWT
  domain until the boundary-flux fix exists (i.e. today's warning becomes a
  hard error), or is warn-and-proceed acceptable, given the primary use
  case (ICBUND-style exclusion of a genuinely low/no-flow region) may see
  no practical artifact? Needs a real answer — and probably a real test
  with a non-trivial cross-boundary flux, not just the ad hoc smoke test
  here — before this is proposed upstream.
