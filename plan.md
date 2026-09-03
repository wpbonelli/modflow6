# Plan: correct the PRT IDOMAIN-subset implementation

- Status: agreed, not yet implemented
- Branch: `prt-idomain-subset`
- Related: issue [#2962](https://github.com/MODFLOW-ORG/modflow6/issues/2962),
  discussion [#2420](https://github.com/MODFLOW-ORG/modflow6/discussions/2420),
  `design/idomain-subset-coupling.md`
- Context: the branch already relaxes the GWF-PRT exchange check from "IDOMAIN
  identical" to "PRT IDOMAIN is a subset of GWF IDOMAIN", and builds
  `noder_gwf2prt` / `noder_prt2gwf` / `ipos_prt2gwf` maps in
  `check_and_map_domains`. Two defects remain (below), plus a weak test.

## Defect 1 — inconsistent numbering for the wired GWF arrays

`exg_ar` sets `prtmodel%fmi%gwfhead => gwfmodel%x`, `gwfsat => npf%sat`,
`gwfspdis => npf%spdis`, `gwfceltyp => npf%icelltype`,
`gwfstrgss/gwfstrgsy => sto%strg*`. These are GWF **reduced** numbering.
The branch translates the index at *some* read sites (`prtfmi_ad` dry/rewet
loop, `accumulate_flows` storage gather, `get_gwfflowja`) but not others:

- `src/Solution/ParticleTracker/Method/MethodDis.f90:251,367,370`
- `src/Solution/ParticleTracker/Method/MethodDisv.f90:323,370`
- `src/Solution/ParticleTracker/Method/MethodModel.f90:107`
- `src/Model/ParticleTracking/prt-prp.f90:704,709`

all read `this%fmi%gwfsat(ic)` with `ic` a **PRT** reduced node. When PRT
excludes any cell that is not the last reduced node, PRT and GWF numbering
diverge and these read the wrong cell's saturation (wrong cell top, wrong
z-remap between cells, wrong PRP release-cell weighting). The shipped test
only excludes the last node, so nodes `1..N-1` keep identical numbers and
the bug is masked.

### Fix (decided): FMI owns PRT-sized copies

Matches the "Array transfer" bullet of #2962: in mapped mode FMI owns
`GWFHEAD/GWFSAT/GWFSPDIS/GWFFLOWJA/GWFSTRGSS/GWFSTRGSY` sized on the PRT
grid and refills them from the GWF model each step. Every consumer then
indexes with native PRT node / `ipos` numbers and needs no per-site
translation. This is also what base `FlowModelInterface` already does in
budget-file mode (`allocate_arrays`, `flows_from_file` branch).

**`src/Model/ParticleTracking/prt-fmi.f90`**

1. Add source pointers to `PrtFmiType` (borrowed, not owned):
   `gwfhead_src`, `gwfsat_src`, `gwfspdis_src`, `gwfceltyp_src`,
   `gwfstrgss_src`, `gwfstrgsy_src`. Nullify in `fmi_cr`; nullify (do not
   deallocate) in `prtfmi_da`.
2. New `prtfmi_ad` step, run **before** `accumulate_flows` and the
   dry/rewet loop, guarded by `associated(this%noder_prt2gwf)`:
   ```fortran
   do n = 1, this%dis%nodes
     ng = this%noder_prt2gwf(n)
     this%gwfhead(n)     = this%gwfhead_src(ng)
     this%gwfsat(n)      = this%gwfsat_src(ng)
     this%gwfceltyp(n)   = this%gwfceltyp_src(ng)
     this%gwfspdis(:, n) = this%gwfspdis_src(:, ng)
   end do
   do ipos = 1, this%dis%con%nja
     this%gwfflowja(ipos) = this%gwfflowja_src(this%ipos_prt2gwf(ipos))
   end do
   if (this%igwfstrgss /= 0) then
     do n = 1, this%dis%nodes
       this%gwfstrgss(n) = this%gwfstrgss_src(this%noder_prt2gwf(n))
     end do
   end if
   ! ...gwfstrgsy likewise, guarded by igwfstrgsy
   ```
   Refill `gwfspdis`/`gwfceltyp` too even though PRT tracking does not read
   them today — uniform treatment, no latent trap if a consumer is added.
3. Revert the `ng` branching in the dry/rewet loop to plain
   `this%gwfsat(n)` / `this%gwfhead(n)`.
4. Revert the `noder_prt2gwf` branch in `accumulate_flows` storage
   accumulation to the simple `this%StorageFlows + this%gwfstrgss` form.
5. Delete `get_gwfflowja`; it is no longer needed.
6. `prtfmi_da`: when `associated(this%noder_prt2gwf)`, explicitly
   `mem_deallocate` the arrays that base `fmi_da` only frees under
   `flows_from_file` — `gwfstrgss`, `gwfstrgsy`, `gwfceltyp`. Base `fmi_da`
   already frees `gwfhead/gwfsat/gwfspdis/gwfflowja` via its "could be from
   mem_checkin" path, which also covers an owned allocation.

**`src/Exchange/exg-gwfprt.f90` (`exg_ar`)**

- Identical-domain case (maps not built): unchanged — `=>` assignment plus
  `mem_checkin` as today.
- Subset case (`check_and_map_domains` built the maps): instead of
  `=>`/`mem_checkin`, `mem_allocate` each array on
  `prtmodel%fmi%memoryPath` sized on the PRT grid
  (`gwfflowja(nja)`, `gwfsat(nodes)`, `gwfhead(nodes)`,
  `gwfspdis(3,nodes)`, `gwfceltyp(nodes)`, and `gwfstrgss/gwfstrgsy(nodes)`
  under the existing `inmst`/`insto`/`iusesy` guards), and set the matching
  `*_src` pointer to the GWF array. Set `igwfstrgss`/`igwfstrgsy`/
  `igwfceltyp` exactly as now.
- Have `check_and_map_domains` return `differs` (or expose whether maps
  were built) so `exg_ar` can branch. Simplest: make it a `logical`
  function, or set a flag on the exchange.

**`src/Solution/ParticleTracker/Method/MethodDis.f90`,
`MethodDisv.f90`**

- Revert `q = this%fmi%get_gwfflowja(ipos)` to
  `q = this%fmi%gwfflowja(ipos)` and drop the added `ipos` local where it
  was only introduced for that call.
- `MethodModel.f90`, `prt-prp.f90`: no change needed — the untranslated
  `gwfsat(ic)` reads become correct once `gwfsat` is PRT-sized.

## Defect 2 — flow across the cut boundary is dropped

A PRT boundary cell `C` active in both models may have a GWF neighbor `E`
that is active in GWF but excluded from PRT. GWF moves water across the
`C`-`E` face; PRT has no connection there, so `facenbr` for that face is 0,
`load_cell_face_flows` contributes nothing, and the flow vanishes. `C` is
then not mass-balanced: its interior velocity field is distorted and a
particle may hit a spurious `TERM_NO_EXITS` in `C` instead of being
advected to the cut and terminating `TERM_BOUNDARY` there.

`gwfflowja` is "positive into a cell" (`gwf.f90:755`); `faceflow > 0` and
`distflow` ("net distributed flow into cell", `CellDefn.f90:45`) share that
sign convention.

### Fix (decided): Option B — diagnostic + fold into `distflow`

**Build the dropped-connection list** (in `check_and_map_domains`, after
the maps, subset already guaranteed):

- For each PRT node `pn`, `gn = noder_prt2gwf(pn)`; walk GWF connections
  `jpos` of `gn`; any GWF neighbor `gm` with `noder_gwf2prt(gm) == 0` is a
  dropped connection. Store as a PRT-node-keyed CSR structure on the FMI
  object: `iadrop(pn+1)` offsets into `jadrop(:)`, where `jadrop` holds GWF
  `ja` positions (`jpos`). `mem_allocate` on the FMI memoryPath;
  `mem_deallocate` in `prtfmi_da` under `associated` guard.

**Per-step accumulation** (in the new `prtfmi_ad` refill block):

```fortran
do pn = 1, this%dis%nodes
  this%dropflow(pn) = DZERO
  do k = this%iadrop(pn), this%iadrop(pn + 1) - 1
    this%dropflow(pn) = this%dropflow(pn) + this%gwfflowja_src(this%jadrop(k))
  end do
end do
```

`dropflow` is a PRT-sized FMI array, positive = net flow into the PRT cell
from excluded cells.

**Diagnostic**: after accumulation, for any `pn` where
`abs(dropflow(pn))` exceeds `frac * sum(abs(kept face flows of pn))`
(propose `frac = 1.0d-2`), emit a **one-time** warning (guard with a
logical so it fires once per run) naming the cell and the percentage, e.g.
`"PRT cell <nodestr> on the active-domain boundary of exchange <name>
carries <pct>% of its flow across the boundary into cells excluded from
the PRT domain; place the boundary where this flow is negligible or
interpret nearby results with care."` The "kept face flow sum" can be
taken from `gwfflowja_src` over `gn`'s non-dropped connections, computed in
the same walk.

**Fold into the tracking field** (`MethodDis.f90` and `MethodDisv.f90`,
`load_cell_flows`):

```fortran
defn%distflow = this%fmi%SourceFlows(defn%icell) + &
                this%fmi%SinkFlows(defn%icell) + &
                this%fmi%StorageFlows(defn%icell)
if (associated(this%fmi%dropflow)) &
  defn%distflow = defn%distflow + this%fmi%dropflow(defn%icell)
```

This makes the boundary cell mass-consistent: the cross-cut flow is
represented as a uniformly distributed source/sink over the cell volume
(the PRT analog of GWT's `FLOW_IMBALANCE_CORRECTION`). Approximate — a
particle that should reach the cut face may instead be drawn toward a
distributed sink and end in the cell interior — but in the regime #2962
requires (flux across the cut negligible) the effect is small, and it
removes the spurious hard failures when the domain is sized imperfectly.
Note `distflow` does not set `iweaksink`, so such a termination is not
labelled a weak sink.

## Minor items

- **Map build timing.** #2962 says build maps in `exg_df`; the branch uses
  `exg_ar`. Fine for PRT — the maps size nothing, they only translate
  indices, and `dis%con` is available at `ar`. Add a one-line comment
  noting the deliberate choice.
- **Input definitions.** Update the PRT DIS and DISV `IDOMAIN`
  descriptions in `doc/mf6io` / the dfn files to state the PRT active
  domain may be a subset of the coupled GWF domain.
- **PRP release into an excluded cell.** Already rejected: a release point
  in a cell with `idomain = 0` in the PRT model fails PRT's own idomain
  check before the exchange is consulted. No new check needed; confirm
  with a note in the test.
- **Release note.** Extend the existing `develop.toml` feature entry to
  state the boundary of the reduced PRT domain behaves as an effective
  no-flow boundary and must be placed where cross-boundary flow is
  negligible (mirror jdhughes-dev's GWT sizing guidance).

## Test rework (`autotest/test_prt_exg.py`)

Current `prtexg01idmu` excludes only `idomain[-1,-1,-1]` (the last node),
so PRT and GWF numbering stay identical for every other cell and neither
defect is exercised.

- **`prtexg01idmu`**: exclude an **interior block** away from the last
  node — e.g. a 2x2 patch of cells in the middle of the 10x10 grid, or a
  rectangular sub-region the release-point pathlines from the top-left
  cell pass near but not through. This forces PRT reduced numbering to be
  offset from GWF's for all cells after the excluded block, exercising the
  node and connection maps and the refill.
  - Assert the run completes.
  - Assert the excluded cells produce **no** track output (absent, not
    present-with-zeros).
  - Compare pathlines/travel times for the released particles against the
    full-domain `prtexg01` run: they should match to tracking tolerance
    for cells the two runs share, provided the excluded block does not lie
    on a pathline. Pick the block accordingly.
- **`prtexg01idmn`**: keep as-is (PRT active where GWF inactive) — still an
  `exg_ar` hard error, still `xfail`.
- Optionally add `prtexg01idmb`: an excluded block placed so that real GWF
  flow crosses the cut, asserting the diagnostic warning appears in the
  listing file and the run still completes (exercises Defect 2 Option B).

## Why this problem exists at all

PRT holds pointers to the *complete* GWF arrays, so the flow data is
visible. What it lacks is **addressability**. The entire tracking engine is
written against `this%fmi%dis` — PRT's own discretization, whose `nodes`,
`con%ia/ja`, `con%nja` and node-user maps are all derived from PRT's own
`idomain`. A cell or connection PRT's `idomain` excludes has no reduced
index in PRT's `dis`, so no engine loop visits it and no PRT-sized array
has a slot for it. The GWF flow across a dropped connection is reachable in
principle (its GWF `ja` position is computable) but there is no PRT-side
structure that expects it — hence the dropped-connection list and the
`distflow` lumping.

Two implementation choices, not physics, cause this:

1. **`idomain` is overloaded.** One `idomain = 0` flag means three things
   at once: "produce no results here" (the only thing the user wants to
   control), "renumber the grid as if this cell were absent", and "present
   a no-flow boundary to neighbours, so their `ja` is built shorter".
   MT3DMS's `ICBUND` was a separate array carrying only the first meaning;
   the grid stayed whole. MODFLOW 6 folded all three into `idomain` plus a
   per-model `dis`.

2. **PRT depends only on its own `dis` + FMI arrays**, by design, because
   it must also run standalone from a budget/grb file with no live GWF
   model. The GWF-PRT exchange grafts live GWF pointers onto that
   self-contained structure rather than letting PRT borrow GWF's
   connectivity. When the two `idomain`s diverge the two `dis%con` graphs
   diverge, and the exchange maps can only bridge connections present on
   both sides.

The "ideal" treatment — track the particle through excluded cells on GWF's
field and merely stop *recording* it — would require PRT to iterate GWF's
connectivity (or a shared `dis`) behind a separate "track/report here?"
mask: essentially re-introducing `ICBUND` semantics. That is a much larger
change and complicates the standalone-from-file path, which has no GWF grid
to borrow.

It is also probably not what a user wants. For "I only care about this
sub-region", terminating a particle at the boundary (`TERM_BOUNDARY`,
which PRT already does) is the desired behaviour — pathlines should not
wander through excluded cells and end somewhere arbitrary. The only real
defect is that the *last step*, inside the boundary cell, is computed from
a mass-inconsistent velocity field. Option B fixes exactly that and
nothing more, which is the right scope rather than a compromise.

## Implementation order

1. Defect 1 refactor (FMI owns copies) + revert the partial per-site
   translations. Run `autotest/test_prt_exg.py` — still passes with the
   last-node exclusion.
2. Rework `prtexg01idmu` to an interior block; confirm it now fails
   (proves the refactor is actually exercised), then confirm it passes.
3. Defect 2: dropped-connection list, `dropflow`, diagnostic, `distflow`
   fold-in. Add `prtexg01idmb`.
4. Minor items: dfn/doc text, release note, comments.
5. `ruff` the test file; `fprettify` touched Fortran per repo style.
