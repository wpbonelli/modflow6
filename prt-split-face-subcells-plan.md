# PRT tracking through split faces: subcell splits for DISU and nested grids

Plan for tracking particles through cells whose faces are shared with more than one
neighbor. This arises on DISU grids and at non-matching PRT-PRT exchange boundaries
(nested grids). It expands the "Tracking on DISU grids" item of the PRT development plan
(claude.ai doc "MODFLOW 6 PRT development plan"), which left open how hard internal
regridding would be.

## Terms

- Split face: a cell face shared with two or more neighbors. A lateral face is split when
  neighbors cover different elevation ranges (DISU `IHC=2`, vertically staggered). A top
  or bottom face is split when several cells lie above or below it.
- Sub-face: the part of a split face shared with one neighbor. Each sub-face carries that
  connection's flow from GWF.
- Subcell split: a partition of one cell into subcells such that every sub-face belongs to
  exactly one subcell, with flows assigned to the faces between subcells. Built per cell
  when a particle enters it, not as a refined grid for the whole model.

## Problem

PRT has tracking methods only for DIS and DISV. `prt.f90` selects one by discretization
type (`src/Model/ParticleTracking/prt.f90:313`) and DISU falls through. The cell
definition holds one neighbor and one flow per face, including top and bottom
(`src/Solution/ParticleTracker/Domain/CellDefn.f90`, `facenbr`, `faceflow`). A split face
cannot be represented.

The simplest DISU method sums the sub-face flows onto the face, applies the existing cell
method with uniform flux on that face, and lets the exit point choose the neighbor. That
is wrong in two ways where neighbors' flows differ:

- Particles cross the face in proportion to each sub-face's area, not its flow. Advective
  mass is misallocated between neighbors.
- Sub-face flows of opposite sign cancel. Water entering through one sub-face and leaving
  through another on the same face is invisible, and particles can neither enter nor
  leave through those sub-faces as the flow solution says.

PRT-PRT exchanges on nested grids (draft #2751, branch `prt-exg`) have the same problem:
a coarse cell at the exchange boundary has a face shared with several fine cells in the
other model. DISV models coupled this way need it even without DISU.

DISU cells are vertical prisms (`TOP`, `BOT`, `AREA`, optional `VERTICES`/`CELL2D`), so
subcells are prisms too. A fully 3D polyhedral discretization (issue #3207,
`full3d-disu-feasibility.md` on branch `full3d-disu-feasibility`) would need a different
tracking method and is not addressed here.

## Precedent in the code

DISV quad refinement already does a subcell split. A rectangular cell with one or more
lateral faces split in two (`can_be_quad`) is divided into four rectangular subcells
(`src/Solution/ParticleTracker/Domain/CellUtil.f90:114`, `cell_poly_to_quad`). External
subcell flows come from the half-faces. Top, bottom and distributed flows are spread
evenly, a quarter to each subcell. The four internal flows come from mass balance plus a
fixed closure. `MethodCellPollockQuad` tracks Pollock within each subcell and passes the
particle between subcells inside the cell. The general split described below has this as
a special case.

## Why per-cell splits suffice

The alternative is to propagate refinement into a conforming grid, where every face is
shared by exactly one cell on each side. Converting DISU to DISV requires this: plan-view
refinement propagates vertically, because all layers share one mesh, and vertical
refinement propagates laterally, because lateral connections exist only within a layer.
A per-cell split stops at the cell boundary instead. The neighbor on an unsplit face sees
the split cell as one cell, and that face's flow is shared among the subcells touching it
by area. There is no correctness reason to prefer propagation when the goal is tracking
the GWF solution:

- GWF is the source of truth. It provides one flow per connection and storage and source
  terms per cell. Any velocity field inside a cell is a reconstruction, and a valid one
  reproduces every connection flow, balances every cell, and keeps normal flux continuous
  across faces. PRT should track the same flows GWT and GWE use.
- Propagation adds no information. The extra faces it creates carry flows GWF never
  computed. They come either from sharing by area, the per-cell assumption applied in
  more places, or from re-solving flow on the refined grid, which gives a different model
  that no longer reproduces the connection flows or matches GWT's budget. A re-solve
  constrained to the GWF connection flows couples cells only through how each face's flow
  is shared among its parts. Fixing that sharing decouples it into per-cell problems.
- Uniform flux per connection is already the model's assumption. GWF's two-point flux
  treats a connection as one flow over one area, and the existing tracking methods assume
  uniform flux on each face. The per-cell split keeps this on every unsplit face and every
  sub-face. Both sides of a face then see the same uniform flux, so normal flux is
  continuous even where the two cells are cut at different elevations. Propagation would
  make flux non-uniform at a resolution GWF never solved.

What per-cell splits do change is the path within a split cell, through the cut positions
and the internal flows. That is a deviation bounded by the cell, of the same order as the
within-cell interpolation the existing methods make everywhere. The neighbor a particle
exits to follows the sub-face flows either way.

Propagation would matter only to resolve flow more finely than GWF did. That is
downscaling, a different goal, and belongs outside PRT.

## Approach

### Face partition (static)

At setup, for each cell, list each face's sub-faces: neighbor connection index and
extent.

- Lateral faces: extent is the elevation interval the two cells share,
  `[max(bot_n, bot_m), min(top_n, top_m)]`. This is the interval NPF uses for a staggered
  connection (`src/Model/GroundWaterFlow/gwf-npf.f90:855`), so PRT's geometry matches
  GWF's conductance.
- Top and bottom faces: extent is the intersection of the two cells' footprints in plan.
  This needs polygon clipping on `CELL2D`/`VERTICES`. For rectangular cells with
  rectangular neighbors (quadtree), intersection reduces to interval overlap in x and y.
- A DISU package without vertices is an input error for PRT.
- Storage is one record per connection, of the same order as the existing connection
  arrays.

Boundary flows assigned to a split face (FMI `BoundaryFlows` by `IFLOWFACE`) are spread
over its sub-faces by area. FMI sizes `BoundaryFlows` by the maximum polygon vertex count
of a DISV grid (`src/Model/ParticleTracking/prt-fmi.f90:164`). DISU needs the same from
its own vertex data.

### Subcell split (dynamic)

Built when a particle enters a cell with a split face, as the DISV method already
rebuilds a cell definition on every entry. Cells without split faces use the existing
methods unchanged. "Dynamic" here means two things:

- Flows change every time step.
- In water-table cells the saturated thickness changes, so lateral sub-face extents are
  clipped at the water table and slabs above it vanish.

There are two cases, with different numerics.

Vertical stacking (lateral faces split only). Cut the cell at every neighbor's top and
bottom that falls inside it, giving n stacked slabs. Each slab gets its lateral sub-face
flows, a share of the distributed flow by saturated volume, and the cell's top flow
(top slab) or bottom flow (bottom slab). The n-1 flows between slabs are the unknowns.
There are n slab balance equations, one redundant given the cell balance, so the
internal flows follow uniquely by summing from the bottom. No closure assumption is
needed. Each slab is a prism with one neighbor per face, so Pollock (rectangular) or
ternary (polygonal) applies unchanged. The one assumption: where a neighbor spans more
than one slab (because another face's neighbors made the cuts), its flow is shared by
overlap thickness.

Plan-view splits (top or bottom split). Partition the footprint by the neighbors'
footprints, giving columns of subcells, then stack slabs within columns as above if
lateral faces are also split. The internal flows are underdetermined: internal faces
outnumber the independent balance equations (for the 2x2 quad case, four internal faces
and three independent equations, one free value). A closure is needed. Options:

- Local finite-volume solve on the subcell graph: two-point conductances from subcell
  geometry, Neumann conditions from sub-face flows, sources by volume. The heads are
  unique up to a constant and the flows are unique. The system is small (tens of
  unknowns) and dense LU suffices. PRT receives no hydraulic conductivity through FMI
  (budget files carry none), so this uses isotropic K within the cell unless K is passed
  in. Within one cell only the anisotropy ratio matters. External sub-face flows must be
  imposed as uniform flux per sub-face, so normal flux stays continuous with the
  neighbor.
- An algebraic closure in the style of the quad method. Cheap, but has no physical basis
  beyond symmetry and does not generalize past regular splits.

Subcells from clipping convex footprints are convex. A non-convex parent cell can yield
non-convex subcells, which ternary does not handle without further splitting.

### Tracking in a split cell

- The particle enters the subcell that owns its entry sub-face.
- It is tracked in that subcell with the existing subcell methods (Pollock on rectangular
  subcells, ternary on polygonal ones) and passed between subcells inside the cell, as
  `MethodCellPollockQuad` already does.
- On exit through a sub-face, the neighbor is the sub-face's connection. No exit-point
  rule is needed.
- Elevation maps from the source sub-face's extent onto the same extent in the neighbor,
  which is the overlap. On split faces this is the "overlap-relative" mapping from the
  "Vertically offset cells" item, obtained as a consequence rather than an option.

### Data structures

The cell definition needs a variable number of sub-faces per face: for each face, a
list of (connection, flow, extent). DIS and DISV cells have one sub-face per face, so
their behavior and results are unchanged. The subcell split itself is scratch storage on
the method object, sized to the largest split seen.

## Costs

Not measured; estimates from the code structure.

- Memory: the static face partition, about one record per connection, plus per-method
  scratch for one cell's subcells (per thread once tracking is threaded). With no
  refined grid for the whole model there is no per-cell cache. Caching splits per cell
  per time step is possible but optional.
- Runtime: only cells with split faces pay. Vertical stacking is a short sweep. Plan-view
  splits cost a small dense solve per cell entry. Clipping is static and done at setup.
- Maintenance: vertical stacking is small and determinate, and a subcell mass-balance
  check tests it directly. Plan-view splits carry the lasting cost: robust clipping,
  tolerances for degenerate overlaps, and non-convex cells. Expect something comparable
  to the ternary method (`MethodCellTernary`, `MethodSubcellTernary`,
  `TernarySolveUtils`). Against that, the general split can replace `CellRectQuad` and
  `MethodCellPollockQuad`.

## Rejected alternatives

Refined grid for the whole model, built at setup. A refined DISU could be built in code
the way interface models build theirs. But tracking on it needs heads or flows on the
refined cells, which means solving a flow problem on the refined grid: a second flow
model. Per-cell splits redistribute GWF's flows locally and need nothing more.

Reusing interface models. `GwfInterfaceModel` (`src/Model/Connection/`) builds a DISU from
a stencil of existing cells in two models (`GridConnection`), gives it an NPF instance,
and adds its terms to the global matrix during the flow solve. It never creates new cells
and never solves a local problem. A model object with NPF per cell entry would be far too
heavy (memory manager allocations, model setup). What carries over: `GridConnection`
for PRT-PRT exchanges in parallel runs, where the neighbor model's cell geometry is remote.

Preprocessing with imod-python. Its regridding
(https://deltares.github.io/imod-python/user-guide/07-regridding.html):

- regrids model input packages only, then the user re-solves, producing a different flow
  model rather than tracking through the original one;
- supports DIS and DISV, and lists DISU as unsupported;
- regrids only in the x-y plane, not across layers;
- does not claim conservation of face flows.

Converting a vertically refined DISU to DISV means refining every column, which removes
the reason for DISU, and the re-solved model is no longer the calibrated one.

Offline flow downscaling. The same subcell split could run in Python, writing a refined
grid and budget file for PRT to read through FMI. It loses the live GWF-PRT exchange,
inflates output, and still needs a DISV target. It is useful as a prototype (next steps),
not as the delivery.

## Validation

- Subcell mass balance: internal flows close each subcell's balance to round-off.
- Split invariance: a cell split artificially (neighbors at equal elevations, uniform
  flow) gives the same paths as the unsplit cell.
- DIS, DISV and DISU versions of one grid give the same paths (the shared invariant in
  the PRT plan).
- A quad-refined grid as DISV (existing quad method) and as DISU (general split): the
  local-solve closure should match the quad closure on regular splits, or the difference
  should be understood and documented.
- Sub-face flux weighting: many particles released uniformly upstream of a split face
  cross each sub-face in proportion to its flow. This is the metric the summed-face
  method fails.
- Comparison with mod-PATH3DU on nested and vertically refined DISU grids.
- PRT-PRT on nested grids: a nested pair of models gives the same paths as the equivalent
  single DISU model.

## Open questions

- Plan-view closure: local solve, or leave plan-view splits to the "Velocity
  reconstruction on general grids" item (a Waterloo-style fit per cell, as in
  mod-PATH3DU, takes sub-face flows directly) and implement only vertical stacking?
- Is isotropic K in the local solve acceptable, or should FMI pass K (live exchange only,
  since budget files have none)?
- Non-convex parent cells: reject on PRT DISU input, or split further?
- Should the split depend on saturation only through clipping, or also re-cut slabs at
  the water table?

## Next steps

1. Prototype the face partition and subcell flows in Python against a DISU binary grid
   file and budget file (FloPy `MfGrdFile`, `CellBudgetFile`). Measure the summed-face
   method's sub-face flux-weighting error on vertically refined and nested test grids.
   This decides whether vertical stacking and plan-view splits are each worth building.
2. DISU method with summed faces, as in the PRT plan's "Tracking on DISU grids" item,
   including the sub-face lists in the cell definition. This is the baseline and the
   fallback.
3. Vertical stacking in Fortran, with the subcell mass-balance unit test.
4. Plan-view splits by local solve, or defer to velocity reconstruction (open question).
5. Use the same face partition for non-matching faces in the PRT-PRT exchange (#2751),
   with neighbor geometry from the other model, through `GridConnection` in parallel
   runs.
