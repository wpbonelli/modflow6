# Feasibility analysis: fully 3D polyhedral DISU (issue #3207)

Review of the `FULL_3D` DISU proposal in
[MODFLOW-ORG/modflow6#3207](https://github.com/MODFLOW-ORG/modflow6/issues/3207) and its
implementation on `stebond/modflow6@full3d-disu` (commits `239bff4709` feature,
`9c3dc60cae` CMake build), based on `develop` at `6565dc81`.

**Method: static review only.** The fork was not built and its tests were not run. Claims
below about runtime behavior are inferred from reading the source and are marked
*(unverified)* where it matters.

## TL;DR

- **The numerical idea is sound and well established.** A two-point flux finite-volume
  scheme on a polyhedral mesh whose input is "cell volume + cell center + (face area, two
  center-to-face distances) per connection" is exactly the TOUGH2/TOUGH3 `MESH` model
  (`ELEME` volumes, `CONNE` `D1/D2/AREA/BETAX`). For *confined, fully saturated* flow and
  transport on K-orthogonal (Voronoi-type) meshes there is nothing new to invent.
- **The implementation does not reformulate MODFLOW for 3D cells. It disguises them as
  prisms.** Each `FULL_3D` cell is stored as a prism with `top = 1`, `bot = 0`,
  `area = volume` (`Disu.f90:1179`). Most of the model (STO, GWT MST/DSP/ADV, GWE EST/CND,
  FMI) "works" because a thickness of 1 makes `area * (top - bot) == volume` and
  `hwva * thksat == hwva`. Correctness rests on coincidence rather than on the geometry
  being represented.
- **Your instinct is right that real support needs reformulation and refactoring,** but
  only for the parts the fork excludes. Confined flow needs no new math. Anything with a
  water table (convertible cells, Newton, UZF, CSUB, drying/rewetting, staggered
  connections), buoyancy, PRT, and faithful post-processing does. None of those can be
  added later without polyhedral geometry (vertices/faces) in the *input*, which this
  design deliberately omits.
- **As written it is not mergeable.** The capability is fenced off with ~60 scattered
  `is_full_3d()`/`C3D_FULL` checks, while ~250 references to `dis%top/bot/area` across 37
  source files remain unguarded. Several code paths are unguarded and likely wrong (parallel
  interface models, exchanges between models with different origins, GWT/GWE package sets,
  binary grid file). The tests never exercise a non-hexahedral cell. The PR also bundles
  unrelated changes (CMake build, VTU writer for a third-party format, PRT/UZF tweaks).
- **Recommendation:** interested in principle, but not as a DISU option. Ask for (1) a
  separate discretization type, (2) preceded by a small framework refactor that routes
  geometry through `DisBaseType` methods instead of raw `top/bot/area` arrays, (3) scoped to
  confined GWF first, (4) with polyhedral geometry (vertices + faces) as input or at least
  as an optional input that is used for validation, and (5) VTU output dropped in favor of
  FloPy-side export. Details and a suggested response are at the end.

---

## 1. What "fully 3D" means here (and what the issue's first paragraph says)

MODFLOW 6 "fully unstructured" DISU is unstructured in *connectivity*, not in *cell
shape*. Every cell is a vertical prism (`TOP`, `BOT`, `AREA`), and every connection is
vertical (`IHC=0`, `HWVA` = area) or horizontal (`IHC=1/2`, `HWVA` = width, later
multiplied by a saturated thickness). The issue's first paragraph says that users asking for
"fully unstructured" grids (discussion #1232) still end up with layered prisms, because
that is all the available generators (GRIDGEN, FloPy Voronoi/triangle) produce. The
proposal relaxes the prism constraint entirely: cells are arbitrary polyhedra (3D Voronoi,
Delaunay duals) with no top, no bottom, no layers and no horizontal/vertical distinction.

## 2. What the fork actually does

| Area | Change |
|---|---|
| Input (DFN) | `FULL_3D` keyword; `VOLUME`, `XCELL`, `YCELL`, `ZCELL` griddata; `AREA` reinterpreted as optional "exposed area" for RCH/EVT; `IHC=3`; `HWVA` = face area. Applied to GWF, GWT and GWE DISU. |
| Internal geometry | `top1d=1`, `bot1d=0`, `area1d=volume` placeholders; new `zc`, `volume`, `surfarea` arrays (`Disu.f90:1179-1229`). |
| Direction | `connection_normal`/`connection_vector` for `IHC=3` return the center-to-center unit vector (`Disu.f90:1672`, `:1730`). There is no face-normal input. |
| Conductance | `C = A / (CL1/K1 + CL2/K2)` (`fcond_sat`); `K` = `hyeff` of the full tensor along the center-to-center direction (`gwf-npf.f90:2677`). |
| Specific discharge | Area-weighted least-squares fit of face fluxes per cell (`gwf-npf.f90:3014`). |
| XT3D | ANGLDEGX check bypassed. `IHC=3` treated like vertical for the area term and as "not all horizontal". |
| Guards | `check_full_3d` in `gwf.f90:218` (package allow-list); NPF/STO option and ICELLTYPE/ICONVERT checks; MAW/LAK/SFR special cases; GWT-GWT and GWE-GWE exchanges rejected. |
| Exchanges | GWF-GWF `IHC=3` with `center_vector`; XT3D on exchange rejected. |
| Output | Binary grid file writes `TOP/BOT` as a cube of equal volume around `ZCELL` (`Disu.f90:1586`). New `VtuOutput.f90` reads a `.v3d` polyhedron file "written by VORO2MESH" and writes VTU/PVD per time step. |
| Unrelated | `CellUtil.f90` dummy-argument `pointer` attribute removed; a `use SparseModule` added in `gwf-uzf-input.f90`; separate CMake/`build_vs.bat` commit. |

About 2,900 added lines, almost all in one commit.

## 3. Numerical formulation

### 3.1 Confined flow: sound, with one important caveat

The TPFA conductance is the standard one and reduces to the existing vertical-connection
formula. Taking K along the center-to-center direction via `hyeff` is the same approximation
MF6 already makes for rotated anisotropy.

**Caveat: the face normal is not an input.** For `IHC=3` both TPFA and XT3D use the
center-to-center vector as the face normal. That is exact only on orthogonal meshes
(Voronoi with generating points as centers, or Delaunay duals built from circumcenters). It
is not exact for centroid-based duals, CVT meshes whose "centers" are reported as centroids
rather than generators, or any mesh where `XCELL/YCELL/ZCELL` is the centroid. Classic DISU
has `ANGLDEGX` precisely so the face normal can differ from the center line. `FULL_3D` has
no equivalent, so:

- TPFA is inconsistent on non-orthogonal meshes and with anisotropy (non-K-orthogonal),
  which is the usual case.
- XT3D, the tool MF6 offers to fix exactly that, is given the wrong normal on
  non-orthogonal meshes. It does correct anisotropy-induced non-K-orthogonality on
  orthogonal meshes, which is the case the fork tests.
- The only consistency check is a warning when `CL1+CL2` differs from the center distance
  by more than 25% (`Disu.f90`, `disu_ck`). That is loose, and it does not check
  orthogonality at all.

A face-normal input (3 components per connection), or face vertices from which it can be
computed, is needed before this can be called general polyhedral support.

### 3.2 Confined storage and transport: correct only by coincidence

STO, GWT MST, GWE EST and FMI compute `area * (top - bot)`. DSP/CND fall through to the
horizontal branch for `IHC=3` and use `thksat = (top - bot) * sat = 1` (`gwt-dsp.f90:832`).
The central-difference advection weights use `cl1/cl2` for any non-vertical connection.
All of these give the right answer only because of the placeholders. Nothing in the code
records that these call sites depend on it, so any future change to them that uses `top` or
`bot` as an elevation (rather than inside a difference) silently breaks `FULL_3D`.

### 3.3 Unconfined flow: this is the real reformulation

MF6's unconfined machinery assumes a prism:

- `sat = (h - bot) / (top - bot)` (quadratic-smoothed), linear in head;
- horizontal conductance ∝ saturated thickness × width, upstream-weighted;
- Newton derivatives of storage and conductance in terms of `top - bot`;
- convertible storage switching between `Ss·V` and `Sy·A`.

For a general polyhedron the equivalents are a per-cell *hypsometric* function
`V_n(h)` (volume below elevation `h`, piecewise cubic in `h`), its derivative (the
horizontal cross-section area at `h`, which replaces `Sy·A`), and a per-face wetted-area
function `A_f(h)` for conductance. All of these can be computed from vertices and faces,
for example by clipping each polyhedron against a horizontal plane and tabulating. But:

1. They require geometry that `FULL_3D` input does not contain.
2. Newton–Raphson in MF6 is built around the prism derivative forms in NPF, STO and the
   advanced packages. Each would need a geometric derivative from the discretization.
3. Upstream weighting assumes "horizontal" connections with a saturated-thickness
   multiplier. For arbitrary faces the wetted-area fraction is not monotone in the simple
   way the existing smoothing functions assume.

This is feasible (reservoir and vadose-zone codes do it), but it is a project in its own
right and not an incremental extension of this fork. The issue's target users (geothermal,
reservoir, fractured rock) are mostly fully saturated, so the confined-only scope is
defensible. The scope should be stated as permanent, unless geometry goes into the input
now.

### 3.4 Other physics

- **BUY:** needs elevation differences between cell centers and face midpoints along the
  connection. `ZCELL` gives the former; face z would be needed for the latter (or the
  center-line approximation, as TOUGH's `BETAX`). This is feasible without a water table.
- **HFB:** a barrier on a face generalizes naturally (multiply the face conductance), but
  the current input is layer/cell-pair based and the code assumes horizontal connections.
  This is modest work.
- **CSUB:** vertical compaction of layered interbeds. Not meaningful without a vertical
  structure.
- **UZF:** 1D vertical columns. Not meaningful.
- **PRT:** Pollock-family methods need prism/quad cells. Polyhedral tracking would need a
  new method. Not addressed by the fork and not guarded (see §4).

## 4. Integration and implementation findings

These are static findings. "Likely" means inferred from the code and not run.

1. **Parallel / interface-model path is unguarded (likely crash).** Exchanges whose models
   live on different processes always use an interface model
   (`DisConnExchange.f90:481`). The interface model's DISU is built by
   `GridConnection%getDiscretization` (`GridConnection.f90:1007` onward) with
   `ifull3d = 0`, placeholder `top/bot/area = -huge`, and no `zc`. Connections are copied
   with `IHC=3`. With anisotropic K, `hy_eff` → `connection_normal(IHC=3)` dereferences
   `zc`, which is unallocated. Even with isotropic K, `GwfGwfConnection.f90:651` calls
   `calcSatThickness(.., ihc=3)` and multiplies the face area by it (horizontal branch).
   None of the `FULL_3D` guards run in the interface model. The issue's claim that the
   full suite passes excludes `--parallel`.
2. **GWF-GWF exchange direction ignores model origins.** `center_vector`
   (`exg-gwfgwf.f90:1478`) uses model-local `xc/yc` from `get_cell_center`. Classic
   exchange/XT3D code transforms to global coordinates (`dis_transform_xy`). Two `FULL_3D`
   models with different `XORIGIN/YORIGIN/ANGROT` get the wrong K direction. There is also
   no z origin concept. Exchange faces are skipped in the specific-discharge fit
   (`exg-gwfgwf.f90:791`), so model-boundary cells get a biased velocity, which is also
   what DSP uses.
3. **Package allow-list only covers GWF.** `check_full_3d` runs in GWF only. GWT/GWE
   models on a `FULL_3D` DISU have no package allow-list: every GWT/GWE package is
   accepted, whether or not it was considered. SWF-GWF exchanges and PRT models reading a
   `FULL_3D` DISU are not checked either.
4. **Binary grid file is misleading.** The DFN text says the GRB file "is not written", but
   the code writes it with a synthesized cube `TOP/BOT` and no marker that the grid is
   `FULL_3D` (`Disu.f90:1586`). FloPy's `MfGrdFile`, ModelMuse, PRT post-processing and
   any API/BMI client reading `TOP/BOT/AREA` will reconstruct a different grid without
   warning. A new grid type in the GRB header (and a BMI-visible flag) is needed.
5. **API/BMI exposes the placeholders.** `get_value("…/DIS/TOP")` returns 1,
   `AREA` returns the volume. Coupled codes doing area-based exchanges will be silently
   wrong.
6. **Redundant, unvalidated geometry.** The VTU geometry file contains the actual
   polyhedra (points + faces), but they are used only for output. Volume, centers, face
   areas and distances are not cross-checked against them. If the geometry is available,
   the model should *compute* these quantities (as DISV does from vertices), or at least
   validate the user values. That would remove an entire class of input errors and provide
   the face normals that §3.1 needs.
7. **Two ways to give cell centers.** `XCELL/YCELL` griddata arrays in addition to
   `CELL2D` for classic DISU. Minor, but it is DFN surface area that persists.
8. **Tests do not exercise polyhedra.** Every test grid is a randomly stretched
   *hexahedral* grid (`test_gwf_disu_full3d.py::build_grid`), compared against classic
   DISU on the same grid. That verifies the placeholder mapping is self-consistent but
   never runs a non-prismatic cell, a non-orthogonal face, or anything a prism grid could
   not already represent. The Voronoi benchmark mentioned in the issue is not in the PR.
   The XT3D test (rotated tensor, linear head) is the most meaningful one.
9. **Scattered capability checks.** About 60 `is_full_3d()`/`ifull3d`/`C3D_FULL`
   branches in 14 files fence the feature off, while `dis%top/bot/area` is referenced
   ~250 times in 37 files (GWF 15, SWF 5, GWT/GWE/TSP 9, PRT 3, export 2, exchanges 3, …).
   Every future feature touching those arrays has to remember the exemption. This is the
   maintenance cost that matters most.
10. **Out-of-scope content.** The VTU writer targets a format from a specific external
    tool (VORO2MESH `.v3d`) and duplicates what FloPy's VTK export does for other grids.
    MF6 has no other native VTK output. The CMake/`build_vs.bat` commit, the PRT
    `CellUtil` change and the UZF import change are unrelated to the feature. The single
    ~2,900-line commit should be split.

The audit found nothing suspicious in the build scripts or VTU reader (local paths only, no
network or `execute_process`). The `.v3d` reader validates point indices and fails cleanly
on malformed input.

## 5. Framework implications

Your hunch about framework refactoring is right, but the work is more contained than a
rewrite. The core issue is that `DisBaseType` exposes `top`, `bot` and `area` as public
arrays that packages use directly as geometry. Clean polyhedral support needs packages to
ask the discretization geometric questions instead:

| Query | Prism default | Polyhedral |
|---|---|---|
| `get_cell_volume(n)` (exists) | `A·(top−bot)` | `V_n` |
| `get_sat_volume(n, h)` / `dVdh` | prism formulas | hypsometry (or error if confined-only) |
| `get_face_area(ipos, h)` | `hwva·thksat` or `hwva` | `hwva` (confined) / `A_f(h)` |
| `get_face_normal(ipos)` | from `ihc`/`anglex` | input or from faces |
| `get_cell_center(n)` incl. z | `(xc, yc, (top+bot)/2)` | `(xc, yc, zc)` |
| `get_exposed_area(n)` (RCH/EVT) | `area` | input |
| capability flags (`has_vertical_extent`, `supports_unconfined`) | true | false |

Most call sites in STO, MST/EST, DSP/CND, FMI, ADV and the BND packages could move to
these with no behavior change for DIS/DISV/DISU. That refactor is useful on its own: it
cleans up the existing `ihc==0` versus horizontal branching in DSP/CND/NPF/XT3D. After it,
a polyhedral discretization becomes a new `DisBaseType` subclass plus a short list of
packages that check `supports_unconfined()` once, instead of 60 `is_full_3d()` branches.
The interface-model `GridConnection` would also need to carry `zc` (and normals) in
virtual model data, which fixes finding 1 properly.

**Option vs. new package.** Under `FULL_3D` the required arrays change, `TOP/BOT` become
forbidden, `AREA` and `HWVA` change meaning, `IHC` gets a value that is legal only in this
mode, and `NVERT/CELL2D` are forbidden. That makes it a different discretization in all but
name. A new type (working name `DISP` or `DIS3D`) gives a clean DFN, its own GRB grid type,
and a clear "not supported with DISP" story per package. It can still reuse
`ConnectionsType` and much of `DisuType` by extension.

## 6. Feasibility summary

| Scope | Feasibility | Effort | Notes |
|---|---|---|---|
| Confined GWF, isotropic or K-orthogonal anisotropy | High | Low (fork already does it) | Needs the §5 refactor to be maintainable |
| + XT3D on non-orthogonal meshes | Medium | Medium | Needs face normals in input |
| + GWT/GWE (ADV, DSP/CND, MST/EST, SSM) | High | Low–medium | Works now by placeholder; should go through geometry API |
| + GWF-GWF / GWT-GWT / GWE-GWE exchanges, serial and parallel | Medium | Medium | Interface model and virtual data need z and normals; global coordinates |
| + BUY, HFB | Medium | Medium | Face elevations / face-based barrier input |
| + water table: convertible cells, Newton, drying/rewetting | Low–medium | High | Hypsometric functions from vertices; new derivative forms throughout NPF/STO/advanced packages |
| UZF, CSUB | Not meaningful | – | Inherently vertical/layered |
| PRT | Low | High | New tracking method for polyhedra |
| Post-processing (GRB, FloPy, ModelMuse) | High | Medium | New GRB grid type; FloPy grid class; VTK export in FloPy |

## 7. Suggested response to the issue author

1. **Interest:** yes for confined polyhedral flow and transport. The TOUGH-style
   formulation is credible and the verification so far is encouraging.
2. **Exposure:** prefer a new discretization package over a DISU option, preceded by (or
   bundled with) a geometry-accessor refactor of `DisBaseType` so consumers stop reading
   `top/bot/area` directly. Avoid placeholder geometry.
3. **Input:** accept face geometry (vertices + faces, or at minimum face normals) so the
   model can compute or validate volume, centers, areas and distances, and so XT3D gets
   correct normals. This also keeps the door open to unconfined support later.
4. **First PR scope:** discretization + NPF/STO (confined) + CHD/WEL/GHB/DRN/RIV, with
   tests on actual polyhedral meshes (a small committed 3D Voronoi mesh and the analytical
   benchmark), and explicit errors for everything else, including GWT/GWE, PRT, SWF, and
   parallel until supported. Follow-ups: transport, exchanges (serial and parallel), RCH/EVT,
   advanced packages.
5. **Drop from the PR:** VTU writer (do this in FloPy), CMake build files, and the
   unrelated PRT/UZF edits.
6. **Output:** new GRB grid type, and FloPy support, before release.
