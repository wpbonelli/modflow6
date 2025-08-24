# PRT cyclic pathline case 2 patch plan

game plan for fixing the second vertical-motion cycle (tracking entering a loop) case identified recently.

## problems

1. at tracking time, PRT sees the water table as calculated for the end of the time step. it doesn't see change
   from the previous step.
2. therefore, interpolating vertical velocities only valid when >0 flow thru both faces, i.e. saturated cells.
3. with newton, dry/active cells can get infinitesimal flows, e.g. up thru top of partially saturated cells.
   these are meaningless numerical noise but they can currently produce unrealistic exit face results like
   a particle "leaping" from the water table into a dry cell.

## solutions

no easy fixes for 1 and 2. deal with the slight unphysicality for now.
later, account for water table change over time step. modify velocity calculation accordingly.

for 3, clamp particle to water table. don't let particles follow upwards flow in partially saturated cells.
if `DRY_TRACKING_METHOD DROP` the water table is like a one-way membrane permeable only from above.
if it recedes it will will drag the particle down with it.
if `DRY_TRACKING_METHOD STAY` the water table can "strand" the particle, leaving it active but stationary.
if `DRY_TRACKING_METHOD STOP`, stop.
this can all be done by a single well-chosen procedure override in the cell-level method.
in the current first draft it is done separately in the two subcell-level methods.

## gameplan

solve 3. punt 1-2 for now.

also, introduce two new particle events, "drop" and "strand", occurring when a particle enters the water
with `DROP` configured, or when the water falls below it with `STAY`, respectively.

this general direction is indicated in code comments.

