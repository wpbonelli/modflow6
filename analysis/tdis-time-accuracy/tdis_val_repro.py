"""Write small GWF models that show the check_tdis_timing zero-PERLEN bug.

check_tdis_timing returns at the first period with PERLEN = 0, so later periods
are not validated. Each case has a zero-length steady-state first period
followed by an invalid second period, except the control, whose first period
has PERLEN = 1 and so is validated normally.

Usage: python -I tdis_val_repro.py <workspace> <path to mf6>
See README.md for running the models and the expected results.
"""

import sys
from pathlib import Path

import flopy

root = Path(sys.argv[1])
mf6 = sys.argv[2]
cases = {
    "neg_perlen_after_zero": [(0.0, 1, 1.0), (-10.0, 1, 1.0)],
    "neg_perlen_control": [(1.0, 1, 1.0), (-10.0, 1, 1.0)],
    "neg_tsmult_after_zero": [(0.0, 1, 1.0), (10.0, 3, -2.0)],
    "zero_nstp_after_zero": [(0.0, 1, 1.0), (10.0, 0, 1.0)],
}
for name, pd in cases.items():
    ws = root / name
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=mf6)
    flopy.mf6.ModflowTdis(sim, nper=len(pd), perioddata=pd)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=2)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfsto(gwf, steady_state={0: True}, transient={1: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 1.0]])
    sim.write_simulation(silent=True)
