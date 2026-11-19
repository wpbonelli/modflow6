"""
Ternary exit times in subcells where the particle's beta coordinate is
constant.

Flow runs along a single row of DISV cells, with injection (or extraction)
wells in the interior cells, so the x velocity changes linearly across each
cell and the y velocity is zero. In the ternary method, the subcells on the
top and bottom of each cell then have their outer face parallel to the
flow, so the particle's beta coordinate stays constant there while its
alpha velocity varies. Those subcells take a closed-form exit time, which
had misplaced parentheses and gave a wrong or negative time ("texit is
negative (unexpected)").

The velocity is linear in x, so Pollock's method is exact. Particles are
released off the cell's centerline, in the top subcell, and their ternary
exit times must match Pollock's.
"""

from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest
from flopy.utils.gridutil import get_disv_kwargs

nlay, nrow, ncol = 1, 1, 8
top, botm = 1.0, [0.0]
porosity = 0.1
# release points in cell 2 (x 1..2), above the centerline y = 0.5
RELEASE_POINTS = [(1.6, y, 0.5) for y in (0.55, 0.7, 0.85, 0.95)]
EXIT = 1
# well rate in each interior cell: positive accelerates the flow, negative
# decelerates it
RATES = {"accelerating": 0.5, "decelerating": -0.5}


def disv():
    return get_disv_kwargs(nlay, nrow, ncol, 1.0, 1.0, top, botm)


def run_gwf(ws, mf6, rate):
    sim = flopy.mf6.MFSimulation(sim_name="gwf", exe_name=mf6, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **disv())
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=0, save_specific_discharge=True, save_saturation=True
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0), 10.0], [(0, ncol - 1), 0.0]]
    )
    flopy.mf6.ModflowGwfwel(
        gwf, stress_period_data=[[(0, j), rate] for j in range(1, ncol - 1)]
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord="gwf.bud",
        head_filerecord="gwf.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowIms(sim)
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]


def run_prt(ws, gwf_ws, mf6, ternary):
    sim = flopy.mf6.MFSimulation(sim_name="prt", exe_name=mf6, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim)
    prt = flopy.mf6.ModflowPrt(sim, modelname="prt")
    flopy.mf6.ModflowPrtdisv(prt, **disv())
    flopy.mf6.ModflowPrtmip(prt, porosity=porosity)
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename="prt.prp",
        nreleasepts=len(RELEASE_POINTS),
        packagedata=[
            [i, (0, 1), x, y, z] for i, (x, y, z) in enumerate(RELEASE_POINTS)
        ],
        nreleasetimes=1,
        releasetimes=[(0.0,)],
        dev_forceternary=ternary,
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(prt, pname="oc", trackcsv_filerecord=["prt.trk.csv"])
    gwf_ws = Path(gwf_ws).absolute()
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_ws / "gwf.hds"),
            ("GWFBUDGET", gwf_ws / "gwf.bud"),
        ],
    )
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename="prt.ems")
    sim.register_solution_package(ems, [prt.name])
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    assert success, "\n".join(buff[-5:])
    trk = pd.read_csv(Path(ws) / "prt.trk.csv")
    return trk[trk.ireason == EXIT].reset_index(drop=True)


@pytest.mark.parametrize("flow", list(RATES))
def test_ternary_const_beta(function_tmpdir, targets, flow):
    mf6 = targets["mf6"]
    gwf_ws = function_tmpdir / "gwf"
    run_gwf(gwf_ws, mf6, RATES[flow])
    pollock = run_prt(function_tmpdir / "pollock", gwf_ws, mf6, ternary=False)
    ternary = run_prt(function_tmpdir / "ternary", gwf_ws, mf6, ternary=True)

    # premise: particles cross the row along straight lines, y unchanged
    assert pollock.irpt.nunique() == len(RELEASE_POINTS)
    for irpt, p in pollock.groupby("irpt"):
        assert np.allclose(p.y, RELEASE_POINTS[irpt - 1][1], rtol=0, atol=1e-12)

    # ternary exits the same cells as Pollock, at the same times
    cols = ["irpt", "icell"]
    assert ternary[cols].equals(pollock[cols]), (
        f"pollock:\n{pollock}\nternary:\n{ternary}"
    )
    assert np.allclose(ternary.t, pollock.t, rtol=0, atol=1e-9), pd.DataFrame(
        {"irpt": ternary.irpt, "icell": ternary.icell, "dt": ternary.t - pollock.t}
    )
