"""
Tracking times and releases at a computed end of tracking.

Some times at which tracking ends are sums computed by PRT, which can
come out 1 ulp below the time a user would type for them. A given time
within roundoff after such an end is treated as the end. All cases use
ten stress periods of 0.1, which end at 0.9999999999999999, not 1.0,
and run without EXTEND_TRACKING.

    end:      a tracking time at the end of the simulation (1.0). It's
              reported once, in the last time step, and the final time
              step and termination events follow at the same time. A
              tracking time well after the end (1.5) isn't reported,
              and a warning says so.
    stoptt:   a tracking time at the stop time from STOPTRAVELTIME. A
              release at 0.35 with travel time 0.1 stops at
              0.44999999999999996; the tracking time 0.45 is reported,
              and the termination follows at the same time.
    release:  a release time at the end of the simulation (1.0). It
              isn't skipped: the particle is released and terminated at
              1.0, as it would be if the simulation ended exactly
              there. A release time well after the end (1.5) is
              skipped, with a warning.

In each case the particle is still moving when tracking ends. Each case
runs on a DIS grid with Pollock's method, and on a DISV grid with
Pollock's method and with the ternary method.
"""

from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "endtt"
nrow = ncol = 10
top, botm = 1.0, [0.0]
porosity = 0.1
nper, perlen = 10, 0.1
PERIODDATA = [(perlen, 1, 1.0)] * nper
TEND = 1.0  # the end of the simulation, as a user would type it
TRELEASE, STOPTRAVELTIME, TSTOP = 0.35, 0.1, 0.45
CASES = {
    # case: (release times, tracking times, stop travel time)
    "end": ([0.0], [0.5, TEND, 1.5], None),
    "stoptt": ([TRELEASE], [0.4, TSTOP, 0.5], STOPTRAVELTIME),
    "release": ([0.0, TEND, 1.5], None, None),
}
RELEASE_CELL = (0, 2, 2)
# off the cell center and the grid's line of symmetry, so the particle
# doesn't track exactly along an internal edge between ternary subcells
RELEASE_POINT = (2.4, 7.3, 0.5)
GRIDS = {"dis": "pollock", "disv": "pollock", "disvt": "ternary"}
cases = [f"{simname}{case[0]}_{grid}" for case in CASES for grid in GRIDS]

RELEASE, EXIT, TIMESTEP, TERMINATE, USERTIME = 0, 1, 2, 3, 5


def parse(name):
    code, grid = name[len(simname) :].split("_")
    case = next(c for c in CASES if c[0] == code)
    return case, grid


def float_exactly(sim):
    sim.simulation_data.float_precision = 17
    sim.simulation_data.float_characters = 26


def add_dis(model, grid):
    gwf = isinstance(model, flopy.mf6.ModflowGwf)
    if grid == "dis":
        cls = flopy.mf6.ModflowGwfdis if gwf else flopy.mf6.ModflowPrtdis
        cls(
            model,
            nlay=1,
            nrow=nrow,
            ncol=ncol,
            delr=1.0,
            delc=1.0,
            top=top,
            botm=botm,
        )
    else:
        cls = flopy.mf6.ModflowGwfdisv if gwf else flopy.mf6.ModflowPrtdisv
        cls(
            model,
            **get_disv_kwargs(
                nlay=1,
                nrow=nrow,
                ncol=ncol,
                delr=1.0,
                delc=1.0,
                tp=top,
                botm=botm,
            ),
        )


def cellid(grid, kij):
    k, i, j = kij
    return (k, i, j) if grid == "dis" else (k, i * ncol + j)


def build_gwf_sim(name, ws, mf6, grid):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6, sim_ws=ws)
    float_exactly(sim)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=PERIODDATA)
    gwfname = get_model_name(name, "gwf")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    add_dis(gwf, grid)
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, save_saturation=True)
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            [cellid(grid, (0, 0, 0)), top],
            [cellid(grid, (0, nrow - 1, ncol - 1)), 0.0],
        ],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowIms(sim)
    return sim


def build_prt_sim(name, gwf_ws, prt_ws, mf6, grid, case):
    releasetimes, tracktimes, stoptraveltime = CASES[case]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6, sim_ws=prt_ws)
    float_exactly(sim)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=PERIODDATA)
    prtname = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)
    add_dis(prt, grid)
    flopy.mf6.ModflowPrtmip(prt, porosity=porosity)
    x, y, z = RELEASE_POINT
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=1,
        packagedata=[[0, cellid(grid, RELEASE_CELL), x, y, z]],
        nreleasetimes=len(releasetimes),
        releasetimes=[(t,) for t in releasetimes],
        stoptraveltime=stoptraveltime,
        dev_forceternary=GRIDS[grid] == "ternary",
        extend_tracking=False,
    )
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        trackcsv_filerecord=[f"{prtname}.trk.csv"],
        ntracktimes=len(tracktimes) if tracktimes else None,
        tracktimes=[(t,) for t in tracktimes] if tracktimes else None,
    )
    gwfname = get_model_name(name, "gwf")
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", Path(gwf_ws) / f"{gwfname}.hds"),
            ("GWFBUDGET", Path(gwf_ws) / f"{gwfname}.bud"),
        ],
    )
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{prtname}.ems")
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_models(test, grid, case):
    mf6 = test.targets["mf6"]
    gwf_sim = build_gwf_sim(test.name, test.workspace, mf6, grid)
    prt_sim = build_prt_sim(
        test.name, test.workspace, test.workspace / "prt", mf6, grid, case
    )
    return gwf_sim, prt_sim


def read_lst(test):
    lst = (test.workspace / "prt" / "mfsim.lst").read_text().upper()
    return " ".join(lst.split())  # messages are wrapped


def check_end(trk, test):
    # the tracking time at the typed end is reported once, in the last step
    at = trk[(trk.ireason == USERTIME) & (trk.t == TEND)]
    assert len(at) == 1, f"tracking time {TEND} reported {len(at)} times:\n{trk}"
    assert at.kper.iloc[0] == nper

    # the particle was still moving: the final time step and termination
    # events are reported at the tracking time, not before it
    last = trk.tail(2)
    assert list(last.ireason) == [TIMESTEP, TERMINATE], trk.tail(4)
    assert (last.t == TEND).all(), trk.tail(4)

    # the tracking time well after the end isn't reported, with a warning
    assert (trk.t <= TEND).all(), trk.tail(4)
    assert "1 TRACKING TIME(S) AFTER THE END OF THE SIMULATION" in read_lst(test)


def check_stoptt(trk, test):
    # the premise: the computed stop time is just below the typed one
    tstop = TRELEASE + STOPTRAVELTIME
    assert tstop < TSTOP and np.nextafter(tstop, np.inf) == TSTOP

    # the tracking time at the typed stop time is reported once, and the
    # termination follows at the same time; the later one isn't reported
    at = trk[(trk.ireason == USERTIME) & (trk.t == TSTOP)]
    assert len(at) == 1, f"tracking time {TSTOP} reported {len(at)} times:\n{trk}"
    last = trk.iloc[-1]
    assert last.ireason == TERMINATE and last.t == TSTOP, trk.tail(4)
    assert (trk.t <= TSTOP).all(), trk.tail(4)


def check_release(trk, test):
    # the release at the typed end isn't skipped: the particle is
    # released and terminated there, in the last time step, with a time
    # step event in between, as for a release exactly at the end
    late = trk[trk.trelease == TEND]
    assert list(late.ireason) == [RELEASE, TIMESTEP, TERMINATE], f"\n{late}"
    assert (late.t == TEND).all() and (late.kper == nper).all(), late

    # the first particle ends at the end of the simulation as computed
    first = trk[trk.trelease == 0.0]
    assert first.ireason.iloc[-1] == TERMINATE
    assert first.t.iloc[-1] == np.nextafter(TEND, -np.inf), first.tail(2)

    # the release well after the end is skipped, with a warning
    assert not (trk.trelease == 1.5).any()
    lst = read_lst(test)
    # only the release at 1.5 is skipped
    assert lst.count("SKIPPING RELEASE TIME") == 1 and "(T=1.5" in lst


def check_output(test, grid, case):
    # the premise: TDIS's computed end is just below the typed end
    tsim = 0.0
    for _ in range(nper):
        tsim += perlen
    assert tsim < TEND and np.nextafter(tsim, np.inf) == TEND

    prtname = get_model_name(test.name, "prt")
    trk = pd.read_csv(test.workspace / "prt" / f"{prtname}.trk.csv")
    assert trk.ireason.iloc[0] == RELEASE
    for _, p in trk.groupby("trelease"):
        assert (np.diff(p.t) >= 0).all(), f"times decrease:\n{p.tail(6)}"
    {"end": check_end, "stoptt": check_stoptt, "release": check_release}[case](
        trk, test
    )


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    case, grid = parse(name)
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(t, grid, case),
        check=lambda t: check_output(t, grid, case),
        targets=targets,
        compare=None,
    )
    test.run()
