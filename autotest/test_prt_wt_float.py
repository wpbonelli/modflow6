"""
Tests particle tracking across the water-table surface.

A 2-layer, 1-row, 3-column GWF model with:
  - Newton solver and convertible cells (icelltype=1)
  - CHD on the left (col 0, h=4.5) and right (col 2, h=3.5)
  - EVT on the middle cell (col 1) with IFLOWFACE=-1 to route
    evapotranspiration through the top face, creating upward
    velocity that drives particles to the water-table surface
  - The right CHD uses IFLOWFACE=3 to mark the east face as a
    PRT boundary so that floating particles terminate there

Three cases exercise each dry_tracking_method at the water table:

  DROP  — particle "floats" to the water table and continues
           tracking laterally until it exits as TERM_BOUNDARY
  STOP  — particle terminates (TERM_INACTIVE) at the water table
  STAY  — particle beaches at the water table for the remainder
           of the simulation and terminates as TERM_NO_EXITS
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "prtwt"
cases = [
    simname + "_drop",
    simname + "_stop",
    simname + "_stay",
]

# grid
nlay, nrow, ncol = 2, 1, 3
delr, delc = 10.0, 10.0

# layer tops/bots: layer 0 dry upper; layer 1 water-table lower
top = [10.0, 5.0]
bot = [[5.0] * ncol, [0.0] * ncol]
tops_arr = np.array([[10.0] * ncol] * nrow)
botm_arr = np.array([[[5.0] * ncol] * nrow, [[0.0] * ncol] * nrow])

# hydraulics
K = 0.1
K33 = 0.01
porosity = 0.1
strt = 4.0

# EVT
et_surface = 5.0
et_rate = 0.001  # m/day; drives upward velocity ~0.008 m/day in col 1
et_depth = 5.0

# CHD heads
h_left = 4.5
h_right = 3.5

# stress period: single steady-state period long enough for particle to exit
nper = 1
perlen = 400.0
nstp = 1
period_data = [(perlen, nstp, 1.0)]


def build_gwf_sim(name, gwf_ws, mf6):
    sim = flopy.mf6.MFSimulation(
        sim_name=name, exe_name=mf6, version="mf6", sim_ws=gwf_ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=nper, perioddata=period_data)

    gwf_name = get_model_name(name, "gwf")
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_name,
        newtonoptions="NEWTON UNDER_RELAXATION",
        save_flows=True,
    )

    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1e-6,
        outer_maximum=100,
        under_relaxation="DBD",
        under_relaxation_theta=0.7,
        under_relaxation_kappa=0.001,
        under_relaxation_gamma=0.0,
        under_relaxation_momentum=0.0,
        inner_maximum=100,
        inner_dvclose=1e-7,
        rcloserecord=0.001,
        linear_acceleration="BICGSTAB",
        relaxation_factor=0.97,
        pname=gwf_name,
    )

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=tops_arr,
        botm=botm_arr,
        pname="dis",
    )

    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,
        k=K,
        k33=K33,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
    )

    # Left CHD: h=4.5, no IFLOWFACE (not a PRT boundary)
    chd_left = [((1, 0, 0), h_left)]
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: chd_left},
        save_flows=True,
        pname="chd_left",
    )

    # Right CHD: h=3.5, IFLOWFACE=3 (east face = PRT boundary)
    chd_right = [((1, 0, 2), h_right, 3)]
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: chd_right},
        auxiliary=["iflowface"],
        save_flows=True,
        pname="chd_right",
    )

    # EVT on middle cell (col 1), IFLOWFACE=-1 routes ET through top face
    # stress_period_data format: (cellid, surface, rate, depth, [aux...])
    evt_spd = [((1, 0, 1), et_surface, et_rate, et_depth, -1)]
    flopy.mf6.ModflowGwfevt(
        gwf,
        stress_period_data={0: evt_spd},
        auxiliary=["iflowface"],
        save_flows=True,
        pname="evt",
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=gwf_name + ".cbb",
        head_filerecord=gwf_name + ".hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


def build_prt_sim(name, gwf, prt_ws, mf6, dry_tracking_method):
    sim = flopy.mf6.MFSimulation(
        sim_name=name, exe_name=mf6, version="mf6", sim_ws=prt_ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="days", nper=nper, perioddata=period_data)

    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

    flopy.mf6.ModflowPrtdis(
        prt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=tops_arr,
        botm=botm_arr,
        pname="dis",
    )

    flopy.mf6.ModflowPrtmip(prt, porosity=porosity, pname="mip")

    # Release a single particle in layer 1 (0-based), row 0, col 1
    # at z=3.5m — below the expected water table (~4m) so the upward
    # EVT velocity drives it to the water-table surface.
    x_mid = gwf.modelgrid.xyzcellcenters[0][0, 1]
    y_mid = gwf.modelgrid.xyzcellcenters[1][0, 1]
    z_rel = 3.5
    prp_data = [[0, (1, 0, 1), x_mid, y_mid, z_rel]]

    flopy.mf6.ModflowPrtprp(
        prt,
        nreleasepts=1,
        packagedata=prp_data,
        nreleasetimes=1,
        releasetimes=[(0.0,)],
        dry_tracking_method=dry_tracking_method,
        pname="prp",
    )

    prt_bud = prt_name + ".bud"
    prt_trk = prt_name + ".trk"
    prt_csv = prt_name + ".csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        budget_filerecord=[prt_bud],
        track_filerecord=[prt_trk],
        trackcsv_filerecord=[prt_csv],
        saverecord=[("BUDGET", "ALL")],
        pname="oc",
    )

    import os

    gwf_ws = gwf.model_ws
    rel_gwf = os.path.relpath(gwf_ws, start=prt_ws)
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", f"{rel_gwf}/{gwf.name}.hds"),
            ("GWFBUDGET", f"{rel_gwf}/{gwf.name}.cbb"),
        ],
    )

    flopy.mf6.ModflowEms(sim, pname="ems", filename=prt_name + ".ems")

    return sim


def build_models(idx, test):
    name = cases[idx]
    dry_tracking_method = name.split("_")[-1]  # "drop", "stop", or "stay"
    gwf_sim = build_gwf_sim(name, test.workspace / "gwf", test.targets["mf6"])
    prt_sim = build_prt_sim(
        name,
        gwf_sim.get_model(),
        test.workspace / "prt",
        test.targets["mf6"],
        dry_tracking_method=dry_tracking_method,
    )
    return gwf_sim, prt_sim


def check_output(idx, test, snapshot):
    name = cases[idx]
    prt_name = get_model_name(name, "prt")
    prt_ws = test.workspace / "prt"
    csv_path = prt_ws / (prt_name + ".csv")

    pls = pd.read_csv(csv_path)
    method = name.split("_")[-1]

    # one particle, so exactly one termination event
    term = pls[pls.ireason == 3]
    assert len(term) == 1, f"expected 1 termination, got {len(term)}"

    if method == "stop":
        # STOP: particle terminated at water-table surface (TERM_INACTIVE=7)
        assert term.iloc[0].istatus == 7, (
            f"STOP: expected istatus=7, got {term.iloc[0].istatus}"
        )
        assert len(pls[pls.ireason == 0]) == 1  # exactly one release
        # particle exits at the water table (ireason=1 cell-exit) then terminates
    elif method == "stay":
        # STAY: particle beached at water table, terminates as TERM_NO_EXITS=5
        assert term.iloc[0].istatus == 5, (
            f"STAY: expected istatus=5, got {term.iloc[0].istatus}"
        )
    elif method == "drop":
        # DROP: particle floats laterally and exits as TERM_BOUNDARY=2
        assert term.iloc[0].istatus == 2, (
            f"DROP: expected istatus=2 (TERM_BOUNDARY), got {term.iloc[0].istatus}"
        )
        # multiple records: release + tracking + termination
        assert len(pls) > 2, f"DROP: expected >2 records, got {len(pls)}"

    actual = pls.drop(["name", "icell"], axis=1, errors="ignore").round(4)
    assert snapshot == actual.to_records(index=False)


@pytest.mark.snapshot
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, array_snapshot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, array_snapshot),
        targets=targets,
        compare=None,
    )
    test.run()
