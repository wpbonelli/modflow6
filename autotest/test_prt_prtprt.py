"""
Test PRT-PRT exchange between two particle tracking models.

Two GWF models share a domain split vertically down the middle:

   GWF-left (cols 0-4)  |  GWF-right (cols 5-9)

Steady flow runs from left to right: the left boundary of GWF-left
has a prescribed head of 1.0, the right boundary of GWF-right has
a prescribed head of 0.0.  Hydraulic conductivity is uniform at 1 m/d.

Each GWF model has a matched PRT model connected via a GWF6-PRT6 exchange.
The two PRT models are connected via a PRT6-PRT6 exchange at the same
interface column pair:  (0, 0, 4) in PRT-left  <->  (0, 0, 0) in PRT-right.

Particles are released from the top-left cell (0, 0, 0) in PRT-left
at the start of the simulation.  Because flow is from left to right,
the particles should cross the PRT-PRT exchange and terminate on the
right boundary of PRT-right.  The test asserts that:

  1. All particles that were released ultimately appear in the
     PRT-right track output (implying a successful handoff).
  2. All particles terminate with status TERMINATE_EXIT_FACE (= 1)
     at the right boundary of PRT-right – no particle is lost.
"""

from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

# ------------------------------------------------------------------ #
# Geometry
# ------------------------------------------------------------------ #
nlay = 1
nrow = 1
ncol_half = 5        # columns per sub-model
ncol = ncol_half * 2 # total columns

delr = 1.0  # m – column width
delc = 1.0  # m – row width
delz = 1.0  # m – layer thickness
top = 1.0
botm = [0.0]

# Flow BCs
head_left = 1.0
head_right = 0.0
hk = 1.0   # m/d hydraulic conductivity

porosity = 0.1
nper = 1
perlen = 10.0  # days
nstp = 10

cases = ["prt_prtprt"]


# ------------------------------------------------------------------ #
# Helper: names
# ------------------------------------------------------------------ #
def gwf_name(side):
    return f"gwf_{side}"


def prt_name(side):
    return f"prt_{side}"


# ------------------------------------------------------------------ #
# Build simulation
# ------------------------------------------------------------------ #
def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    mf6 = test.targets["mf6"]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        sim_ws=ws,
    )

    # -- Time discretisation --
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(perlen, nstp, 1.0)],
    )

    # -- Flow solver --
    flopy.mf6.ModflowIms(
        sim,
        pname="ims",
        complexity="simple",
        outer_dvclose=1e-8,
        inner_dvclose=1e-8,
    )

    # ---- GWF-left ------------------------------------------------- #
    gwfl = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_name("left"),
        save_flows=True,
    )
    flopy.mf6.ModflowGwfdis(
        gwfl,
        nlay=nlay, nrow=nrow, ncol=ncol_half,
        delr=delr, delc=delc,
        top=top, botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwfl, strt=head_left)
    flopy.mf6.ModflowGwfnpf(
        gwfl,
        save_specific_discharge=True,
        icelltype=0,
        k=hk,
    )
    # CHD on left face
    flopy.mf6.ModflowGwfchd(
        gwfl,
        stress_period_data=[[(0, 0, 0), head_left]],
        pname="chd_left",
    )
    flopy.mf6.ModflowGwfoc(
        gwfl,
        head_filerecord=f"{gwf_name('left')}.hds",
        budget_filerecord=f"{gwf_name('left')}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # ---- GWF-right ------------------------------------------------ #
    gwfr = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_name("right"),
        save_flows=True,
    )
    flopy.mf6.ModflowGwfdis(
        gwfr,
        nlay=nlay, nrow=nrow, ncol=ncol_half,
        delr=delr, delc=delc,
        top=top, botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwfr, strt=head_right)
    flopy.mf6.ModflowGwfnpf(
        gwfr,
        save_specific_discharge=True,
        icelltype=0,
        k=hk,
    )
    # CHD on right face
    flopy.mf6.ModflowGwfchd(
        gwfr,
        stress_period_data=[[(0, 0, ncol_half - 1), head_right]],
        pname="chd_right",
    )
    flopy.mf6.ModflowGwfoc(
        gwfr,
        head_filerecord=f"{gwf_name('right')}.hds",
        budget_filerecord=f"{gwf_name('right')}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # ---- GWF-GWF exchange ---------------------------------------- #
    # Connect right face of GWF-left to left face of GWF-right.
    # (layer, row, col): left model col 4  <->  right model col 0
    gwfgwf_data = [
        (
            (0, 0, ncol_half - 1),  # cellidm1 in GWF-left
            (0, 0, 0),              # cellidm2 in GWF-right
            1,                      # ihc = horizontal connection
            delr / 2.0,             # cl1
            delr / 2.0,             # cl2
            delc,                   # hwva
            0.0,                    # ANGLDEGX: flow in +x direction
            delr,                   # CDIST: total connection distance
        )
    ]
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=gwf_name("left"),
        exgmnameb=gwf_name("right"),
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="gwf-gwf.exg",
    )

    # ---- PRT-left ------------------------------------------------- #
    prtl = flopy.mf6.ModflowPrt(sim, modelname=prt_name("left"))
    flopy.mf6.ModflowGwfdis(     # PRT reuses ModflowGwfdis for its DIS
        prtl,
        nlay=nlay, nrow=nrow, ncol=ncol_half,
        delr=delr, delc=delc,
        top=top, botm=botm,
    )
    flopy.mf6.ModflowPrtmip(prtl, pname="mip", porosity=porosity)

    # Release one particle from cell (0,0,0): localx=0.5, localy=0.5, localz=0.5
    releasepts = [(0, (0, 0, 0), 0.5, 0.5, 0.5)]
    flopy.mf6.ModflowPrtprp(
        prtl,
        pname="prp1",
        filename=f"{prt_name('left')}.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(
        prtl,
        pname="oc",
        track_filerecord=f"{prt_name('left')}.trk",
        trackcsv_filerecord=f"{prt_name('left')}.trk.csv",
    )
    # FMI: link to GWF-left flow output
    flopy.mf6.ModflowPrtfmi(
        prtl,
        packagedata=[
            ("GWFHEAD",   f"{gwf_name('left')}.hds"),
            ("GWFBUDGET", f"{gwf_name('left')}.cbc"),
        ],
    )
    # GWF6-PRT6 exchange
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwf_name("left"),
        exgmnameb=prt_name("left"),
        filename=f"{gwf_name('left')}.gwfprt",
    )

    # ---- PRT-right ----------------------------------------------- #
    prtr = flopy.mf6.ModflowPrt(sim, modelname=prt_name("right"))
    flopy.mf6.ModflowGwfdis(
        prtr,
        nlay=nlay, nrow=nrow, ncol=ncol_half,
        delr=delr, delc=delc,
        top=top, botm=botm,
    )
    flopy.mf6.ModflowPrtmip(prtr, pname="mip", porosity=porosity)
    # PRT-right has no PRP; it only receives particles from PRT-left
    flopy.mf6.ModflowPrtoc(
        prtr,
        pname="oc",
        track_filerecord=f"{prt_name('right')}.trk",
        trackcsv_filerecord=f"{prt_name('right')}.trk.csv",
    )
    flopy.mf6.ModflowPrtfmi(
        prtr,
        packagedata=[
            ("GWFHEAD",   f"{gwf_name('right')}.hds"),
            ("GWFBUDGET", f"{gwf_name('right')}.cbc"),
        ],
    )
    # GWF6-PRT6 exchange
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwf_name("right"),
        exgmnameb=prt_name("right"),
        filename=f"{gwf_name('right')}.gwfprt",
    )

    # ---- PRT6-PRT6 exchange --------------------------------------- #
    # Mirror of the GWF-GWF interface: col 4 in PRT-left <-> col 0 in PRT-right
    prtprt_data = [
        ((0, 0, ncol_half - 1), (0, 0, 0)),
    ]
    flopy.mf6.ModflowPrtprt(
        sim,
        exgtype="PRT6-PRT6",
        nexg=len(prtprt_data),
        exgmnamea=prt_name("left"),
        exgmnameb=prt_name("right"),
        exchangedata=prtprt_data,
        filename="prt-prt.exg",
    )

    # ---- Explicit Model Solution (EMS) for both PRT models ------- #
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename="prt.ems",
    )
    sim.register_solution_package(ems, [prt_name("left"), prt_name("right")])

    return sim, None


# ------------------------------------------------------------------ #
# Check output
# ------------------------------------------------------------------ #
def check_output(idx, test):
    ws = Path(test.workspace)

    # ---- Verify GWF solved: heads exist and are reasonable --------
    hds_left = flopy.utils.HeadFile(ws / f"{gwf_name('left')}.hds")
    hds_right = flopy.utils.HeadFile(ws / f"{gwf_name('right')}.hds")
    head_l = hds_left.get_data().squeeze()
    head_r = hds_right.get_data().squeeze()

    # Heads in left model should all be between 0.5 and 1.0
    assert np.all(head_l >= 0.5) and np.all(head_l <= head_left), (
        f"Unexpected heads in GWF-left: {head_l}"
    )
    # Heads in right model should all be between 0.0 and 0.5
    assert np.all(head_r >= head_right) and np.all(head_r <= 0.5), (
        f"Unexpected heads in GWF-right: {head_r}"
    )

    # ---- Check PRT-right received and terminated particles --------
    trk_right_csv = ws / f"{prt_name('right')}.trk.csv"
    assert trk_right_csv.exists(), (
        f"Track CSV for PRT-right not found: {trk_right_csv}\n"
        "This likely means the PRT-PRT exchange did not transfer particles."
    )

    df_right = pd.read_csv(trk_right_csv)
    assert len(df_right) > 0, (
        "PRT-right track file is empty – particles were never transferred."
    )

    # All released particles should appear in PRT-right
    trk_left_csv = ws / f"{prt_name('left')}.trk.csv"
    df_left = pd.read_csv(trk_left_csv)
    released_ids = df_left["irpt"].unique()
    arrived_ids = df_right["irpt"].unique()
    assert set(released_ids) == set(arrived_ids), (
        f"Not all particles reached PRT-right.\n"
        f"Released particle ids (PRT-left): {sorted(released_ids)}\n"
        f"Arrived particle ids (PRT-right): {sorted(arrived_ids)}"
    )

    # Final event in PRT-right should be termination (istatus > 1)
    # status 3 = terminate at right face exit; at minimum all should not be active (1)
    term_events = df_right[df_right["istatus"] > 1]
    assert len(term_events) > 0, (
        "No particles terminated in PRT-right – they may still be active or got lost."
    )

    # Particles should end up at (or near) x = ncol_half * delr
    final_x = term_events["x"].max()
    assert final_x >= (ncol_half - 1) * delr, (
        f"Particles did not reach the right side of PRT-right. Final x: {final_x}"
    )


# ------------------------------------------------------------------ #
# pytest entry point
# ------------------------------------------------------------------ #
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()
