"""
Test PRT model with ATS-enabled GWF model.

This test combines a simple GWF model with adaptive time stepping (ATS)
and a PRT model that reads flow results via FMI. The purpose is to
determine if PRT crashes when used with ATS, and to understand what
changes are needed to support ATS in particle tracking.
"""

import os
from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework

simname = "prt_ats"
cases = [simname]

# Grid dimensions
nlay, nrow, ncol = 1, 10, 10
delr = 10.0
delc = 10.0
top = 10.0
botm = [0.0]
strt = 9.0
laytyp = 1
hk = 1.0
sy = 0.1
ss = 0.0
porosity = 0.1

# Time discretization for GWF (ATS will override this)
perlen = [100.0]
nper = len(perlen)
nstp = [1]  # ATS will determine actual number of steps
tsmult = [1.0]

# ATS parameters
dt0 = 10.0
dtmin = 1.0e-3
dtmax = 20.0
dtadj = 2.0
dtfailadj = 5.0


def get_model_name(name, model_type):
    return f"{name}_{model_type}"


def build_gwf_sim(name, ws, mf6):
    """Build GWF simulation with ATS enabled."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    # Create TDIS package with ATS
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    # Add ATS settings
    ats_filerecord = f"{name}.ats"
    atsperiod = [(0, dt0, dtmin, dtmax, dtadj, dtfailadj)]
    tdis.ats.initialize(
        maxats=len(atsperiod),
        perioddata=atsperiod,
        filename=ats_filerecord,
    )

    # Create GWF model
    gwf_name = get_model_name(name, "gwf")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, save_flows=True)

    # Create IMS package
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1e-6,
        outer_maximum=100,
        inner_maximum=100,
        inner_dvclose=1e-6,
        rcloserecord=1e-6,
        linear_acceleration="CG",
        filename=f"{gwf_name}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    # Create DIS package
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # Initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # Node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=laytyp,
        k=hk,
    )

    # Storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        steady_state={0: False},
        transient={0: True},
    )

    # Well in middle of domain
    wel_spd = {0: [[(0, 5, 5), -5.0]]}
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=wel_spd,
        save_flows=True,
    )

    # Constant head on left and right boundaries
    chd_spd = []
    for i in range(nrow):
        chd_spd.append([(0, i, 0), 9.0])  # left side
        chd_spd.append([(0, i, ncol - 1), 8.0])  # right side
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: chd_spd},
        save_flows=True,
    )

    # Output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwf_name}.bud",
        head_filerecord=f"{gwf_name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    """Build PRT simulation to track particles through ATS flow field."""
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # Create TDIS package for PRT
    # NOTE: This uses the same discretization as the GWF model's
    # original TDIS (before ATS adjustments). This is likely the
    # problem - PRT doesn't know about the ATS time steps.
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    # Create PRT model
    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name, save_flows=True)

    # Create DIS package
    flopy.mf6.ModflowGwfdis(
        prt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # Create MIP package (porosity)
    flopy.mf6.ModflowPrtmip(prt, porosity=porosity)

    # Create PRP package (particle release)
    # Release 9 particles from the well cell
    releasepts = []
    for i in range(9):
        x = (5 + 0.5) * delr + (i % 3 - 1) * delr * 0.2
        y = (5 + 0.5) * delc + (i // 3 - 1) * delc * 0.2
        z = top - 0.5
        releasepts.append((i, 0, 5, 5, x, y, z))

    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    flopy.mf6.ModflowPrtprp(
        prt,
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
    )

    # Create output control package
    prt_budget_file = f"{prt_name}.bud"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        budget_filerecord=[prt_budget_file],
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
        saverecord=[("BUDGET", "ALL")],
    )

    # Create FMI package to read GWF results
    gwf_name = get_model_name(name, "gwf")
    gwf_budget_file = gwf_ws / f"{gwf_name}.bud"
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
    grb_file = gwf_ws / f"{gwf_name}.dis.grb"
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFGRID", grb_file),
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )

    # Add explicit model solution
    ems = flopy.mf6.ModflowEms(
        sim,
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_models(idx, test):
    """Build both GWF and PRT simulations."""
    gwf_sim = build_gwf_sim(test.name, test.workspace, test.targets["mf6"])
    prt_sim = build_prt_sim(
        test.name, test.workspace, test.workspace / "prt", test.targets["mf6"]
    )
    return gwf_sim, prt_sim


def check_output(idx, test):
    """Check that both models ran successfully and examine output."""
    from flopy.utils.binaryfile import HeadFile

    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    # Check that GWF produced output
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
    gwf_budget_file = gwf_ws / f"{gwf_name}.bud"
    assert gwf_head_file.is_file(), "GWF head file not found"
    assert gwf_budget_file.is_file(), "GWF budget file not found"

    # Load and check GWF budget to see ATS time steps
    fpth = gwf_ws / f"{name}.lst"
    if fpth.is_file():
        mflist = flopy.utils.Mf6ListBudget(fpth)
        inc = mflist.get_incremental()
        print(f"\nATS created {len(inc['totim'])} time steps")
        print(f"Time steps: {inc['totim']}")

    # Check that PRT produced output (or failed as expected)
    prt_track_csv_file = prt_ws / f"{prt_name}.trk.csv"
    prt_budget_file = prt_ws / f"{prt_name}.bud"

    # If we get here, see if PRT output exists
    if prt_track_csv_file.is_file():
        import pandas as pd

        tracks = pd.read_csv(prt_track_csv_file)
        print(f"\nPRT tracked {len(tracks)} pathline points")
        print(f"Tracking completed successfully!")
    else:
        print("\nPRT track file not found - model may have failed")


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    """Test PRT with ATS-enabled GWF model."""
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
