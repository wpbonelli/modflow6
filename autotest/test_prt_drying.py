"""
Tests particle tracking through dry/rewet cells.

The grid is a minimal cross-section with 1 row,
2 columns and 2 layers, with particles released
from the center of the upper left cell, and the
lower left cell inactive. The lower right cell
is a constant head boundary.
 
    z
    |
 (0,2)_____________(2,2)
    |       |        |
    |   * release    |
    |_______|________|
    |       |        |
    |       |        |
    |_______|________|____x
 (0,0)             (2,0) 

"""

from pathlib import Path

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils import PathlineFile
from flopy.utils.binaryfile import HeadFile
from framework import TestFramework
from prt_test_utils import (
    all_equal,
    check_budget_data,
    check_track_data,
    get_model_name,
    get_partdata,
    has_default_boundnames,
    DEFAULT_EXIT_SOLVE_TOL,
)

simname = "prtdry"
cases = [simname, f"{simname}newt"]


nlay = 2
nrow = 1
ncol = 2
tdis_pd = [(1.0, 1, 1.0)]
top = 2
botm = [
    [
        1,
        1,
    ],
    [0, 0],
]
idomain = np.ones((nlay, nrow, ncol), dtype=int)
idomain[1, 0, 0] = 0


def build_gwf_sim(name, ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    # create tdis package
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=len(tdis_pd),
        perioddata=tdis_pd,
    )

    # create gwf model
    gwfname = f"{name}_gwf"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        newtonoptions="NEWTON" if "newt" in name else None,
    )

    # create gwf discretization

    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwf,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        top=top,
        botm=botm,
        idomain=idomain,
    )

    # create gwf initial conditions package
    flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic")

    # create gwf node property flow package
    flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
        gwf,
        pname="npf",
        save_saturation=True,
        save_specific_discharge=True,
    )

    # create wel package
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={0: [[(1, 0, 1), -2.0, -2]]},
        auxiliary=["IFLOWFACE"],
    )

    # create gwf chd package
    chd_pd = {0: [[(0, 0, 0), 1.2], [(0, 0, 1), 0.8]]}
    if False:
        chd_pd.append([(1, 0, 1), 0.2])
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        pname="CHD-1",
        stress_period_data=chd_pd,
    )

    # create gwf output control package
    # output file names
    gwf_budget_file = f"{gwfname}.bud"
    gwf_head_file = f"{gwfname}.hds"
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=gwf_budget_file,
        head_filerecord=gwf_head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # create iterative model solution for gwf model
    ims = flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="BICGSTAB" if "newt" in name else "CG",
    )

    return sim


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # create tdis package
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=len(tdis_pd),
        perioddata=tdis_pd,
    )

    # create prt model
    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name, save_flows=True)

    # create prt discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        top=top,
        botm=botm,
        idomain=idomain,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=0.1)

    # create prp package
    releasepts = [(0, 0, 0, 0, 0.5, 0.5, 0.5)]
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        exit_solve_tolerance=DEFAULT_EXIT_SOLVE_TOL,
        extend_tracking=True,
        local_z=True,
    )

    # create output control package
    prt_budget_file = f"{prt_name}.bud"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[prt_budget_file],
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
        saverecord=[("BUDGET", "ALL")],
    )

    # create the flow model interface
    gwf_name = get_model_name(name, "gwf")
    gwf_budget_file = gwf_ws / f"{gwf_name}.bud"
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )

    # add explicit model solution
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_models(idx, test):
    gwf_ws = test.workspace / "gwf"
    prt_ws = test.workspace / "prt"
    gwf_ws.mkdir()
    prt_ws.mkdir()
    gwf_sim = build_gwf_sim(test.name, gwf_ws, test.targets["mf6"])
    prt_sim = build_prt_sim(test.name, gwf_ws, prt_ws, test.targets["mf6"])
    return gwf_sim, prt_sim


def check_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    gwf_sim, prt_sim = test.sims[:2]
    gwf = gwf_sim.get_model()
    prt = prt_sim.get_model()
    hds = gwf.output.head().get_data()
    mg = gwf.modelgrid
    pathlines = pd.read_csv(prt_ws / f"{prt.name}.trk.csv")

    # setup plot
    plot_data = True
    if plot_data:
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
        for a in ax:
            a.set_aspect("equal")

        # plot pathline in map view
        pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0])
        pmv.plot_grid()
        pmv.plot_array(hds[0], alpha=0.1, cmap="jet")
        pathlines.plot(
            title="MF6 pathlines",
            kind="line",
            x="x",
            y="y",
            ax=ax[0],
            legend=False,
            color="blue",
        )

        # plot pathline in cross section
        pxs = flopy.plot.PlotCrossSection(
            modelgrid=mg, ax=ax[1], line={"row": 0}
        )
        pxs.plot_grid()
        pa = pxs.plot_array(hds, alpha=0.1, cmap="jet")
        plt.colorbar(pa, shrink=0.25)
        pathlines.plot(
            title="MF6 pathlines",
            kind="line",
            x="x",
            y="z",
            ax=ax[1],
            legend=False,
            color="blue",
        )

        # view/save plot
        plt.show()
        plt.savefig(gwf_ws / f"test_{simname}.png")


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        compare=None,
    )
    test.run()
