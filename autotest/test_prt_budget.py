"""
Tests particle mass budget tracking with a very
simple horizontal steady-state flow system. The
grid is a 1x1x10 horizontal line with 10 columns.
Particles are released from the left-most cell.
Pathlines are compared against a MODPATH 7 model.
"""

from pathlib import Path

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.mf6.utils.postprocessing import get_structured_faceflows
from flopy.utils import PathlineFile
from flopy.utils.binaryfile import HeadFile
from framework import TestFramework
from prt_test_utils import (
    HorizontalCase,
    all_equal,
    check_budget_data,
    check_track_data,
    get_model_name,
    get_partdata,
    has_default_boundnames,
)

simname = "prtbud"
cases = [
    simname,  # default (no budget options)
    f"{simname}_2",  # test budget_boundary
    f"{simname}_3",  # test budget_weaksink  
    f"{simname}_5_9",  # test budget_no_exits (covers both cell and subcell exits)
    f"{simname}_6",  # test budget_stopzone
    f"{simname}_7",  # test budget_inactive
    f"{simname}_10",  # test budget_timeout
    f"{simname}_all",  # test all budget options enabled
]


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    # Determine which budget options to enable based on test case name
    budget_boundary = "boundary" in name or "all" in name
    budget_weaksink = "weaksink" in name or "all" in name
    budget_no_exits = "no_exits" in name or "all" in name  
    budget_stopzone = "stopzone" in name or "all" in name
    budget_inactive = "inactive" in name or "all" in name
    budget_timeout = "timeout" in name or "all" in name
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
        nper=HorizontalCase.nper,
        perioddata=[
            (HorizontalCase.perlen, HorizontalCase.nstp, HorizontalCase.tsmult)
        ],
    )

    # create prt model
    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name, save_flows=True)

    # create prt discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=HorizontalCase.nlay,
        nrow=HorizontalCase.nrow,
        ncol=HorizontalCase.ncol,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=HorizontalCase.porosity)

    # convert mp7 to prt release points and check against expectation
    partdata = get_partdata(prt.modelgrid, HorizontalCase.releasepts_mp7)
    coords = partdata.to_coords(prt.modelgrid)
    releasepts = [(i, 0, 0, 0, c[0], c[1], c[2]) for i, c in enumerate(coords)]

    # create prp package
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
        stop_at_weak_sink=False,
        boundnames=True,
        extend_tracking=True,
    )

    # create output control package with budget options
    prt_budget_file = f"{prt_name}.bud"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    oc_kwargs = {
        "pname": "oc",
        "budget_filerecord": [prt_budget_file],
        "track_filerecord": [prt_track_file],
        "trackcsv_filerecord": [prt_track_csv_file],
        "saverecord": [("BUDGET", "ALL")],
    }
    
    # Add budget options if enabled
    if budget_boundary:
        oc_kwargs["budget_boundary"] = True
    if budget_weaksink:
        oc_kwargs["budget_weaksink"] = True
    if budget_no_exits:
        oc_kwargs["budget_no_exits"] = True
    if budget_stopzone:
        oc_kwargs["budget_stopzone"] = True
    if budget_inactive:
        oc_kwargs["budget_inactive"] = True
    if budget_timeout:
        oc_kwargs["budget_timeout"] = True
        
    flopy.mf6.ModflowPrtoc(prt, **oc_kwargs)

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


def build_mp7_sim(name, ws, mp7, gwf):
    partdata = get_partdata(gwf.modelgrid, HorizontalCase.releasepts_mp7)
    mp7_name = get_model_name(name, "mp7")
    pg = flopy.modpath.ParticleGroup(
        particlegroupname="G1",
        particledata=partdata,
        filename=f"{mp7_name}.sloc",
    )
    mp = flopy.modpath.Modpath7(
        modelname=mp7_name,
        flowmodel=gwf,
        exe_name=mp7,
        model_ws=ws,
        headfilename=f"{gwf.name}.hds",
        budgetfilename=f"{gwf.name}.bud",
    )
    mpbas = flopy.modpath.Modpath7Bas(
        mp,
        porosity=HorizontalCase.porosity,
    )
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="pathline",
        trackingdirection="forward",
        budgetoutputoption="summary",
        stoptimeoption="extend",
        particlegroups=[pg],
    )

    return mp


def build_models(idx, test):
    gwf_sim = HorizontalCase.get_gwf_sim(test.name, test.workspace, test.targets["mf6"])
    gwf = gwf_sim.get_model()
    sto = flopy.mf6.ModflowGwfsto(gwf)
    prt_sim = build_prt_sim(
        test.name, test.workspace, test.workspace / "prt", test.targets["mf6"]
    )
    mp7_sim = build_mp7_sim(
        test.name,
        test.workspace / "mp7",
        test.targets["mp7"],
        gwf_sim.get_model(),
    )
    return gwf_sim, prt_sim, mp7_sim


def check_output(idx, test):
    from flopy.plot.plotutil import to_mp7_pathlines

    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    mp7_ws = test.workspace / "mp7"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    mp7_name = get_model_name(name, "mp7")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid

    # check mf6 output files exist
    gwf_budget_file = f"{gwf_name}.bud"
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    assert (gwf_ws / gwf_budget_file).is_file()
    assert (gwf_ws / gwf_head_file).is_file()
    assert (prt_ws / prt_track_file).is_file()
    assert (prt_ws / prt_track_csv_file).is_file()
    assert (prt_ws / prp_track_file).is_file()
    assert (prt_ws / prp_track_csv_file).is_file()

    # check mp7 output files exist
    assert (mp7_ws / mp7_pathline_file).is_file()

    # make sure pathline df has "name" (boundname) column and default values
    assert "name" in mf6_pls
    assert has_default_boundnames(mf6_pls)

    # make sure all mf6 pathline data have correct model and PRP index (1)
    assert all_equal(mf6_pls["imdl"], 1)
    assert all_equal(mf6_pls["iprp"], 1)
    
    # report termination statistics for debugging
    terminated_particles = mf6_pls[mf6_pls["ireason"] > 1]
    if len(terminated_particles) > 0:
        termination_counts = terminated_particles["ireason"].value_counts()
        print(f"Test {name} termination summary: {dict(termination_counts)}")
    else:
        print(f"Test {name}: No particles terminated")

    # check budget data
    check_budget_data(
        prt_ws / f"{name}_prt.lst", HorizontalCase.perlen, HorizontalCase.nper
    )

    # check particle mass budget behavior based on configuration
    check_mass_budget_behavior(name, "PRP", prt_ws, prt_name, mf6_pls)

    # check cell-by-cell flows
    prt_budget_file = prt_ws / f"{prt_name}.bud"
    prt_bud = flopy.utils.CellBudgetFile(prt_budget_file, precision="double")
    prt_bud_data = prt_bud.get_data(kstpkper=(0, 0))
    assert len(prt_bud_data) == 2
    flowja = prt_bud.get_data(text="FLOW-JA-FACE")[0][0, 0, :]
    prp = prt_bud.get_data(text="PRP")[0].squeeze()
    assert flowja.shape == (28,)
    assert prp.shape == (9,)
    frf, fff, flf = get_structured_faceflows(
        flowja,
        grb_file=gwf_ws / f"{gwf_name}.dis.grb",
        verbose=True,
    )
    assert not fff.any()
    assert not flf.any()
    assert frf.any()
    assert all(v == 9 for v in frf[:-1])

    # check mf6 prt particle track data were written to binary/CSV files
    # and that different formats are equal
    for track_csv in [prt_ws / prt_track_csv_file, prt_ws / prp_track_csv_file]:
        check_track_data(
            track_bin=prt_ws / prt_track_file,
            track_hdr=prt_ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
            track_csv=track_csv,
        )

    # convert mf6 pathlines to mp7 format
    mf6_pls = to_mp7_pathlines(mf6_pls)

    # sort both dataframes by particleid and time
    mf6_pls.sort_values(by=["particleid", "time"], inplace=True)
    mp7_pls.sort_values(by=["particleid", "time"], inplace=True)

    # drop columns for which there is no direct correspondence between mf6 and mp7
    del mf6_pls["sequencenumber"]
    del mf6_pls["particleidloc"]
    del mf6_pls["xloc"]
    del mf6_pls["yloc"]
    del mf6_pls["zloc"]
    del mp7_pls["sequencenumber"]
    del mp7_pls["particleidloc"]
    del mp7_pls["xloc"]
    del mp7_pls["yloc"]
    del mp7_pls["zloc"]

    # compare mf6 / mp7 pathline data
    assert mf6_pls.shape == mp7_pls.shape
    assert np.allclose(mf6_pls, mp7_pls, atol=1e-3)


def check_mass_budget_behavior(name, text, prt_ws, prt_name, mf6_pls):
    """Check that particle mass is correctly included/excluded from budget based on termination status."""
    import flopy.utils
    
    # Load budget file
    prt_budget_file = prt_ws / f"{prt_name}.bud"
    prt_bud = flopy.utils.CellBudgetFile(prt_budget_file, precision="double")
    
    # Get budget terms
    try:
        bud_data = prt_bud.get_data(text=text)
        if bud_data:
            bud = bud_data[0].squeeze()
            total_mass_in_budget = bud.sum() if bud.size > 0 else 0.0
        else:
            total_mass_in_budget = 0.0
    except:
        total_mass_in_budget = 0.0
    
    # Count total particles released
    total_particles = len(mf6_pls["irpt"].unique())
    
    # Group particles by termination reason (ireason maps to status codes)
    terminated_particles = mf6_pls[mf6_pls["ireason"] > 1]  # ireason > 1 means terminated
    termination_counts = terminated_particles["ireason"].value_counts()
    
    # Expected mass in budget based on termination configuration
    expected_mass = 0.0
    
    if "boundary" in name and 2 in termination_counts:
        expected_mass += termination_counts[2]  # TERM_BOUNDARY
    if "weaksink" in name and 3 in termination_counts:
        expected_mass += termination_counts[3]  # TERM_WEAKSINK
    if "no_exits" in name:
        if 5 in termination_counts:
            expected_mass += termination_counts[5]  # TERM_NO_EXITS
        if 9 in termination_counts:
            expected_mass += termination_counts[9]  # TERM_NO_EXITS_SUB
    if "stopzone" in name and 6 in termination_counts:
        expected_mass += termination_counts[6]  # TERM_STOPZONE
    if "inactive" in name and 7 in termination_counts:
        expected_mass += termination_counts[7]  # TERM_INACTIVE
    if "timeout" in name and 10 in termination_counts:
        expected_mass += termination_counts[10]  # TERM_TIMEOUT
    if "all" in name:
        # All terminated particles should have mass in budget
        expected_mass = len(terminated_particles)
    
    # For default case (no budget options), no mass should remain in budget
    if name == simname:
        expected_mass = 0.0
    
    # Check that actual mass in budget matches expected
    assert abs(total_mass_in_budget - expected_mass) < 1e-6, \
        f"Test {name}: Expected {expected_mass} mass in budget, got {total_mass_in_budget}. "\
        f"Termination counts: {dict(termination_counts)}"


def plot_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    mp7_ws = test.workspace / "mp7"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    mp7_name = get_model_name(name, "mp7")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid

    # check mf6 output files exist
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # set up plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    for a in ax:
        a.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0])
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mf6_plines = mf6_pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        pl.plot(
            title="MF6 pathlines",
            kind="line",
            x="x",
            y="y",
            ax=ax[0],
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
        )

    # plot mp7 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[1])
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mp7_plines = mp7_pls.groupby(["particleid"])
    for ipl, (pid, pl) in enumerate(mp7_plines):
        pl.plot(
            title="MP7 pathlines",
            kind="line",
            x="x",
            y="y",
            ax=ax[1],
            legend=False,
            color=cm.plasma(ipl / len(mp7_plines)),
        )

    # view/save plot
    plt.show()
    plt.savefig(prt_ws / f"test_{simname}.png")


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()
