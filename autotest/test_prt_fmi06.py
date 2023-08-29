"""
Tests particle tracking on a vertex (DISV) grid.
"""

from pathlib import Path
from typing import Optional

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils import PathlineFile
from flopy.utils.binaryfile import HeadFile

from prt_test_utils import (
    all_equal,
    check_budget_data,
    check_track_data,
    get_partdata,
    has_default_boundnames,
    to_mp7_format,
)
from test_gwf_disv_uzf import create_disv_mesh

# simulation name
simname = "prtfmi06"

# test cases
ex = [
    f"{simname}a",
]

# model info
nlay = 1
nper = 1
perlen = 10
nstp = 5
tsmult = 1.0
tdis_rc = [(perlen, nstp, tsmult)]
botm = [20.0]
strt = 20
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97
porosity = 0.1

# vertex grid properties
verts, cell2d = create_disv_mesh()

# release points in mp7 format (using local coordinates)
releasepts_mp7 = [
    # node number, localx, localy, localz
    # (0-based indexing converted to 1-based for mp7 by flopy)
    (i * 10, 0.5, 0.5, 0.5)
    for i in range(10)
]


def build_gwf_sim(idx, dir, mf6):
    # model name
    gwfname = f"{ex[idx]}_gwf"

    # build MODFLOW 6 files
    ws = dir
    sim = flopy.mf6.MFSimulation(
        sim_name=gwfname, version="mf6", exe_name=mf6, sim_ws=ws
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=tdis_rc
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, newtonoptions="NEWTON", save_flows=True
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="MODERATE",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    ncpl = len(cell2d)
    nvert = len(verts)
    disv = flopy.mf6.ModflowGwfdisv(
        gwf,
        nlay=nlay,
        ncpl=ncpl,
        nvert=nvert,
        top=25.0,
        botm=botm,
        vertices=verts,
        cell2d=cell2d,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
    )

    # constant head boundary
    spd = {
        0: [[(0, 0), 1.0, 1.0], [(0, 99), 0.0, 0.0]],
        # 1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        pname="CHD-1",
        stress_period_data=spd,
        auxiliary=["concentration"],
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord="{}.cbc".format(gwfname),
        head_filerecord="{}.hds".format(gwfname),
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename="{}.oc".format(gwfname),
    )

    # Print human-readable heads
    obs_lst = []
    for k in np.arange(0, 1, 1):
        for i in np.arange(40, 50, 1):
            obs_lst.append(["obs_" + str(i + 1), "head", (k, i)])

    obs_dict = {f"{gwfname}.obs.csv": obs_lst}
    obs = flopy.mf6.ModflowUtlobs(
        gwf, pname="head_obs", digits=20, continuous=obs_dict
    )

    return sim


def build_prt_sim(idx, ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=ex[idx],
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=tdis_rc
    )

    # create prt model
    prtname = f"{ex[idx]}_prt"
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)

    # create prt discretization
    ncpl = len(cell2d)
    nvert = len(verts)
    disv = flopy.mf6.ModflowGwfdisv(
        prt,
        nlay=nlay,
        ncpl=ncpl,
        nvert=nvert,
        top=25.0,
        botm=botm,
        vertices=verts,
        cell2d=cell2d,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # convert mp7 particledata to prt release points
    partdata = get_partdata(prt.modelgrid, releasepts_mp7)
    coords = partdata.to_coords(prt.modelgrid)
    releasepts = [
        (i, (0, r[0]), c[0], c[1], c[2])
        for i, (r, c) in enumerate(zip(releasepts_mp7, coords))
    ]

    # create prp package
    prp_track_file = f"{prtname}.prp.trk"
    prp_track_csv_file = f"{prtname}.prp.trk.csv"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        stop_at_weak_sink="saws" in prtname,
        boundnames=True,
    )

    # create output control package
    prt_track_file = f"{prtname}.trk"
    prt_track_csv_file = f"{prtname}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )

    # create the flow model interface
    gwfname = f"{ex[idx]}_gwf"
    gwf_budget_file = f"{gwfname}.cbc"
    gwf_head_file = f"{gwfname}.hds"
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
        filename=f"{prtname}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_mp7_sim(idx, ws, mp7, gwf):
    # convert mp7 particledata to prt release points
    partdata = get_partdata(gwf.modelgrid, releasepts_mp7)

    # create modpath 7 simulation
    mp7name = f"{ex[idx]}_mp7"
    pg = flopy.modpath.ParticleGroup(
        particlegroupname="G1",
        particledata=partdata,
        filename=f"{mp7name}.sloc",
    )
    mp = flopy.modpath.Modpath7(
        modelname=mp7name,
        flowmodel=gwf,
        exe_name=mp7,
        model_ws=ws,
    )
    mpbas = flopy.modpath.Modpath7Bas(
        mp,
        porosity=porosity,
    )
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="pathline",
        trackingdirection="forward",
        budgetoutputoption="summary",
        stoptimeoption="total",
        particlegroups=[pg],
    )

    return mp


@pytest.mark.parametrize("idx, name", list(enumerate(ex)))
def test_prt_fmi06(idx, name, function_tmpdir, targets):
    # workspace
    ws = function_tmpdir

    # test case name
    name = ex[idx]

    # model names
    gwfname = f"{ex[idx]}_gwf"
    prtname = f"{ex[idx]}_prt"
    mp7name = f"{ex[idx]}_mp7"

    # build mf6 models
    gwfsim = build_gwf_sim(idx, ws, targets.mf6)
    prtsim = build_prt_sim(idx, ws, targets.mf6)

    # run mf6 models
    for sim in [gwfsim, prtsim]:
        sim.write_simulation()
        success, _ = sim.run_simulation()
        assert success

    # extract mf6 models
    gwf = gwfsim.get_model(gwfname)
    prt = prtsim.get_model(prtname)

    # extract model grid
    mg = gwf.modelgrid

    # todo build mp7 model
    mp7sim = build_mp7_sim(idx, ws, targets.mp7, gwf)

    # todo run mp7 model
    mp7sim.write_input()
    success, _ = mp7sim.run_model()
    assert success

    # check mf6 output files exist
    gwf_budget_file = f"{gwfname}.cbc"
    gwf_head_file = f"{gwfname}.hds"
    prt_track_file = f"{prtname}.trk"
    prt_track_csv_file = f"{prtname}.trk.csv"
    prp_track_file = f"{prtname}.prp.trk"
    prp_track_csv_file = f"{prtname}.prp.trk.csv"
    assert (ws / gwf_budget_file).is_file()
    assert (ws / gwf_head_file).is_file()
    assert (ws / prt_track_file).is_file()
    assert (ws / prt_track_csv_file).is_file()
    assert (ws / prp_track_file).is_file()
    assert (ws / prp_track_csv_file).is_file()

    # check mp7 output files exist
    mp7_pathline_file = f"{mp7name}.mppth"
    assert (ws / mp7_pathline_file).is_file()

    # load mp7 pathline results
    plf = PathlineFile(ws / mp7_pathline_file)
    mp7_pldata = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pldata["particleid"] = mp7_pldata["particleid"] + 1
    mp7_pldata["particlegroup"] = mp7_pldata["particlegroup"] + 1
    mp7_pldata["node"] = mp7_pldata["node"] + 1
    mp7_pldata["k"] = mp7_pldata["k"] + 1

    # load mf6 pathline results
    mf6_pldata = pd.read_csv(ws / prt_track_csv_file, na_filter=False)

    # make sure pathline df has "name" (boundname) column and default values
    assert "name" in mf6_pldata
    assert has_default_boundnames(mf6_pldata)

    # make sure all mf6 pathline data have correct model and PRP index (1)
    assert all_equal(mf6_pldata["imdl"], 1)
    assert all_equal(mf6_pldata["iprp"], 1)

    # check budget data were written to mf6 prt list file
    check_budget_data(ws / f"{name}_prt.lst", perlen, nper, nstp)

    # check mf6 prt particle track data were written to binary/CSV files
    # and that different formats are equal
    for track_csv in [ws / prt_track_csv_file, ws / prp_track_csv_file]:
        check_track_data(
            track_bin=ws / prt_track_file,
            track_hdr=ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
            track_csv=track_csv,
        )

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # setup plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    for a in ax:
        a.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0])
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.2)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mf6_plines = mf6_pldata.groupby(["iprp", "irpt", "trelease"])
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
    pmv.plot_array(hds[0], alpha=0.2)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mp7_plines = mp7_pldata.groupby(["particleid"])
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
    # plt.show()
    plt.savefig(ws / f"test_{simname}.png")

    # convert mf6 pathlines to mp7 format
    mf6_pldata_mp7 = to_mp7_format(mf6_pldata)

    # sort both dataframes by particleid and time
    mf6_pldata_mp7.sort_values(by=["particleid", "time"], inplace=True)
    mp7_pldata.sort_values(by=["particleid", "time"], inplace=True)

    # drop columns for which there is no direct correspondence between mf6 and mp7
    del mf6_pldata_mp7["sequencenumber"]
    del mf6_pldata_mp7["particleidloc"]
    del mf6_pldata_mp7["xloc"]
    del mf6_pldata_mp7["yloc"]
    del mf6_pldata_mp7["zloc"]
    del mf6_pldata_mp7["node"]  # node numbers reversed in y direction in mp7
    del mp7_pldata["sequencenumber"]
    del mp7_pldata["particleidloc"]
    del mp7_pldata["xloc"]
    del mp7_pldata["yloc"]
    del mp7_pldata["zloc"]
    del mp7_pldata["node"]

    # compare mf6 / mp7 pathline data
    assert mf6_pldata_mp7.shape == mp7_pldata.shape
    assert np.allclose(mf6_pldata_mp7, mp7_pldata, atol=1e-3)
