"""
This test is similar to test_prt_fmi01.py, except
particles are split across two release packages,
and the grid has an inactive region. This tests
that cell numbers recorded in pathline data have
been converted from reduced to user node numbers.
This is verified by using FloPy to intersect path
points with the grid, then compute node numbers.

GWF and PRT models run in separate simulations
via flow model interface.

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme,
except for 2 inactive cells in the bottom left
and top right corners.

Particles are released from the top left cell.

Results are compared against a MODPATH 7 model.

This test case also configures recording events
to check that they can be explicitly specified.
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

from prt_test_utils import check_budget_data, check_track_data, get_event, to_mp7_format


# simulation/model names
simname = "prtfmi02"
gwfname = f"{simname}_gwf"
prtname = f"{simname}_prt"
mp7name = f"{simname}_mp7"

# test cases
ex = [
    f"{simname}all",
    f"{simname}rel",
    f"{simname}trst",
    f"{simname}tstp",
    f"{simname}wksk",
]

# output file names
gwf_budget_file = f"{gwfname}.bud"
gwf_head_file = f"{gwfname}.hds"
prt_track_file = f"{prtname}.trk"
prt_track_csv_file = f"{prtname}.trk.csv"
mp7_pathline_file = f"{mp7name}.mppth"

# model info
nlay = 1
nrow = 10
ncol = 10
top = 1.0
botm = [0.0]
nper = 1
perlen = 1.0
nstp = 1
tsmult = 1.0
porosity = 0.1

# release points
# todo: define for mp7 first, then use flopy utils to convert to global coords for mf6 prt
releasepts_a = [
    # index, k, i, j, x, y, z
    # (0-based indexing converted to 1-based for mf6 by flopy)
    [i, 0, 0, 0, float(f"0.{i + 1}"), float(f"9.{i + 1}"), 0.5]
    for i in range(4)
]
releasepts_b = [
    # index, k, i, j, x, y, z
    # (0-based indexing converted to 1-based for mf6 by flopy)
    [i, 0, 0, 0, float(f"0.{i + 5}"), float(f"9.{i + 5}"), 0.5]
    for i in range(5)
]
releasepts_mp7_a = [
    # node number, localx, localy, localz
    # (0-based indexing converted to 1-based for mf6 by flopy)
    (0, float(f"0.{i + 1}"), float(f"0.{i + 1}"), 0.5)
    for i in range(4)
]
releasepts_mp7_b = [
    # node number, localx, localy, localz
    # (0-based indexing converted to 1-based for mf6 by flopy)
    (0, float(f"0.{i + 5}"), float(f"0.{i + 5}"), 0.5)
    for i in range(5)
]

# idomain
idomain = np.ones((nlay, nrow, ncol), dtype=int)
idomain[0, 0, 9] = 0
idomain[0, 9, 0] = 0
# idomain = idomain.ravel()


def build_gwf_sim(idx, ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=ex[idx],
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    # create tdis package
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=nper,
        perioddata=[(perlen, nstp, tsmult)],
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

    # create gwf discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwf,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
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

    # create gwf chd package
    spd = {
        0: [[(0, 0, 0), 1.0, 1.0], [(0, 9, 9), 0.0, 0.0]],
        1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        pname="CHD-1",
        stress_period_data=spd,
        auxiliary=["concentration"],
    )

    # create gwf output control package
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=gwf_budget_file,
        head_filerecord=gwf_head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # create iterative model solution for gwf model
    ims = flopy.mf6.ModflowIms(sim)

    return sim


def build_prt_sim(idx, ws, mf6):
    # create simulation
    name = ex[idx]
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
        nper=nper,
        perioddata=[(perlen, nstp, tsmult)],
    )

    # create prt model
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)

    # create prt discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        idomain=idomain,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # create prp packages
    event = get_event(name)
    rpts_a = [r + [event] for r in releasepts_a]
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp_a",
        filename=f"{prtname}_a.prp",
        nreleasepts=len(rpts_a),
        packagedata=rpts_a,
        perioddata={0: ["FIRST"]},
    )
    rpts_b = [r + [event] for r in releasepts_b]
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp_b",
        filename=f"{prtname}_b.prp",
        nreleasepts=len(rpts_b),
        packagedata=rpts_b,
        perioddata={0: ["FIRST"]},
    )

    # create output control package
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )

    # create the flow model interface
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


def build_mp7_sim(ws, mp7, gwf):
    pd_a = flopy.modpath.ParticleData(
        partlocs=[p[0] for p in releasepts_mp7_a],
        localx=[p[1] for p in releasepts_mp7_a],
        localy=[p[2] for p in releasepts_mp7_a],
        localz=[p[3] for p in releasepts_mp7_a],
        timeoffset=0,
        drape=0,
    )
    pd_b = flopy.modpath.ParticleData(
        partlocs=[p[0] for p in releasepts_mp7_b],
        localx=[p[1] for p in releasepts_mp7_b],
        localy=[p[2] for p in releasepts_mp7_b],
        localz=[p[3] for p in releasepts_mp7_b],
        timeoffset=0,
        drape=0,
    )
    pg_a = flopy.modpath.ParticleGroup(
        particlegroupname="GA",
        particledata=pd_a,
        filename=f"{mp7name}_a.sloc",
    )
    pg_b = flopy.modpath.ParticleGroup(
        particlegroupname="GB",
        particledata=pd_b,
        filename=f"{mp7name}_b.sloc",
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
        stoptimeoption="extend",
        particlegroups=[pg_a, pg_b],
    )

    return mp


@pytest.mark.parametrize("idx, name", enumerate(ex))
def test_prt_fmi02(idx, name, function_tmpdir, targets):
    ws = function_tmpdir

    # build mf6 simulations
    gwfsim = build_gwf_sim(idx, ws, targets.mf6)
    prtsim = build_prt_sim(idx, ws, targets.mf6)

    # write mf6 simulation input files
    gwfsim.write_simulation()
    
    # run mf6 gwf simulation
    success, _ = gwfsim.run_simulation()
    assert success

    # write mf6 prt simulation input files
    prtsim.write_simulation()

    # run mf6 prt simulation
    success, _ = prtsim.run_simulation()
    assert success

    # extract models
    gwf = gwfsim.get_model(gwfname)
    prt = prtsim.get_model(prtname)

    # extract model grid
    mg = gwf.modelgrid

    # build mp7 model
    mp7sim = build_mp7_sim(ws, targets.mp7, gwf)

    # run mp7 model
    mp7sim.write_input()
    success, _ = mp7sim.run_model()
    assert success

    # check mf6 output files exist
    assert (ws / gwf_budget_file).is_file()
    assert (ws / gwf_head_file).is_file()
    assert (ws / prt_track_file).is_file()
    assert (ws / prt_track_csv_file).is_file()

    # check mp7 output files exist
    assert (ws / mp7_pathline_file).is_file()

    # load mp7 pathline results
    plf = PathlineFile(ws / mp7_pathline_file)
    mp7_pldata = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based
    mp7_pldata["particleid"] = mp7_pldata["particleid"] + 1
    mp7_pldata["particlegroup"] = mp7_pldata["particlegroup"] + 1
    mp7_pldata["node"] = mp7_pldata["node"] + 1
    mp7_pldata["k"] = mp7_pldata["k"] + 1

    # load mf6 pathline results
    mf6_pldata = pd.read_csv(ws / prt_track_csv_file)

    # if event is ALL, output should be the same as MODPATH 7,
    # so continue with comparisons.
    # if event is RELEASE, expect 1 location for each particle.
    # if event is TRANSIT, expect full results minus start loc.
    # if event is TIMESTEP or WEAKSINK, output should be empty.
    # in either case, return early and skip MP7 comparison.
    event = get_event(ex[idx])
    if event == "RELEASE":
        assert len(mf6_pldata) == len(releasepts_a) + len(releasepts_b)
        return
    elif event == "TRANSIT":
        assert len(mf6_pldata) == (len(mp7_pldata) - 2 * (len(releasepts_a) + len(releasepts_b)))
        return
    elif event == "TIMESTEP" or event == "WEAKSINK":
        assert len(mf6_pldata) == 0
        return

    # make sure mf6 pathline data have correct
    #   - model index (1)
    #   - PRP index (1 or 2, depending on release point index)
    def all_equal(col, val):
        a = col.to_numpy()
        return a[0] == val and (a[0] == a).all()

    assert all_equal(mf6_pldata["imdl"], 1)
    assert set(mf6_pldata[mf6_pldata["iprp"] == 1]["irpt"].unique()) == set(
        range(1, 5)
    )
    assert set(mf6_pldata[mf6_pldata["iprp"] == 2]["irpt"].unique()) == set(
        range(1, 6)
    )

    # check budget data were written to mf6 prt list file
    check_budget_data(ws / f"{simname}_prt.lst", perlen, nper)

    # check mf6 prt particle track data were written to binary/CSV files
    check_track_data(
        track_bin=ws / prt_track_file,
        track_hdr=ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
        track_csv=ws / prt_track_csv_file,
    )

    # check that particle names are particle indices
    # assert len(mf6_pldata) == len(mf6_pldata[mf6_pldata['irpt'].astype(str).eq(mf6_pldata['name'])])

    # get head, budget, and spdis results from GWF model
    hds = HeadFile(ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # setup plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(13, 13))
    for a in ax:
        a.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0])
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
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
    pmv.plot_array(hds[0], alpha=0.1)
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

    # drop duplicate locations
    # (mp7 includes a duplicate location at the end of each pathline??)
    cols = ["x", "y", "z", "time"]
    mp7_pldata = mp7_pldata.drop_duplicates(subset=cols)
    mf6_pldata_mp7 = mf6_pldata_mp7.drop_duplicates(subset=cols)

    # drop columns for which there is no direct correspondence between mf6 and mp7
    del mf6_pldata_mp7["particleid"]
    del mf6_pldata_mp7["sequencenumber"]
    del mf6_pldata_mp7["particleidloc"]
    del mf6_pldata_mp7["xloc"]
    del mf6_pldata_mp7["yloc"]
    del mf6_pldata_mp7["zloc"]
    del mp7_pldata["particleid"]
    del mp7_pldata["sequencenumber"]
    del mp7_pldata["particleidloc"]
    del mp7_pldata["xloc"]
    del mp7_pldata["yloc"]
    del mp7_pldata["zloc"]

    # sort both dataframes by particleid and time
    mf6_pldata_mp7 = mf6_pldata_mp7.sort_values(by=cols)
    mp7_pldata = mp7_pldata.sort_values(by=cols)

    # compare mf6 / mp7 pathline data
    assert mf6_pldata_mp7.shape == mp7_pldata.shape
    assert np.allclose(mf6_pldata_mp7, mp7_pldata, atol=1e-3)

    # check that cell numbers are correct
    for i, row in list(mf6_pldata.iterrows()):
        x, y, z, t, ilay, icell = (
            row.x,
            row.y,
            row.z,
            row.t,
            row.ilay,
            row.icell,
        )
        k, i, j = mg.intersect(x, y, z)
        nn = mg.get_node([k, i, j]) + 1
        neighbors = mg.neighbors(nn)
        assert np.isclose(nn, icell, atol=1) or any(
            (nn - 1) == n for n in neighbors
        )
        assert ilay == (k + 1) == 1
