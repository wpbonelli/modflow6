"""
This test checks that zones defined in the model
input package (MIP) are correctly applied to the
particle tracking simulation, i.e. particles are
terminated when they enter the appropriate zone.

GWF and PRT models run in separate simulations
via flow model interface.

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme.

There are two stop zones in the top right and
bottom left of the grid.

Particles are released from the top left cell
and are either captured by either of the stop
zones, or continue to the bottom right cell.

Results are compared against a MODPATH 7 model.
"""


import os
from pathlib import Path

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from flopy.utils import PathlineFile
from flopy.utils.binaryfile import HeadFile
from matplotlib.collections import LineCollection

from prt_test_utils import check_budget_data, check_track_data, to_mp7_format

# model names
name = "prtfmi02"
gwfname = f"{name}"
prtname = f"{name}_prt"
mp7name = f"{name}_mp7"

# output file names
gwf_budget_file = f"{gwfname}.bud"
gwf_head_file = f"{gwfname}.hds"
prt_track_file = f"{prtname}.trk"
prt_track_csv_file = f"{prtname}.trk.csv"
mp7_pathline_file = f"{mp7name}.mppth"

# problem info
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
releasepts = [
    # particle index, k, i, j, x, y, z
    # (0-based indexing converted to 1-based for mf6 by flopy)
    (i, 0, 0, 0, float(f"0.{i + 1}"), float(f"9.{i + 1}"), 0.5)
    for i in range(9)
]
releasepts_mp7 = [
    # node number, localx, localy, localz
    # (0-based indexing converted to 1-based for mp7 by flopy)
    (0, float(f"0.{i + 1}"), float(f"0.{i + 1}"), 0.5)
    for i in range(9)
]

# stop zones
izone = np.zeros((nlay, nrow, ncol), dtype=int)
izones = [(0, 1, 8), (0, 8, 1)]
for iz in izones:
    izone[iz] = 1


def build_gwf_sim(ws, mf6):
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


def build_prt_sim(ws, mf6):
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
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(
        prt,
        pname="mip",
        porosity=porosity,
        izone=izone,  # todo reinstate when izone implemented
    )

    # create prp package
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        istopzone=1,
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
    partdata = flopy.modpath.ParticleData(
        partlocs=[p[0] for p in releasepts_mp7],
        localx=[p[1] for p in releasepts_mp7],
        localy=[p[2] for p in releasepts_mp7],
        localz=[p[3] for p in releasepts_mp7],
        timeoffset=0,
        drape=0,
    )
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
        stoptimeoption="extend",
        stopzone=1,
        zones=izone,
        zonedataoption="on",
        particlegroups=[pg],
    )

    return mp


def test_prt_fmi03(function_tmpdir, targets):
    ws = function_tmpdir

    # build mf6 simulations
    gwfsim = build_gwf_sim(ws, targets.mf6)
    prtsim = build_prt_sim(ws, targets.mf6)

    # run mf6 simulations
    for sim in [gwfsim, prtsim]:
        sim.write_simulation()
        success, _ = sim.run_simulation()
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

    # check mf6 cell budget file
    check_budget_data(ws / f"{name}_prt.lst", perlen, nper)

    # check mf6 track data
    check_track_data(
        track_bin=ws / prt_track_file,
        track_hdr=ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
        track_csv=ws / prt_track_csv_file,
    )

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

    def sort_square_verts(verts):
        """Sort 4 or more points on a square in clockwise order, starting with the top-left point"""

        # sort by y coordinate
        verts.sort(key=lambda v: v[1], reverse=True)

        # separate top and bottom rows
        y0 = verts[0][1]
        t = [v for v in verts if v[1] == y0]
        b = verts[len(t) :]

        # sort top and bottom rows by x coordinate
        t.sort(key=lambda v: v[0])
        b.sort(key=lambda v: v[0])

        # return vertices in clockwise order
        return t + list(reversed(b))

    def plot_stop_zone(nn, ax):
        ifaces = []
        iverts = mg.iverts[nn]

        # sort vertices of well cell in clockwise order
        verts = [tuple(mg.verts[v]) for v in iverts]
        sorted_verts = sort_square_verts(list(set(verts.copy())))
        for i in range(len(sorted_verts) - 1):
            if i == 0:
                p0 = sorted_verts[-1]
                p1 = sorted_verts[i]
                ifaces.append([p0, p1])
            p0 = sorted_verts[i]
            p1 = sorted_verts[(i + 1)]
            ifaces.append([p0, p1])

        lc = LineCollection(ifaces, color="red", lw=4)
        ax.add_collection(lc)

    # plot stop zones
    for iz in izones:
        for a in ax:
            plot_stop_zone(mg.get_node([iz])[0], a)

    # view/save plot
    # plt.show()
    plt.savefig(ws / f"test_{name}.png")

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
