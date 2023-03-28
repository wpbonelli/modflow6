"""
This test checks that zones defined in the model
input package (MIP) are correctly applied to the
particle tracking simulation, i.e. particles are
terminated when they enter the appropriate zone.

GWF and PRT models run in separate simulations
via flow model interface.

The grid is a 10x10 square, based on the flow
system from the FloPy readme. Two test cases
are defined with 1 and 2 layers respectively
(to test zone data can be read as 2D or 3D).

There are two stop zones in the top right and
bottom left of the grid.

Particles are released from the top left cell
and are either captured by either of the stop
zones, or continue to the bottom right cell.

Results are compared against a MODPATH 7 model.
"""


from itertools import repeat
from pathlib import Path
from pprint import pformat

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.plot.plotutil import to_mp7_pathlines
from flopy.utils import PathlineFile
from flopy.utils.binaryfile import HeadFile
from matplotlib.collections import LineCollection
from prt_test_utils import (check_budget_data, check_track_data, get_gwf_sim,
                            get_model_name)

simname = "prtfmi03"
ex = [f"{simname}_l1", f"{simname}_l2"]
stopzone_cells = [(0, 1, 8), (0, 8, 1)]


def create_izone(nlay, nrow, ncol):
    izone = np.zeros((nlay, nrow, ncol), dtype=int)
    for iz in stopzone_cells:
        izone[iz] = 1
    return izone


def build_prt_sim(ctx, name, ws, mf6):
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
        nper=ctx.nper,
        perioddata=[(ctx.perlen, ctx.nstp, ctx.tsmult)],
    )

    # create prt model
    prtname = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)

    # create prt discretization
    ctx.nlay = int(name[-1])
    botm = [ctx.top - (k + 1) for k in range(ctx.nlay)]
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=ctx.nlay,
        nrow=ctx.nrow,
        ncol=ctx.ncol,
        top=ctx.top,
        botm=botm,
    )

    # create mip package
    izone = create_izone(ctx.nlay, ctx.nrow, ctx.ncol)
    flopy.mf6.ModflowPrtmip(
        prt,
        pname="mip",
        porosity=ctx.porosity,
        izone=izone,  # if ctx.nlay == 1 else np.array([izone, izone]),
    )

    # create prp package
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(ctx.releasepts_prt),
        packagedata=ctx.releasepts_prt,
        perioddata={0: ["FIRST"]},
        istopzone=1,
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
    gwfname = get_model_name(name, "gwf")
    gwf_budget_file = f"{gwfname}.bud"
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

    return sim, ctx


def build_mp7_sim(ctx, name, ws, mp7, gwf):
    partdata = flopy.modpath.ParticleData(
        partlocs=[p[0] for p in ctx.releasepts_mp7],
        localx=[p[1] for p in ctx.releasepts_mp7],
        localy=[p[2] for p in ctx.releasepts_mp7],
        localz=[p[3] for p in ctx.releasepts_mp7],
        timeoffset=0,
        drape=0,
    )
    mp7name = get_model_name(name, "mp7")
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
        porosity=ctx.porosity,
    )
    ctx.nlay = int(name[-1])
    izone = create_izone(ctx.nlay, ctx.nrow, ctx.ncol)
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="pathline",
        trackingdirection="forward",
        budgetoutputoption="summary",
        stoptimeoption="extend",
        stopzone=1,
        zones=izone,  # if ctx.nlay == 1 else np.array([izone, izone]),
        zonedataoption="on",
        particlegroups=[pg],
    )

    return mp


@pytest.mark.parametrize("name", ex)
def test_mf6model(name, function_tmpdir, targets):
    ws = function_tmpdir

    # define model names
    gwfname = get_model_name(name, "gwf")
    prtname = get_model_name(name, "prt")
    mp7name = get_model_name(name, "mp7")

    # build mf6 simulations
    gwfsim, ctx = get_gwf_sim(name, ws, targets.mf6)
    gwf = gwfsim.get_model()
    dis = gwf.get_package("DIS")
    ctx.nlay = int(name[-1])
    botm = [ctx.top - (k + 1) for k in range(ctx.nlay)]
    botm_data = np.array(
        [list(repeat(b, ctx.nrow * ctx.ncol)) for b in botm]
    ).reshape((ctx.nlay, ctx.nrow, ctx.ncol))
    dis.nlay = ctx.nlay
    dis.botm.set_data(botm_data)
    prtsim, ctx = build_prt_sim(ctx, name, ws, targets.mf6)

    # run mf6 simulations
    for sim in [gwfsim, prtsim]:
        sim.write_simulation()
        success, buff = sim.run_simulation(report=True)
        assert success, pformat(buff)

    # extract mf6 models
    gwf = gwfsim.get_model(gwfname)
    prt = prtsim.get_model(prtname)

    # extract model grid
    mg = gwf.modelgrid

    # build mp7 model
    mp7sim = build_mp7_sim(ctx, name, ws, targets.mp7, gwf)

    # run mp7 model
    mp7sim.write_input()
    success, buff = mp7sim.run_model(report=True)
    assert success, pformat(buff)

    # check mf6 output files exist
    gwf_budget_file = f"{gwfname}.bud"
    gwf_head_file = f"{gwfname}.hds"
    prt_track_file = f"{prtname}.trk"
    prt_track_csv_file = f"{prtname}.trk.csv"
    assert (ws / gwf_budget_file).is_file()
    assert (ws / gwf_head_file).is_file()
    assert (ws / prt_track_file).is_file()
    assert (ws / prt_track_csv_file).is_file()

    # check mp7 output files exist
    mp7_pathline_file = f"{mp7name}.mppth"
    assert (ws / mp7_pathline_file).is_file()

    # load mp7 pathline results
    plf = PathlineFile(ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # load mf6 pathline results
    mf6_pls = pd.read_csv(ws / prt_track_csv_file)

    # check budget data were written to mf6 prt list file
    check_budget_data(ws / f"{name}_prt.lst", ctx.perlen, ctx.nper)

    # check mf6 prt particle track data were written to binary/CSV files
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

    # setup map view plot
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
    for iz in stopzone_cells:
        for a in ax:
            plot_stop_zone(mg.get_node([iz])[0], a)

    # view/save plot
    # plt.show()
    plt.savefig(ws / f"test_{name}_map.png")

    # check that cell numbers are correct
    for i, row in list(mf6_pls.iterrows()):
        # todo debug final cell number disagreement
        if row.ireason == 3:  # termination
            continue

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

    # convert mf6 pathlines to mp7 format
    mf6_pls = to_mp7_pathlines(mf6_pls)

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

    # drop duplicates and sort both dataframes
    # todo debug why necessary to drop dupes
    cols = ["particleid", "time"]
    mp7_pls = mp7_pls.drop_duplicates(subset=cols)
    mf6_pls = mf6_pls.drop_duplicates(subset=cols)
    mf6_pls = mf6_pls.sort_values(by=cols)
    mp7_pls = mp7_pls.sort_values(by=cols)

    # compare mf6 / mp7 pathline data
    assert mf6_pls.shape == mp7_pls.shape
    assert np.allclose(mf6_pls, mp7_pls, atol=1e-3)