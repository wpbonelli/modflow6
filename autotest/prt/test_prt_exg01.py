"""
Test GWF and PRT models in the same simulation
with an exchange.

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme.
Particles are released from the top left cell.

Results are compared against a MODPATH 7 model.

This test includes two cases, one which gives
boundnames to particles and one which does not.
"""


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
from prt_test_utils import check_budget_data, check_track_data, get_gwf_sim

from framework import TestFramework

simname = "prtexg01"
ex = [simname, f"{simname}bnms"]


# model names
def get_model_name(idx, mdl):
    return f"{ex[idx]}_{mdl}"


def build_sim(idx, ws, mf6):
    # create simulation
    name = ex[idx]
    sim, ctx = get_gwf_sim(name, ws, mf6)

    # create prt model
    prtname = get_model_name(idx, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)

    # create prt discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=ctx.nlay,
        nrow=ctx.nrow,
        ncol=ctx.ncol,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=ctx.porosity)

    # create prp package
    rpts = (
        [r + [str(r[0] + 1)] for r in ctx.releasepts_prt]
        if "bnms" in name
        else ctx.releasepts_prt
    )
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(rpts),
        packagedata=rpts,
        perioddata={0: ["FIRST"]},
        boundnames="bnms" in name,
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

    # create a flow model interface
    # todo Mike Fienen's report (crash when FMI created but not needed)
    # flopy.mf6.ModflowPrtfmi(
    #     prt,
    #     packagedata=[
    #         ("GWFHEAD", gwf_head_file),
    #         ("GWFBUDGET", gwf_budget_file),
    #     ],
    # )

    # create exchange
    gwfname = get_model_name(idx, "gwf")
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwfname,
        exgmnameb=prtname,
        filename=f"{gwfname}.gwfprt",
    )

    # add explicit model solution
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prtname}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim, ctx


def build_mp7_sim(ctx, idx, ws, mp7, gwf):
    partdata = flopy.modpath.ParticleData(
        partlocs=[p[0] for p in ctx.releasepts_mp7],
        localx=[p[1] for p in ctx.releasepts_mp7],
        localy=[p[2] for p in ctx.releasepts_mp7],
        localz=[p[3] for p in ctx.releasepts_mp7],
        timeoffset=0,
        drape=0,
    )
    mp7name = get_model_name(idx, "mp7")
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
    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        simulationtype="pathline",
        trackingdirection="forward",
        budgetoutputoption="summary",
        stoptimeoption="extend",
        particlegroups=[pg],
    )

    return mp


def eval_results(ctx, test):
    print(f"Evaluating results for sim {test.name}")
    simpath = Path(test.workspace)

    # check budget data
    check_budget_data(simpath / f"{test.name}_prt.lst", ctx.perlen, ctx.nper)

    # check particle track data
    prt_track_file = simpath / f"{test.name}_prt.trk"
    prt_track_hdr_file = simpath / f"{test.name}_prt.trk.hdr"
    prt_track_csv_file = simpath / f"{test.name}_prt.trk.csv"
    assert prt_track_file.exists()
    assert prt_track_hdr_file.exists()
    assert prt_track_csv_file.exists()
    check_track_data(
        track_bin=prt_track_file,
        track_hdr=prt_track_hdr_file,
        track_csv=prt_track_csv_file,
    )


@pytest.mark.parametrize("idx, name", enumerate(ex))
def test_mf6model(idx, name, function_tmpdir, targets):
    ws = function_tmpdir
    sim, ctx = build_sim(idx, str(ws), targets.mf6)
    sim.write_simulation()

    test = TestFramework(
        name=name,
        workspace=ws,
        targets=targets,
        check=lambda s: eval_results(ctx, s),
        compare=None,
    )
    test.run()

    # model names
    gwfname = get_model_name(idx, "gwf")
    prtname = get_model_name(idx, "prt")
    mp7name = get_model_name(idx, "mp7")

    # extract model objects
    gwf = sim.get_model(gwfname)
    prt = sim.get_model(prtname)

    # extract model grid
    mg = gwf.modelgrid

    # build mp7 model
    mp7sim = build_mp7_sim(ctx, idx, ws, targets.mp7, gwf)

    # run mp7 model
    mp7sim.write_input()
    success, buff = mp7sim.run_model(report=True)
    assert success, pformat(buff)

    # check mf6 output files exist
    gwf_budget_file = f"{gwfname}.bud"
    gwf_head_file = f"{gwfname}.hds"
    prt_track_file = f"{prtname}.trk"
    prt_track_csv_file = f"{prtname}.trk.csv"
    mp7_pathline_file = f"{mp7name}.mppth"
    assert (ws / gwf_budget_file).is_file()
    assert (ws / gwf_head_file).is_file()
    assert (ws / prt_track_file).is_file()
    assert (ws / prt_track_csv_file).is_file()

    # check mp7 output files exist
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
    mf6_pls = pd.read_csv(ws / prt_track_csv_file).replace(
        r"^\s*$", np.nan, regex=True
    )

    # make sure pathline dataframe has "name" column
    assert "name" in mf6_pls

    # check boundname values
    if "bnms" in name:
        # boundnames should be release point numbers (so pandas parses them as ints)
        assert np.array_equal(
            mf6_pls["name"].to_numpy(), mf6_pls["irpt"].to_numpy()
        )
    else:
        # no boundnames given so check for defaults
        assert pd.isna(mf6_pls["name"]).all()

    # check budget data were written to mf6 prt list file
    check_budget_data(ws / f"{name}_prt.lst", ctx.perlen, ctx.nper)

    # check mf6 prt particle track data were written to binary/CSV files
    check_track_data(
        track_bin=ws / prt_track_file,
        track_hdr=ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
        track_csv=ws / prt_track_csv_file,
    )

    # extract head, budget, and specific discharge results from GWF model
    gwf = sim.get_model(gwfname)
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
    # plt.show()
    plt.savefig(ws / f"test_{name}.png")

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
