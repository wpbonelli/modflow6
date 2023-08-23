"""
Test GWF and PRT models in the same simulation
with an exchange.

The grid is a 10x10 square with a single layer.
The same flow system shown on the FloPy readme.
Particles are released from the top left cell.

Results are compared against a MODPATH 7 model.
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
from prt_test_utils import check_budget_data, check_track_data, to_mp7_format
from simulation import TestSimulation


# simulation/model names
simname = "prtexg1"
gwfname = f"{simname}_gwf"
prtname = f"{simname}_prt"
mp7name = f"{simname}_mp7"

# test cases
ex = [simname]

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
releasepts = [
    # index, k, i, j, x, y, z
    # (0-based indexing converted to 1-based for mf6 by flopy)
    (i, 0, 0, 0, float(f"0.{i + 1}"), float(f"9.{i + 1}"), 0.5, "ALL")
    for i in range(9)
]
releasepts_mp7 = [
    # node number, localx, localy, localz
    # (0-based indexing converted to 1-based for mp7 by flopy)
    (0, float(f"0.{i + 1}"), float(f"0.{i + 1}"), 0.5)
    for i in range(9)
]


def build_sim(idx, ws, mf6):
    from test_prt_fmi01 import build_gwf_sim

    # create simulation
    sim = build_gwf_sim(ex[idx], ws, mf6)

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
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    # create prp package
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
    )

    # create output control package
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )

    # create a flow model interface
    # todo Fienen's report (crash when FMI created but not needed)
    # flopy.mf6.ModflowPrtfmi(
    #     prt,
    #     packagedata=[
    #         ("GWFHEAD", gwf_head_file),
    #         ("GWFBUDGET", gwf_budget_file),
    #     ],
    # )

    # create exchange
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
        particlegroups=[pg],
    )

    return mp


def eval_results(sim):
    print(f"Evaluating results for sim {sim.name}")
    simpath = Path(sim.simpath)

    # check budget data
    check_budget_data(simpath / f"{sim.name}_prt.lst", perlen, nper)

    # check particle track data
    prt_track_file = simpath / f"{sim.name}_prt.trk"
    prt_track_hdr_file = simpath / f"{sim.name}_prt.trk.hdr"
    prt_track_csv_file = simpath / f"{sim.name}_prt.trk.csv"
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
    sim = build_sim(idx, str(ws), targets.mf6)
    sim.write_simulation()

    test = TestFramework()
    test.run(
        TestSimulation(
            name=name,
            exe_dict=targets,
            exfunc=eval_results,
            idxsim=0,
            make_comparison=False,
        ),
        str(ws),
    )

    # extract model objects
    gwf = sim.get_model(gwfname)
    prt = sim.get_model(prtname)

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

    # check budget data were written to mf6 prt list file
    check_budget_data(ws / f"{name}_prt.lst", perlen, nper)

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
    plt.savefig(ws / f"test_{name}.png")

    # convert mf6 pathlines to mp7 format
    mf6_pldata_mp7 = to_mp7_format(mf6_pldata)

    # sort both dataframes by particleid and time
    mf6_pldata_mp7.sort_values(by=["particleid", "time"], inplace=True)
    mp7_pldata.sort_values(by=["particleid", "time"], inplace=True)

    # drop duplicate locations
    # (mp7 includes a duplicate location at the end of each pathline??)
    cols = ["particleid", "x", "y", "z", "time"]
    mp7_pldata = mp7_pldata.drop_duplicates(subset=cols)
    mf6_pldata_mp7 = mf6_pldata_mp7.drop_duplicates(subset=cols)

    # drop columns for which there is no direct correspondence between mf6 and mp7
    del mf6_pldata_mp7["sequencenumber"]
    del mf6_pldata_mp7["particleidloc"]
    del mf6_pldata_mp7["xloc"]
    del mf6_pldata_mp7["yloc"]
    del mf6_pldata_mp7["zloc"]
    del mp7_pldata["sequencenumber"]
    del mp7_pldata["particleidloc"]
    del mp7_pldata["xloc"]
    del mp7_pldata["yloc"]
    del mp7_pldata["zloc"]

    # compare mf6 / mp7 pathline data
    assert mf6_pldata_mp7.shape == mp7_pldata.shape
    assert np.allclose(mf6_pldata_mp7, mp7_pldata, atol=1e-3)
