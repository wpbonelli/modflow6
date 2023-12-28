"""
Tests ability to run a GWF model then a PRT model
in separate simulations via flow model interface,
as well as

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme.

Test cases are defined for the particle release
package (PRP) option STOP_AT_WEAK_SINK, one on
and one with the option off. No effect on results
is expected, because the model has no weak sinks.
(Motivated by an old bug in which particles were
tracked improperly when this option was enabled,
even with no weak sink cells in the vicinity.)

This test also specifies `boundnames=True` for
the PRP package, but does not provide boundnames
values, and checks that the "name" column in the
track output files contain the expected defaults.

Particles are released from the top left cell.

Pathlines are compared against a MODPATH 7 model.
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
from prt_test_utils import (
    BasicDisCase,
    all_equal,
    check_budget_data,
    check_track_data,
    get_model_name,
    get_partdata,
    has_default_boundnames,
)

from framework import TestFramework

simname = "prtfmi01"
cases = [simname, f"{simname}saws", f"{simname}bprp"]


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
        nper=BasicDisCase.nper,
        perioddata=[
            (BasicDisCase.perlen, BasicDisCase.nstp, BasicDisCase.tsmult)
        ],
    )

    # create prt model
    prtname = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)

    # create prt discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=BasicDisCase.nlay,
        nrow=BasicDisCase.nrow,
        ncol=BasicDisCase.ncol,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=BasicDisCase.porosity)

    # convert mp7 to prt release points and check against expectation
    partdata = get_partdata(prt.modelgrid, BasicDisCase.releasepts_mp7)
    coords = partdata.to_coords(prt.modelgrid)
    if "bprp" in name:
        # bad cell indices!
        releasepts = [
            (i, 0, 1, 1, c[0], c[1], c[2]) for i, c in enumerate(coords)
        ]
    else:
        releasepts = [
            (i, 0, 0, 0, c[0], c[1], c[2]) for i, c in enumerate(coords)
        ]
        assert np.allclose(BasicDisCase.releasepts_prt, releasepts)

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
    gwfname = get_model_name(name, "gwf")
    gwf_budget_file = gwf_ws / f"{gwfname}.bud"
    gwf_head_file = gwf_ws / f"{gwfname}.hds"
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


def build_mp7_sim(name, ws, mp7, gwf):
    partdata = get_partdata(gwf.modelgrid, BasicDisCase.releasepts_mp7)
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
        porosity=BasicDisCase.porosity,
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
    gwfsim = BasicDisCase.get_gwf_sim(
        test.name, test.workspace, test.targets.mf6
    )
    prtsim = build_prt_sim(
        test.name, test.workspace, test.workspace / "prt", test.targets.mf6
    )
    return gwfsim, prtsim


def check_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    mp7_ws = test.workspace / "mp7"
    gwfname = get_model_name(name, "gwf")
    prtname = get_model_name(name, "prt")
    mp7name = get_model_name(name, "mp7")

    if "bprp" in name:
        buff = test.buffs[1]
        assert any("Error: release point" in l for l in buff)
        return

    # extract mf6 simulations/models and grid
    gwfsim = test.sims[0]
    prtsim = test.sims[1]
    gwf = gwfsim.get_model(gwfname)
    prt = prtsim.get_model(prtname)
    mg = gwf.modelgrid

    # build/run mp7 model... can't run within framework as
    # flopy needs gwf output files to write mp7 input files
    mp7sim = build_mp7_sim(name, mp7_ws, test.targets.mp7, gwf)

    # run mp7 model
    mp7sim.write_input()
    success, buff = mp7sim.run_model(report=True)
    assert success, pformat(buff)

    # check mf6 output files exist
    gwf_budget_file = f"{gwfname}.bud"
    gwf_head_file = f"{gwfname}.hds"
    prt_track_file = f"{prtname}.trk"
    prt_track_csv_file = f"{prtname}.trk.csv"
    prp_track_file = f"{prtname}.prp.trk"
    prp_track_csv_file = f"{prtname}.prp.trk.csv"
    assert (gwf_ws / gwf_budget_file).is_file()
    assert (gwf_ws / gwf_head_file).is_file()
    assert (prt_ws / prt_track_file).is_file()
    assert (prt_ws / prt_track_csv_file).is_file()
    assert (prt_ws / prp_track_file).is_file()
    assert (prt_ws / prp_track_csv_file).is_file()

    # check mp7 output files exist
    mp7_pathline_file = f"{mp7name}.mppth"
    assert (mp7_ws / mp7_pathline_file).is_file()

    # load mp7 pathline results
    plf = PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # make sure pathline df has "name" (boundname) column and default values
    assert "name" in mf6_pls
    assert has_default_boundnames(mf6_pls)

    # make sure all mf6 pathline data have correct model and PRP index (1)
    assert all_equal(mf6_pls["imdl"], 1)
    assert all_equal(mf6_pls["iprp"], 1)

    # check budget data were written to mf6 prt list file
    check_budget_data(
        prt_ws / f"{name}_prt.lst", BasicDisCase.perlen, BasicDisCase.nper
    )

    # check mf6 prt particle track data were written to binary/CSV files
    # and that different formats are equal
    for track_csv in [
        prt_ws / prt_track_csv_file,
        prt_ws / prp_track_csv_file,
    ]:
        check_track_data(
            track_bin=prt_ws / prt_track_file,
            track_hdr=prt_ws
            / Path(prt_track_file.replace(".trk", ".trk.hdr")),
            track_csv=track_csv,
        )

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
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
    plt.savefig(gwf_ws / f"test_{simname}.png")

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


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        compare=None,
        xfail=[False, "bprp" in name],
    )
    test.run()
