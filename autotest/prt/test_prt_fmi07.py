"""
Tests ability to run a GWF model then a PRT model
in separate simulations via flow model interface,
with release points improperly mapped to cell IDs
(expect failures).

The grid is a 10x10 square with a single layer,
the same flow system shown on the FloPy readme.

Particles are released from the top left cell.
"""


from pprint import pformat

import flopy
import pytest
from prt_test_utils import get_gwf_sim, get_model_name, get_partdata

simname = "prtfmi07"
ex = [simname]


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
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=ctx.nlay,
        nrow=ctx.nrow,
        ncol=ctx.ncol,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=ctx.porosity)

    # convert mp7 to prt release points and check against expectation
    partdata = get_partdata(prt.modelgrid, ctx.releasepts_mp7)
    coords = partdata.to_coords(prt.modelgrid)
    # bad cell indices!
    releasepts = [(i, 0, 1, 1, c[0], c[1], c[2]) for i, c in enumerate(coords)]

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

    return sim


@pytest.mark.parametrize("name", ex)
def test_mf6model(name, function_tmpdir, targets):
    # workspace
    ws = function_tmpdir

    # build mf6 models
    gwfsim, ctx = get_gwf_sim(name, ws, targets.mf6)
    prtsim = build_prt_sim(ctx, name, ws, targets.mf6)

    # run gwf models
    gwfsim.write_simulation()
    success, buff = gwfsim.run_simulation(report=True)
    assert success, pformat(buff)

    # run prt model (expect failure)
    prtsim.write_simulation()
    success, buff = prtsim.run_simulation(report=True)
    assert not success, pformat(buff)
    assert any("Error: release point" in l for l in buff)
