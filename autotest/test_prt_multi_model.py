"""
Test that multiple PRT models can run in the same simulation
as a GWF flow model without crashing.

Regression test for GitHub issue #2811. The bug was that PRT
models shared module-level singleton method objects (MethodPool,
MethodSubcellPool). When two PRT models were in the same simulation,
each model's destructor tried to free the shared singletons, causing
a double-free crash at teardown.

Fixed in commit c80fee3a (refactor(prt): give model ownership of
method hierarchy, #2736), which replaced shared singleton pools with
per-model-instance method object ownership.

The grid is a 10x10 square with a single layer, the same flow system
shown on the FloPy readme. Two PRT models are added to the same
simulation, each with its own particle release points.
"""

from pathlib import Path

import flopy
import pytest
from framework import TestFramework
from prt_test_utils import (
    FlopyReadmeCase,
    check_budget_data,
    check_track_data,
)

simname = "prtmulti"


def build_models(idx, test):
    name = test.name
    sim = FlopyReadmeCase.get_gwf_sim(name, test.workspace, test.targets["mf6"])
    gwf_name = f"{name}_gwf"

    for i in (1, 2):
        prt_name = f"{name}_prt{i}"
        prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

        flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
            prt,
            pname="dis",
            nlay=FlopyReadmeCase.nlay,
            nrow=FlopyReadmeCase.nrow,
            ncol=FlopyReadmeCase.ncol,
            top=FlopyReadmeCase.top,
            botm=FlopyReadmeCase.botm,
        )

        flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)

        flopy.mf6.ModflowPrtprp(
            prt,
            pname="prp1",
            filename=f"{prt_name}.prp",
            nreleasepts=len(FlopyReadmeCase.releasepts_prt),
            packagedata=FlopyReadmeCase.releasepts_prt,
            perioddata={0: ["FIRST"]},
            extend_tracking=True,
        )

        flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            track_filerecord=[f"{prt_name}.trk"],
            trackcsv_filerecord=[f"{prt_name}.trk.csv"],
        )

        flopy.mf6.ModflowGwfprt(
            sim,
            exgtype="GWF6-PRT6",
            exgmnamea=gwf_name,
            exgmnameb=prt_name,
            filename=f"{gwf_name}_prt{i}.gwfprt",
        )

        ems = flopy.mf6.ModflowEms(sim, pname=f"ems{i}", filename=f"{prt_name}.ems")
        sim.register_solution_package(ems, [prt_name])

    return sim, None


def check_output(idx, test):
    name = test.name
    ws = Path(test.workspace)

    for i in (1, 2):
        prt_name = f"{name}_prt{i}"
        assert (ws / f"{prt_name}.trk").is_file()
        assert (ws / f"{prt_name}.trk.csv").is_file()
        check_budget_data(
            ws / f"{prt_name}.lst",
            FlopyReadmeCase.perlen,
            FlopyReadmeCase.nper,
        )
        check_track_data(
            track_bin=ws / f"{prt_name}.trk",
            track_hdr=ws / f"{prt_name}.trk.hdr",
            track_csv=ws / f"{prt_name}.trk.csv",
        )


@pytest.mark.parametrize("idx, name", [(0, simname)])
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
