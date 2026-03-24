"""
Test PRT-PRT exchange between two particle tracking models.

The two model domains form a 10-cell line split down the middle:

        left        |        right
  0   1   2   3   4 | 0   1   2   3   4
 _______________________________________
| * | ->| ->| ->| ->| ->| ->| ->| ->|   |
|___|___|___|___|___|___|___|___|___|___|
  *particle release

Steady flow runs from left to right.

PRT models are connected by an exchange at (0, 0, 4) <-> (0, 0, 0).

Particles are released from the left-most cell (0, 0, 0) in the left
model. They should reach and cross the exchange and terminate in the
rightmost cell of the right model.  The test asserts that:
"""

from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

nlay = 1
nrow = 1
ncol_half = 5
ncol = ncol_half * 2

delr = 1.0
delc = 1.0
delz = 1.0
top = 1.0
botm = [0.0]

head_left = 1.0
head_right = 0.0
hk = 1.0

porosity = 0.1
nper = 1
perlen = 10.0
nstp = 10

releasepts = [(0, (0, 0, 0), 0.5, 0.5, 0.5)]

cases = ["prtprtexg"]


def get_model_name(idx, mdl):
    return f"{cases[idx]}_{mdl}"


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    mf6 = test.targets["mf6"]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        sim_ws=ws,
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(perlen, nstp, 1.0)],
    )
    flopy.mf6.ModflowIms(
        sim,
        pname="ims",
        complexity="simple",
        outer_dvclose=1e-8,
        inner_dvclose=1e-8,
    )

    gwfl = flopy.mf6.ModflowGwf(
        sim,
        modelname=get_model_name(idx, "gwfl"),
        save_flows=True,
    )
    flopy.mf6.ModflowGwfdis(
        gwfl,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_half,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwfl, strt=head_left)
    flopy.mf6.ModflowGwfnpf(
        gwfl,
        save_specific_discharge=True,
        icelltype=0,
        k=hk,
    )
    flopy.mf6.ModflowGwfchd(
        gwfl,
        stress_period_data=[[(0, 0, 0), head_left]],
        pname="chd_left",
    )
    flopy.mf6.ModflowGwfoc(
        gwfl,
        head_filerecord=f"{get_model_name(idx, 'gwfl')}.hds",
        budget_filerecord=f"{get_model_name(idx, 'gwfl')}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwfr = flopy.mf6.ModflowGwf(
        sim,
        modelname=get_model_name(idx, "gwfr"),
        save_flows=True,
    )
    flopy.mf6.ModflowGwfdis(
        gwfr,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_half,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwfr, strt=head_right)
    flopy.mf6.ModflowGwfnpf(
        gwfr,
        save_specific_discharge=True,
        icelltype=0,
        k=hk,
    )
    # CHD on right face
    flopy.mf6.ModflowGwfchd(
        gwfr,
        stress_period_data=[[(0, 0, ncol_half - 1), head_right]],
        pname="chd_right",
    )
    flopy.mf6.ModflowGwfoc(
        gwfr,
        head_filerecord=f"{get_model_name(idx, 'gwfr')}.hds",
        budget_filerecord=f"{get_model_name(idx, 'gwfr')}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwfgwf_data = [
        (
            (0, 0, ncol_half - 1),
            (0, 0, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            0.0,
            delr,
        )
    ]
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=get_model_name(idx, "gwfl"),
        exgmnameb=get_model_name(idx, "gwfr"),
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="gwf-gwf.exg",
    )

    prtl = flopy.mf6.ModflowPrt(sim, modelname=get_model_name(idx, "prtl"))
    flopy.mf6.ModflowPrtdis(
        prtl,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_half,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowPrtmip(prtl, pname="mip", porosity=porosity)
    flopy.mf6.ModflowPrtprp(
        prtl,
        pname="prp1",
        filename=f"{get_model_name(idx, 'prtl')}.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(
        prtl,
        pname="oc",
        trackcsv_filerecord=f"{get_model_name(idx, 'prtl')}.trk.csv",
    )
    flopy.mf6.ModflowPrtfmi(
        prtl,
        packagedata=[
            ("GWFHEAD", f"{get_model_name(idx, 'gwfl')}.hds"),
            ("GWFBUDGET", f"{get_model_name(idx, 'gwfl')}.cbc"),
        ],
    )
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=get_model_name(idx, "gwfl"),
        exgmnameb=get_model_name(idx, "prtl"),
        filename=f"{get_model_name(idx, 'gwfl')}.gwfprt",
    )

    prtr = flopy.mf6.ModflowPrt(sim, modelname=get_model_name(idx, "prtr"))
    flopy.mf6.ModflowPrtdis(
        prtr,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_half,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowPrtmip(prtr, pname="mip", porosity=porosity)
    # right model has no PRP, it just receives particles from left
    flopy.mf6.ModflowPrtoc(
        prtr,
        pname="oc",
        trackcsv_filerecord=f"{get_model_name(idx, 'prtr')}.trk.csv",
    )
    flopy.mf6.ModflowPrtfmi(
        prtr,
        packagedata=[
            ("GWFHEAD", f"{get_model_name(idx, 'gwfr')}.hds"),
            ("GWFBUDGET", f"{get_model_name(idx, 'gwfr')}.cbc"),
        ],
    )
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=get_model_name(idx, "gwfr"),
        exgmnameb=get_model_name(idx, "prtr"),
        filename=f"{get_model_name(idx, 'gwfr')}.gwfprt",
    )

    prtprt_data = [
        # nodem1, nodem2, ihc, iflowface1, iflowface2
        ((0, 0, ncol_half - 1), (0, 0, 0), 1, 3, 1),
    ]
    flopy.mf6.ModflowPrtprt(
        sim,
        exgtype="PRT6-PRT6",
        nexg=len(prtprt_data),
        exgmnamea=get_model_name(idx, "prtl"),
        exgmnameb=get_model_name(idx, "prtr"),
        gwfmodelname1=get_model_name(idx, "gwfl"),
        gwfmodelname2=get_model_name(idx, "gwfr"),
        exchangedata=prtprt_data,
        filename="prt-prt.exg",
        auxiliary=["IFLOWFACE1", "IFLOWFACE2"],
    )

    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename="prt.ems",
    )
    sim.register_solution_package(
        ems, [get_model_name(idx, "prtl"), get_model_name(idx, "prtr")]
    )

    return sim


def check_output(idx, test):
    ws = Path(test.workspace)

    hds_l = flopy.utils.HeadFile(ws / f"{get_model_name(idx, 'gwfl')}.hds")
    hds_r = flopy.utils.HeadFile(ws / f"{get_model_name(idx, 'gwfr')}.hds")
    head_l = hds_l.get_data().squeeze()
    head_r = hds_r.get_data().squeeze()

    assert np.all(head_l >= 0.5) and np.all(head_l <= head_left)
    assert np.all(head_r >= head_right) and np.all(head_r <= 0.5)

    pls_l = pd.read_csv(ws / f"{get_model_name(idx, 'prtl')}.trk.csv")
    pls_r = pd.read_csv(ws / f"{get_model_name(idx, 'prtr')}.trk.csv")
    assert len(pls_l) > 0
    assert len(pls_r) > 0

    irpts_l = pls_l["irpt"].unique()
    irpts_r = pls_r["irpt"].unique()
    assert set(irpts_l) == set(irpts_r)

    terminations = pls_r[pls_r.ireason == 3]
    assert len(terminations) == len(irpts_l)
    assert all(terminations.istatus == 5)

    final_x = terminations["x"].max()
    assert np.isclose(final_x, delr * ncol)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()
