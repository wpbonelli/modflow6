"""
Regression test for same_exchange_cells (DisConnExchange.f90), which matches a
GWT-GWT exchange to the GWF-GWF exchange it draws its flows from by comparing
connected cells rather than exchange order alone.

The comparison is only performed on the side(s) of the pair where both models
are local to the current process (see the routine's docstring). In serial,
every model is local, so both sides are always compared. In parallel, whether
a side is compared depends on how the GWF-GWF and GWT-GWT exchanges each order
their two models (exgmnamea/exgmnameb) -- nothing requires a GWT-GWT exchange
to list its two models in the same order as the GWF-GWF exchange between the
same model pair.

This test builds a GWT-GWT exchange whose model order is reversed relative to
the GWF-GWF exchange it should link to, and whose second connection is a
genuine mismatch (row 1 of one model paired with row 0 of the other, which
does not correspond to any real GWF-GWF connection). A correct match should be
rejected with "Cannot find GWF-GWF exchange", since the two exchanges do not
in fact connect the same cells.

Case:
  - gwtgwtorder : mismatched GWT-GWT exchange; expect the match to fail.
"""

from pathlib import Path

import flopy
import pytest
from framework import TestFramework

cases = ["gwtgwtorder"]

nlay, nrow, ncol = 1, 2, 2
delr, delc = 1.0, 1.0
top, botm = 1.0, [0.0]
strt = 1.0
hk = 1.0
porosity = 0.1
perlen, nstp = 1.0, 1

# real GWF-GWF connections: gwfa's rightmost column to gwfb's leftmost column,
# row by row
gwfgwf_exgdata = [
    [(0, i, ncol - 1), (0, i, 0), 1, 0.5 * delr, 0.5 * delr, delc, 0.0, delr]
    for i in range(nrow)
]

# GWT-GWT exchange with reversed model order (gwtb, gwta) relative to the
# GWF-GWF exchange (gwfa, gwfb). Row 0 is a genuine connection; row 1 pairs
# gwtb's row 1 cell with gwta's row 0 cell, which does not correspond to any
# real GWF-GWF connection.
gwtgwt_exgdata = [
    [(0, 0, 0), (0, 0, ncol - 1), 1, 0.5 * delr, 0.5 * delr, delc, 0.0, delr],
    [(0, 1, 0), (0, 0, ncol - 1), 1, 0.5 * delr, 0.5 * delr, delc, 0.0, delr],
]


def add_dis(model, dis, xorigin):
    dis(
        model,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        xorigin=xorigin,
    )


def build_flow_model(sim, name, xorigin, upgradient):
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, model_nam_file=f"{name}.nam"
    )
    add_dis(gwf, flopy.mf6.ModflowGwfdis, xorigin)
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, k33=hk)
    if upgradient:
        flopy.mf6.ModflowGwfwel(
            gwf,
            stress_period_data={0: [[(0, i, 0), 1.0, 1.0] for i in range(nrow)]},
            auxiliary="CONCENTRATION",
            pname="WEL-1",
        )
    else:
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data={0: [[(0, i, ncol - 1), 0.0] for i in range(nrow)]},
            pname="CHD-1",
        )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return gwf


def build_transport_model(sim, name, xorigin, upgradient):
    gwt = flopy.mf6.MFModel(
        sim, model_type="gwt6", modelname=name, model_nam_file=f"{name}.nam"
    )
    gwt.name_file.save_flows = True
    add_dis(gwt, flopy.mf6.ModflowGwtdis, xorigin)
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="upstream")
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    sources = [("WEL-1", "AUX", "CONCENTRATION")] if upgradient else None
    flopy.mf6.ModflowGwtssm(gwt, sources=sources)
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{name}.cbc",
        concentration_filerecord=f"{name}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
    )
    return gwt


def build_models(idx, test):
    name = cases[idx]
    nouter, ninner = 100, 300
    hclose, rclose = 1e-8, 1e-8

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(perlen, nstp, 1.0)]
    )

    gwfnames = ["gwfa", "gwfb"]
    build_flow_model(sim, gwfnames[0], 0.0, True)
    build_flow_model(sim, gwfnames[1], ncol * delr, False)
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_exgdata),
        exgmnamea=gwfnames[0],
        exgmnameb=gwfnames[1],
        exchangedata=gwfgwf_exgdata,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="gwfgwf.gwfgwf",
    )
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        filename="gwf.ims",
    )
    sim.register_ims_package(imsgwf, gwfnames)

    gwtnames = ["gwta", "gwtb"]
    build_transport_model(sim, gwtnames[0], 0.0, True)
    build_transport_model(sim, gwtnames[1], ncol * delr, False)
    # model order reversed relative to the GWF-GWF exchange above, with
    # GWFMODELNAME1/2 kept consistent with that reversal
    flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        gwfmodelname1=gwfnames[1],
        gwfmodelname2=gwfnames[0],
        adv_scheme="upstream",
        nexg=len(gwtgwt_exgdata),
        exgmnamea=gwtnames[1],
        exgmnameb=gwtnames[0],
        exchangedata=gwtgwt_exgdata,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="gwtgwt.gwtgwt",
    )
    for gwfname, gwtname in zip(gwfnames, gwtnames):
        flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=gwfname,
            exgmnameb=gwtname,
            filename=f"{gwtname}.gwfgwt",
        )
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        filename="gwt.ims",
    )
    sim.register_ims_package(imsgwt, gwtnames)

    return sim, None


def check_output(idx, test):
    # serial writes mfsim.lst; parallel writes one mfsim.p{rank}.lst per
    # process, and either rank may be the one that reports the error
    lst_files = sorted(Path(test.workspace).glob("mfsim*.lst"))
    assert lst_files, f"no mfsim*.lst found in {test.workspace}"
    text = " ".join(" ".join(fpth.read_text().split()) for fpth in lst_files)
    srch = "Cannot find GWF-GWF exchange"
    assert srch in text, (
        "expected the mismatched GWT-GWT exchange (row 1 does not connect "
        "the same cells as the GWF-GWF exchange) to be rejected, but no "
        f"'{srch}' error was reported in {[f.name for f in lst_files]} -- "
        "the GWT-GWT exchange was likely linked to the GWF-GWF exchange "
        "without its connections actually being compared"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        xfail=True,
    )
    test.run()
