"""
Test that a package listed more than once in the PACKAGES block of a
model mover or an exchange mover is reported as an input error, and
that packages with the same name in different models are accepted.
"""

import subprocess

import flopy
import pytest

nlay, nrow, ncol = 1, 1, 4
delr, delc = 10.0, 10.0
errmsg = "listed more than once in PACKAGES block."


def add_model(sim, name, mover=True):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=-10.0,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 0.0]])
    flopy.mf6.ModflowGwfwel(
        gwf,
        mover=mover,
        pname="wel",
        stress_period_data=[[(0, 0, ncol - 1), -1.0]],
    )
    flopy.mf6.ModflowGwfdrn(
        gwf,
        mover=mover,
        pname="drn",
        stress_period_data=[[(0, 0, ncol - 1), 1.0, 1.0]],
    )
    return gwf


def build_model_mvr(ws, exe, duplicate):
    name = "mvr"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim)
    gwf = add_model(sim, name)
    packages = [["wel"], ["drn"]]
    if duplicate:
        packages.append(["WEL"])
    flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=1,
        maxpackages=len(packages),
        packages=packages,
    )
    return sim


def build_exchange_mvr(ws, exe, duplicate):
    sim = flopy.mf6.MFSimulation(sim_name="exgmvr", sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim)
    ims = flopy.mf6.ModflowIms(sim)
    left = add_model(sim, "left")
    right = add_model(sim, "right")
    sim.register_ims_package(ims, [left.name, right.name])
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=1,
        exgmnamea="left",
        exgmnameb="right",
        exchangedata=[[(0, 0, ncol - 1), (0, 0, 0), 1, delr / 2, delr / 2, delc]],
    )
    # the same package name in two models is two different packages
    packages = [["left", "wel"], ["right", "wel"]]
    if duplicate:
        packages.append(["left", "WEL"])
    gwfgwf.mvr.initialize(
        modelnames=True,
        maxmvr=1,
        maxpackages=len(packages),
        packages=packages,
        filename="left-right.exg.mvr",
    )
    return sim


@pytest.mark.parametrize("duplicate", [False, True])
@pytest.mark.parametrize("mover_type", ["model", "exchange"])
def test_mvr_duplicate_packages(function_tmpdir, targets, mover_type, duplicate):
    mf6 = targets["mf6"]
    build = build_model_mvr if mover_type == "model" else build_exchange_mvr
    sim = build(function_tmpdir, mf6, duplicate)
    sim.write_simulation()

    proc = subprocess.run([mf6], cwd=function_tmpdir, capture_output=True, text=True)
    print(proc.stdout)
    if duplicate:
        assert proc.returncode != 0
        assert errmsg in proc.stdout
    else:
        assert proc.returncode == 0
        assert errmsg not in proc.stdout
