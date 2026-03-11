"""
Test DISV grids where adjacent cells share a face with an extra vertex
splitting their shared boundary. This situation is arranged and tested
within a single grid, and in separate grids coupled by an exchange. In
the exchange scenario, 2 connections are configured on the split face.

Three test cases are included:
1. GWF-GWF exchange without interface model
2. GWF-GWF exchange with interface model enabled
3. Single DISV grid with 2 adjacent cells (no exchange)

The former two are not allowed, and MF6 should raise an error. The last
case should run successfully, as we currently tolerate such split faces
within a single grid.
"""

import flopy
import pytest
from framework import TestFramework

cases = [
    "gwfgwf_sf",
    "gwfgwf_sf_ifmod",
    "gwf_sf",
]
ifmod = [False, True, None]


def build_models(idx, test):
    nlay = 1
    delr = 1.0
    delc = 1.0
    top = 1.0
    botm = 0.0
    hk = 1.0
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=test.workspace,
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)]
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        outer_dvclose=1.0e-5,
        outer_maximum=100,
        under_relaxation="NONE",
        inner_maximum=100,
        inner_dvclose=1.0e-6,
        rcloserecord=0.1,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=0.99,
    )

    # single grid case with 2 adjacent cells sharing a face with an extra vertex
    if idx == 2:
        gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
        disv = flopy.mf6.ModflowGwfdisv(
            gwf,
            nlay=nlay,
            ncpl=2,
            nvert=7,
            top=top,
            botm=botm,
            vertices=[
                [0, 0.0, 0.0],  # v0: left cell bottom-left
                [1, 1.0, 0.0],  # v1: shared bottom
                [2, 1.0, 0.5],  # v2: shared middle (extra vertex)
                [3, 1.0, 1.0],  # v3: shared top
                [4, 0.0, 1.0],  # v4: left cell top-left
                [5, 2.0, 0.0],  # v5: right cell bottom-right
                [6, 2.0, 1.0],  # v6: right cell top-right
            ],
            cell2d=[
                [0, 0.5, 0.5, 5, 0, 4, 3, 2, 1],  # left
                [1, 1.5, 0.5, 5, 1, 2, 3, 6, 5],  # right
            ],
        )
        ic = flopy.mf6.ModflowGwfic(gwf, strt=1.0)
        npf = flopy.mf6.ModflowGwfnpf(
            gwf,
            save_flows=True,
            icelltype=0,
            k=hk,
        )
        chd = flopy.mf6.ModflowGwfchd(
            gwf, stress_period_data=[[(0, 0), 1.0], [(0, 1), 0.0]]
        )
        oc = flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord="gwf.hds",
            budget_filerecord="gwf.cbc",
            saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        )

        return sim, None

    # left model
    gwf1 = flopy.mf6.ModflowGwf(sim, modelname="gwf1", save_flows=True)
    vertices1 = [
        [0, 0.0, 0.0],  # v0: bottom-left
        [1, 1.0, 0.0],  # v1: bottom-right
        [2, 1.0, 0.5],  # v2: mid-right (extra vertex)
        [3, 1.0, 1.0],  # v3: top-right
        [4, 0.0, 1.0],  # v4: top-left
    ]
    disv1 = flopy.mf6.ModflowGwfdisv(
        gwf1,
        nlay=nlay,
        ncpl=1,
        nvert=5,
        top=top,
        botm=botm,
        vertices=vertices1,
        # [icell2d, xc, yc, nvert, v0, v1, v2, ...]
        cell2d=[[0, 0.5, 0.5, 5, 0, 4, 3, 2, 1]],
    )
    ic1 = flopy.mf6.ModflowGwfic(gwf1, strt=1.0)
    npf1 = flopy.mf6.ModflowGwfnpf(
        gwf1,
        save_flows=True,
        icelltype=0,
        k=hk,
    )
    chd1 = flopy.mf6.ModflowGwfchd(gwf1, stress_period_data=[[(0, 0), 1.0]])
    oc1 = flopy.mf6.ModflowGwfoc(
        gwf1,
        head_filerecord="gwf1.hds",
        budget_filerecord="gwf1.cbc",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # right model
    gwf2 = flopy.mf6.ModflowGwf(sim, modelname="gwf2", save_flows=True)
    disv2 = flopy.mf6.ModflowGwfdisv(
        gwf2,
        nlay=nlay,
        ncpl=1,
        nvert=5,
        top=top,
        botm=botm,
        vertices=[
            [0, 1.0, 0.0],  # v0: bottom-left
            [1, 2.0, 0.0],  # v1: bottom-right
            [2, 2.0, 1.0],  # v2: top-right
            [3, 1.0, 1.0],  # v3: top-left
            [4, 1.0, 0.5],  # v4: mid-left (extra vertex)
        ],
        # [icell2d, xc, yc, nvert, v0, v1, v2, ...]
        cell2d=[[0, 1.5, 0.5, 5, 0, 4, 3, 2, 1]],
    )
    ic2 = flopy.mf6.ModflowGwfic(gwf2, strt=0.0)
    npf2 = flopy.mf6.ModflowGwfnpf(
        gwf2,
        save_flows=True,
        icelltype=0,
        k=hk,
    )
    chd2 = flopy.mf6.ModflowGwfchd(gwf2, stress_period_data=[[(0, 0), 0.0]])
    oc2 = flopy.mf6.ModflowGwfoc(
        gwf2,
        head_filerecord="gwf2.hds",
        budget_filerecord="gwf2.cbc",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # create exchange with 2 connections, 1 for each segment of the split face
    angldegx = 0.0
    cdist = delr
    gwfgwf_data = [
        [
            (0, 0),  # cell in model 1 (layer, cellid)
            (0, 0),  # cell in model 2 (layer, cellid)
            1,  # ihc (horizontal connection)
            0.5,  # cl1 (distance from gwf1 cell center to face)
            0.5,  # cl2 (distance from gwf2 cell center to face)
            delc / 2.0,  # hwva (flow width = half cell height)
            angldegx,  # auxiliary variable
            cdist,  # auxiliary variable
        ],
        [
            (0, 0),  # cell in model 1 (layer, cellid)
            (0, 0),  # cell in model 2 (layer, cellid)
            1,  # ihc (horizontal connection)
            0.5,  # cl1 (distance from gwf1 cell center to face)
            0.5,  # cl2 (distance from gwf2 cell center to face)
            delc / 2.0,  # hwva (flow width = half cell height)
            angldegx,  # auxiliary variable
            cdist,  # auxiliary variable
        ],
    ]
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        save_flows=True,
        print_flows=True,
        nexg=len(gwfgwf_data),
        exgmnamea="gwf1",
        exgmnameb="gwf2",
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        dev_interfacemodel_on=ifmod[idx],
    )

    return sim, None


def check_output(idx, test):
    if idx == 2:
        # single grid case
        fpth = test.workspace / "gwf.hds"
        hds = flopy.utils.HeadFile(fpth)
        heads = hds.get_data()

        head1 = heads[0, 0, 0]
        head2 = heads[0, 0, 1]

        assert 0.9 < head1 < 1.1
        assert -0.1 < head2 < 0.1
        assert head1 > head2
    else:
        # exchange cases
        pass


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail="gwfgwf" in name,
    )
    test.run()
