from math import sqrt
import os
from pathlib import Path
from typing import Optional
import flopy
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from modflow_devtools.case import Case
from flopy.utils.binaryfile import write_budget, write_head
from flopy.utils.gridutil import uniform_flow_field
from flopy.utils.gridintersect import GridIntersect
from pytest_cases import parametrize_with_cases, parametrize
from shapely.geometry import MultiPoint, LineString
from framework import TestFramework
from simulation import TestSimulation


# expected format of particle track data
dtype = {
    "names": ["kper", "kstp", "iprp", "irpt", "icell", "izone", "istatus", "ireason", "trelease", "t", "x", "y", "z"],
    "formats": ["<i4", "<i4", "<i4", "<i4", "<i4", "<i4", "<i4", "<i4", "<f8", "<f8", "<f8", "<f8", "<f8"],
}

def load_track_data(cbb):
    """
    Temporary method to read particle track data from a cell budget file,
    until the dedicated output file method is completed in mf6 and FloPy.
    """

    upn = cbb.get_unique_package_names(decode=True)
    times = cbb.get_times()
    tracks = []
    for prpnam in upn:
        for totim in times:
            data = cbb.get_data(text='DATA-PRTCL', paknam=prpnam, totim=totim)
            for ploc in data[0]:
                kper, kstp, iprp, ip, icell, izone, istatus, ireason, trelease, t, x, y, z, = [ploc[i] for i in range(3, 16)]
                pathpoint = (kper, kstp, iprp, ip, icell, izone, istatus, ireason, trelease, t, x, y, z)
                tracks.append(pathpoint)
        break

    return np.core.records.fromrecords(tracks, dtype=dtype)

def get_track_dtype(path: os.PathLike):
    """Get the dtype of the track data recarray from the ascii header file."""

    hdr_lns = open(path).readlines()
    hdr_lns_spl = [[ll.strip() for ll in l.split(',')] for l in hdr_lns]
    return np.dtype(list(zip(hdr_lns_spl[0], hdr_lns_spl[1])))

def check_track_dtype(data: np.recarray):
    """Check that the dtype of the track data recarray is correct."""

    assert np.array_equal(data.dtype.names, dtype["names"])
    assert np.array_equal([v[0] for v in data.dtype.fields.values()], dtype["formats"])
    assert data.dtype == dtype

def check_track_data(
        track_bin: os.PathLike,
        track_hdr: os.PathLike,
        track_csv: os.PathLike,
        track_bud: Optional[os.PathLike] = None):
    """Check that track data written to binary, CSV, and budget files are equal."""   

    # get dtype from ascii header file
    dt = get_track_dtype(track_hdr)

    # read output files
    data_bin = np.fromfile(track_bin, dtype=dt)
    data_csv = np.genfromtxt(track_csv, dtype=dt, delimiter=',', skip_header=True)
    if track_bud:
        data_bud = load_track_data(flopy.utils.CellBudgetFile(track_bud))

    # check dtypes
    check_track_dtype(data_bin)
    check_track_dtype(data_csv)
    if track_bud:
        check_track_dtype(data_bud)

    assert data_bin.shape == data_csv.shape, f"Binary and CSV track data shapes do not match: {data_bin.shape} != {data_csv.shape}"
    if track_bud:
        assert data_bin.shape == data_bud.shape, f"Binary and budget track data shapes do not match: {data_bin.shape} != {data_bud.shape}"

    # check particle tracks written to all output files are equal
    # check each column separately to avoid:
    # TypeError: The DType <class 'numpy._FloatAbstractDType'> could not be promoted by <class 'numpy.dtype[void]'>
    for k in data_bin.dtype.names:
        assert np.allclose(data_bin[k], data_csv[k], equal_nan=True)
        if track_bud:
            assert np.allclose(data_bin[k], data_bud[k], equal_nan=True)

def check_budget_data(ctx, lst: os.PathLike, cbb: os.PathLike):
    # load PRT model's list file
    mflist = flopy.utils.mflistfile.ListBudget(
        lst,
        budgetkey="MASS BUDGET FOR ENTIRE MODEL"
    )
    names = mflist.get_record_names()
    entries = mflist.entries

    # check timesteps
    inc = mflist.get_incremental()
    v = inc["totim"][-1]
    assert v == ctx.perlen, f"Last time should be {ctx.perlen}.  Found {v}"

    # entries should be a subset of names
    assert all(e in names for e in entries)

    # todo what other record names should we expect?
    expected_entries = [
        "PRP_IN",
        "PRP_OUT",
    ]
    assert all(en in names for en in expected_entries)

    # load and check cell budget file
    mfbud = flopy.utils.binaryfile.CellBudgetFile(cbb)
    assert mfbud.nlay == ctx.nlay
    assert mfbud.nrow == ctx.nrow
    assert mfbud.ncol == ctx.ncol
    assert len(mfbud.times) == 1
    assert mfbud.times[0] == ctx.perlen


def test_fmi_basic(function_tmpdir, targets):
    """
    Tests ability to run a GWF model then a PRT model
    in separate simulations via flow model interface.
    Also tests for correct output files, and compares
    data written to output files of different formats.
    """

    name = "prtfmi0"
    mf6 = targets.mf6
    ws = function_tmpdir

    def run_flow_model():
        sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=mf6)
        pd = [(1.0, 1, 1.0), (1.0, 1, 1.0)]
        tdis = flopy.mf6.ModflowTdis(sim, nper=len(pd), perioddata=pd)
        ims = flopy.mf6.ModflowIms(sim)
        gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
        dis = flopy.mf6.ModflowGwfdis(gwf, nrow=10, ncol=10)
        ic = flopy.mf6.ModflowGwfic(gwf)
        npf = flopy.mf6.ModflowGwfnpf(
            gwf, save_specific_discharge=True, save_saturation=True
        )
        spd = {
            0: [[(0, 0, 0), 1.0, 1.0], [(0, 9, 9), 0.0, 0.0]],
            1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
        }
        chd = flopy.mf6.ModflowGwfchd(
            gwf, pname="CHD-1", stress_period_data=spd, auxiliary=["concentration"]
        )
        budget_file = f"{name}.bud"
        head_file = f"{name}.hds"
        oc = flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=budget_file,
            head_filerecord=head_file,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )
        sim.write_simulation()
        success, _ = sim.run_simulation()
        assert success

        # check output files
        assert (ws / budget_file).is_file()
        assert (ws / head_file).is_file()

    def run_tracking_model():
        sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=mf6)
        pd = [(1.0, 10, 1.0), (1.0, 10, 1.0)]
        tdis = flopy.mf6.ModflowTdis(sim, nper=len(pd), perioddata=pd)
        prt = flopy.mf6.ModflowPrt(sim, modelname=name)
        dis = flopy.mf6.ModflowGwfdis(prt, nrow=10, ncol=10)
        mip = flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=.2)
        releasepts = [
            # particle id, k, i, j, localx, localy, localz
            (0, 0, 0, 0, 0.5, 0.5, 0.5)
        ]
        prp_track_file = f"{name}.prp.trk"
        prp_track_csv_file = f"{name}.prp.trk.csv"
        prp = flopy.mf6.ModflowPrtprp(
            prt, pname="prp1", filename=f"{name}_1.prp",
            track_filerecord=[prp_track_file],
            trackcsv_filerecord=[prp_track_csv_file],
            nreleasepts=len(releasepts), packagedata=releasepts,
            perioddata={0: ["FIRST"]},
        )
        prt_budget_file = f"{name}.cbb"
        prt_track_file = f"{name}.trk"
        prt_track_csv_file = f"{name}.trk.csv"
        oc = flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            budget_filerecord=[prt_budget_file],
            track_filerecord=[prt_track_file],
            trackcsv_filerecord=[prt_track_csv_file],
            saverecord=[("BUDGET", "ALL")],
        )
        gwf_budget_file = f"{name}.bud"
        gwf_head_file = f"{name}.hds"
        pd = [
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ]
        fmi = flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)
        ems = flopy.mf6.ModflowEms(
            sim, pname="ems",
            filename=f"{name}.ems",
        )
        sim.register_solution_package(ems, [prt.name])
        sim.write_simulation()
        success, _ = sim.run_simulation()
        assert success

        # make sure output files exist
        assert (ws / prt_budget_file).is_file()
        assert (ws / prt_track_file).is_file()
        assert (ws / prt_track_csv_file).is_file()
        assert (ws / prp_track_file).is_file()
        assert (ws / prp_track_csv_file).is_file()

        check_track_data(
            track_bin=ws / prt_track_file,
            track_hdr=ws / Path(prt_track_file.replace(".trk", ".trk.hdr")),
            track_csv=ws / prt_track_csv_file,
            track_bud=ws / prt_budget_file
        )

    run_flow_model()
    run_tracking_model()


class PrtCases:
    """
    Test cases for particle tracking models.
    """

    # prt alone, consuming gwf outputs via flow model interface
    cases_fmi = [
        Case(
            name="prtfmi1",
            perlen=10,
            dtmin=1.001e-5,
            nlay=1,
            nrow=1,
            ncol=3,
            nper=1,
            nstp=1,
            tsmult=1.0,
            steady=[True],
            delr=1.0,
            delc=1.0,
            top=1.0,
            laytyp=0,
            ss=0.0,
            sy=0.1,
            botm=[0.0],
            strt=1.0,
            hnoflo=1e30,
            hdry=-1e30,
            hk=1.0,
            porosity=0.1,
        ),
    ]

    @parametrize(ctx=cases_fmi, ids=[c.name for c in cases_fmi])
    def case_fmi(self, ctx, function_tmpdir, targets):
        workspace = Path(function_tmpdir)
        # create a heads file (with head equal top) for the prt model to consume
        headfile_name = f"{ctx.name}.hds"
        with open(workspace / headfile_name, "wb") as fbin:
            for kstp in range(ctx.nstp):
                write_head(fbin, ctx.top * np.ones((ctx.nrow, ctx.ncol)), kstp=kstp + 1)
        
        # create a budget file for the prt model to consume
        qx = 0.0
        qy = 0.0
        qz = 0.0
        shape = (ctx.nlay, ctx.nrow, ctx.ncol)
        spdis, flowja = uniform_flow_field(qx, qy, qz, shape)
        dt = np.dtype(
            [
                ("ID1", np.int32),
                ("ID2", np.int32),
                ("FLOW", np.float64),
                ("SATURATION", np.float64),
            ]
        )
        sat = np.array(
            [(i, i, 0.0, 1.0) for i in range(ctx.nlay * ctx.nrow * ctx.ncol)], dtype=dt
        )
        budgetfile_name = f"{ctx.name}.bud"
        with open(workspace / budgetfile_name, "wb") as fbin:
            for kstp in range(ctx.nstp):
                write_budget(fbin, flowja, kstp=kstp + 1)
                write_budget(
                    fbin, spdis, text="      DATA-SPDIS", imeth=6, kstp=kstp + 1
                )
                write_budget(
                    fbin, sat, text="        DATA-SAT", imeth=6, kstp=kstp + 1
                )

        # create simulation
        sim = flopy.mf6.MFSimulation(
            sim_name=ctx.name,
            exe_name=targets.mf6,
            version="mf6",
            sim_ws=workspace
        )

        # create tdis package
        pd = (ctx.perlen, ctx.nstp, ctx.tsmult)
        flopy.mf6.modflow.mftdis.ModflowTdis(
            sim, pname="tdis", time_units="DAYS", nper=ctx.nper, perioddata=[pd]
        )

        # create prt model
        prt = flopy.mf6.ModflowPrt(sim, modelname=ctx.name)

        # create prt discretization
        flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
            prt,
            pname="dis",
            nlay=ctx.nlay,
            nrow=ctx.nrow,
            ncol=ctx.ncol,
            length_units="FEET",
            delr=ctx.delr,
            delc=ctx.delc,
            top=ctx.top,
            botm=ctx.botm,
        )

        # create mip package
        flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=ctx.porosity)

        # create prp package
        releasepts = [
            # particle id, k, i, j, localx, localy, localz
            (0, 0, 0, 0, 0.5, 0.5, 0.5)
        ]
        flopy.mf6.ModflowPrtprp(
            prt, pname="prp1", filename=f"{ctx.name}_1.prp",
            nreleasepts=len(releasepts), packagedata=releasepts,
            perioddata={0: ["FIRST"]},
        )

        # create output control package
        flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            budget_filerecord=[f"{ctx.name}.cbb"],
            track_filerecord=[f"{ctx.name}.trk"],
            trackcsv_filerecord=[f"{ctx.name}.trk.csv"],
            saverecord=[("BUDGET", "ALL")],
        )

        # create the flow model interface
        pd = [
            ("GWFHEAD", headfile_name),
            ("GWFBUDGET", budgetfile_name),
        ]
        flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)

        # add explicit model solution
        ems = flopy.mf6.ModflowEms(
            sim, pname="ems",
            filename=f"{ctx.name}.ems",
        )
        sim.register_solution_package(ems, [prt.name])

        return ctx, sim, None, self.eval_fmi
    
    def eval_fmi(self, ctx, sim):
        print(f"Evaluating results for sim {sim.name}")
        simpath = Path(sim.simpath)
        
        # check budget data
        check_budget_data(
            ctx,
            simpath / f"{sim.name}.lst",
            simpath / f"{sim.name}.cbb")

        # check particle track data
        prt_track_file = simpath / f"{sim.name}.trk"
        prt_track_hdr_file = simpath / f"{sim.name}.trk.hdr"
        prt_track_csv_file = simpath / f"{sim.name}.trk.csv"
        prt_track_bud_file = simpath / f"{sim.name}.cbb"
        assert prt_track_file.exists()
        assert prt_track_hdr_file.exists()
        assert prt_track_csv_file.exists()
        check_track_data(
            track_bin=prt_track_file,
            track_hdr=prt_track_hdr_file,
            track_csv=prt_track_csv_file,
            track_bud=prt_track_bud_file
        )


    # gwf+prt models in same simulation via exchange
    cases_exg = [
        Case(
            name="prtexg1",
            nlay=1,
            nrow=10,
            ncol=10,
            top=1.0,
            botm=[0.0],
            nper=1,
            perlen=1.0,
            nstp=1,
            tsmult=1.0,
            porosity=0.1
        )
    ]

    @parametrize(ctx=cases_exg, ids=[c.name for c in cases_exg])
    def case_exg(self, ctx, function_tmpdir, targets):
        workspace = function_tmpdir
        gwfname = f"{ctx.name}"
        prtname = f"{ctx.name}_prt"

        # create simulation
        sim = flopy.mf6.MFSimulation(
            sim_name=ctx.name,
            exe_name=targets.mf6,
            version="mf6",
            sim_ws=workspace
        )

        # create tdis package
        pd = (ctx.perlen, ctx.nstp, ctx.tsmult)
        flopy.mf6.modflow.mftdis.ModflowTdis(
            sim,
            pname="tdis",
            time_units="DAYS",
            nper=ctx.nper,
            perioddata=[pd]
        )

        # create gwf model
        gwf = flopy.mf6.ModflowGwf(
            sim,
            modelname=gwfname,
            save_flows=True
        )

        # create gwf discretization
        flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
            gwf,
            pname="dis",
            nlay=ctx.nlay,
            nrow=ctx.nrow,
            ncol=ctx.ncol,
        )

        # create gwf initial conditions package
        flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic")

        # create gwf node property flow package
        flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            gwf,
            pname="npf", 
            save_saturation=True,
            save_specific_discharge=True,
        )

        # create gwf chd package
        spd = {
            0: [[(0, 0, 0), 1.0, 1.0], [(0, 9, 9), 0.0, 0.0]],
            1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
        }
        chd = flopy.mf6.ModflowGwfchd(
            gwf,
            pname="CHD-1",
            stress_period_data=spd,
            auxiliary=["concentration"]
        )

        # create gwf output control package
        budget_file = f"{gwfname}.bud"
        head_file = f"{gwfname}.hds"
        oc = flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=budget_file,
            head_filerecord=head_file,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        # create iterative model solution for gwf model
        ims = flopy.mf6.ModflowIms(sim)

        # create prt model
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
        releasepts = [
            # particle id, k, i, j, localx, localy, localz
            (0, 0, 0, 0, 0.5, 0.5, 0.5)
        ]
        flopy.mf6.ModflowPrtprp(
            prt, pname="prp1", filename=f"{prtname}_1.prp",
            nreleasepts=len(releasepts), packagedata=releasepts,
            perioddata={0: ["FIRST"]},
        )

        # create output control package
        flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            budget_filerecord=[f"{prtname}.cbb"],
            track_filerecord=[f"{prtname}.trk"],
            trackcsv_filerecord=[f"{prtname}.trk.csv"],
            saverecord=[("BUDGET", "ALL")],
        )

        # create the flow model interface
        pd = [
            ("GWFHEAD", f"{gwfname}.hds"),
            ("GWFBUDGET", f"{gwfname}.bud"),
        ]
        flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)

        # create exchange
        flopy.mf6.ModflowGwfprt(
            sim, exgtype="GWF6-PRT6",
            exgmnamea=gwfname, exgmnameb=prtname,
            filename=f"{gwfname}.gwfprt",
        )

        # add explicit model solution
        ems = flopy.mf6.ModflowEms(
            sim, pname="ems",
            filename=f"{prtname}.ems",
        )
        sim.register_solution_package(ems, [prt.name])

        return ctx, sim, None, self.eval_exg
    
    def eval_exg(self, ctx, sim):
        print(f"Evaluating results for sim {sim.name}")
        simpath = Path(sim.simpath)
        
        # check budget data
        check_budget_data(
            ctx,
            simpath / f"{sim.name}_prt.lst",
            simpath / f"{sim.name}_prt.cbb")

        # check particle track data
        prt_track_file = simpath / f"{sim.name}_prt.trk"
        prt_track_hdr_file = simpath / f"{sim.name}_prt.trk.hdr"
        prt_track_csv_file = simpath / f"{sim.name}_prt.trk.csv"
        prt_track_bud_file = simpath / f"{sim.name}_prt.cbb"
        assert prt_track_file.exists()
        assert prt_track_hdr_file.exists()
        assert prt_track_csv_file.exists()
        check_track_data(
            track_bin=prt_track_file,
            track_hdr=prt_track_hdr_file,
            track_csv=prt_track_csv_file,
            track_bud=prt_track_bud_file
        )


    # MODPATH 7 example problem 1
    #
    # Subproblems 1A and 1B are solved by a single MODFLOW 6 run and a single MODPATH 7 run,
    # so they are included under one "scenario"; each of the two subproblems is represented
    # by its own particle release package (for MODFLOW 6) or particle group (for MODPATH 7)
    case_mp7_p01 = Case(
        name="prtmp7p01",
        # discretization
        length_units = "feet",
        time_units = "days",
        nper=1,
        nstp=1,
        perlen=1.0,
        tsmult=1.0,
        nlay=3,
        nrow=21,
        ncol=20,
        delr=500,
        delc=500,
        top=400.0,
        botm=[220.0, 200.0, 0.0],
        laytyp=[1, 0, 0],
        # conductivities
        kh=[50.0, 0.01, 200.0],
        kv=[10.0, 0.01, 20.0],
        # well
        wel_loc=(2, 10, 9),
        wel_q=-150000.0,
        # recharge
        rch=0.005,
        rch_iface=6,
        rch_iflowface=-1,
        # river
        riv_h=320.0,
        riv_z=317.0,
        riv_c=1.0e5,
        riv_iface=6,
        riv_iflowface=-1,
        # particle tracking params
        porosity=0.1,
        def_zone=1,
        wel_zone=2,
    )
    cases_mp7_p01 = [case_mp7_p01]
    
    @parametrize(ctx=cases_mp7_p01, ids=[c.name for c in cases_mp7_p01])
    def case_mp7_p01(self, ctx, function_tmpdir, targets):
        nm_mf6 = ctx.name
        nm_prt = ctx.name + "_prt"

        headfile = "{}.hds".format(nm_mf6)
        budgetfile = "{}.cbb".format(nm_mf6)
        budgetfile_prt = "{}.cbb".format(nm_prt)

        # Well package
        wd = [(ctx.wel_loc, ctx.wel_q)]

        # River package
        rd = [[
            (0, i, ctx.ncol - 1),
            ctx.riv_h,
            ctx.riv_c,
            ctx.riv_z,
            ctx.riv_iface,
            ctx.riv_iflowface
        ] for i in range(ctx.nrow)]

        zones_lay1 = ctx.def_zone
        zones_lay2 = ctx.def_zone
        zones_lay3 = np.full((ctx.nrow, ctx.ncol), ctx.def_zone, dtype=np.int32)
        zones_lay3[ctx.wel_loc[1:]] = ctx.wel_zone

        # Starting location template size for example 1B;
        # in the original example, there is initially a
        # 3x3 array of particles in each cell in layer 1;
        # in this example, there is initially one
        # particle in each cell in layer 1; the original
        # 3x3 particle arrays can be restored simply by
        # setting sloc_tmpl_size below to 3 instead of 1.
        sloc_tmpl_size = 1

        # Zones
        zones = [zones_lay1, zones_lay2, zones_lay3]

        # Default iface
        defaultiface = {"RCH": 6, "EVT": 6}

        # Example 1A release points
        releasepts = {}
        releasepts['1A'] = []
        zrpt = ctx.top
        k = 0
        j = 2
        for i in range(ctx.nrow):
            nrpt = i
            xrpt = (j + 0.5) * ctx.delr
            yrpt = (ctx.nrow - i - 0.5) * ctx.delc
            rpt = [nrpt, k, i, j, xrpt, yrpt, zrpt]
            releasepts['1A'].append(rpt)

        # Example 1B release points
        releasepts['1B'] = []
        ndivc = sloc_tmpl_size
        ndivr = sloc_tmpl_size
        deldivc = ctx.delc / ndivc
        deldivr = ctx.delr / ndivr
        k = 0
        zrpt = ctx.top
        nrpt = -1
        for i in range(ctx.nrow):
            y0 = (ctx.nrow - i - 1) * ctx.delc
            for j in range(ctx.ncol):
                x0 = j * ctx.delr
                for idiv in range(ndivc):
                    dy = (idiv + 0.5) * deldivc
                    yrpt = y0 + dy
                    for jdiv in range(ndivr):
                        dx = (jdiv + 0.5) * deldivr
                        xrpt = x0 + dx
                        nrpt += 1
                        rpt = [nrpt, k, i, j, xrpt, yrpt, zrpt]
                        releasepts['1B'].append(rpt)
        
        # Get well and river cell numbers
        nodes = {}
        k, i, j = ctx.wel_loc
        nodes['well'] = ctx.ncol * (ctx.nrow * k + i) + j
        nodes['river'] = []
        for rivspec in rd:
            k, i, j = rivspec[0]
            node = ctx.ncol * (ctx.nrow * k + i) + j
            nodes['river'].append(node)
        
        # Instantiate the MODFLOW 6 simulation object
        sim = flopy.mf6.MFSimulation(
            sim_name=nm_mf6, exe_name=targets.mf6, version="mf6", sim_ws=function_tmpdir
        )

        # Instantiate the MODFLOW 6 temporal discretization package
        pd = (ctx.perlen, ctx.nstp, ctx.tsmult)
        flopy.mf6.modflow.mftdis.ModflowTdis(
            sim, pname="tdis", time_units="DAYS", nper=ctx.nper, perioddata=[pd]
        )

        # Instantiate the MODFLOW 6 gwf (groundwater-flow) model
        model_nam_file = "{}.nam".format(nm_mf6)
        gwf = flopy.mf6.ModflowGwf(
            sim, modelname=nm_mf6, model_nam_file=model_nam_file, save_flows=True
        )

        # Instantiate the MODFLOW 6 gwf discretization package
        flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
            gwf,
            pname="dis",
            nlay=ctx.nlay,
            nrow=ctx.nrow,
            ncol=ctx.ncol,
            length_units="FEET",
            delr=ctx.delr,
            delc=ctx.delc,
            top=ctx.top,
            botm=ctx.botm,
        )

        # Instantiate the MODFLOW 6 gwf initial conditions package
        flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic", strt=ctx.top)

        # Instantiate the MODFLOW 6 gwf node property flow package
        flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            gwf, pname="npf", icelltype=ctx.laytyp, k=ctx.kh, k33=ctx.kv,
            save_saturation=True, save_specific_discharge=True,
        )

        # Instantiate the MODFLOW 6 gwf recharge package
        flopy.mf6.modflow.mfgwfrcha.ModflowGwfrcha(
            gwf, recharge=ctx.rch,
            auxiliary=["iface", "iflowface"], aux=[ctx.rch_iface, ctx.rch_iflowface],
        )

        # Instantiate the MODFLOW 6 gwf well package
        flopy.mf6.modflow.mfgwfwel.ModflowGwfwel(
            gwf, maxbound=1, stress_period_data={0: wd}
        )

        # Instantiate the MODFLOW 6 gwf river package
        flopy.mf6.modflow.mfgwfriv.ModflowGwfriv(
            gwf, auxiliary=["iface", "iflowface"], stress_period_data={0: rd}
        )

        # Instantiate the MODFLOW 6 gwf output control package
        head_record = [headfile]
        budget_record = [budgetfile]
        saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")]
        flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
            gwf,
            pname="oc",
            saverecord=saverecord,
            head_filerecord=head_record,
            budget_filerecord=budget_record,
        )

        # Instantiate the MODFLOW 6 prt model
        prt = flopy.mf6.ModflowPrt(
            sim, modelname=nm_prt, model_nam_file="{}.nam".format(nm_prt)
        )

        # Instantiate the MODFLOW 6 prt discretization package
        flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
            prt, pname="dis",
            nlay=ctx.nlay, nrow=ctx.nrow, ncol=ctx.ncol,
            length_units="FEET",
            delr=ctx.delr, delc=ctx.delc,
            top=ctx.top, botm=ctx.botm,
        )

        # Instantiate the MODFLOW 6 prt model input package
        flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=ctx.porosity)

        # Instantiate the MODFLOW 6 prt particle release point (prp) package
        # for example 1A
        nreleasepts1a = len(releasepts['1A'])
        pd = {0: ["FIRST"],}
        track_csv_1a = "{}_1a.trk.csv".format(nm_prt)
        flopy.mf6.ModflowPrtprp(
            prt, pname="prp1a", filename="{}_1a.prp".format(nm_prt),
            trackcsv_filerecord=[track_csv_1a],
            nreleasepts=nreleasepts1a, packagedata=releasepts['1A'],
            perioddata=pd,
        )

        # Instantiate the MODFLOW 6 prt particle release point (prp) package
        # for example 1B
        nreleasepts1b = len(releasepts['1B'])
        pd = {0: ["FIRST"],}
        track_csv_1b = "{}_1b.trk.csv".format(nm_prt)
        flopy.mf6.ModflowPrtprp(
            prt, pname="prp1b", filename="{}_1b.prp".format(nm_prt),
            trackcsv_filerecord=[track_csv_1b],
            nreleasepts=nreleasepts1b, packagedata=releasepts['1B'],
            perioddata=pd,
        )

        # Instantiate the MODFLOW 6 prt output control package
        budget_record = [budgetfile_prt]
        track_record = [f"{nm_prt}.trk"]
        trackcsv_record = [f"{nm_prt}.trk.csv"]
        flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            budget_filerecord=budget_record,
            track_filerecord=track_record,
            trackcsv_filerecord=trackcsv_record,
            saverecord=[("BUDGET", "ALL")],
        )

        # Instantiate the MODFLOW 6 prt flow model interface
        flopy.mf6.ModflowPrtfmi(prt)

        # Create the MODFLOW 6 gwf-prt model exchange
        flopy.mf6.ModflowGwfprt(
            sim, exgtype="GWF6-PRT6",
            exgmnamea=nm_mf6, exgmnameb=nm_prt,
            filename="{}.gwfprt".format(nm_mf6),
        )

        # Create an iterative model solution (IMS) for the MODFLOW 6 gwf model
        ims = flopy.mf6.ModflowIms(
            sim, pname="ims",
            complexity="SIMPLE",
            outer_dvclose=1e-6, inner_dvclose=1e-6,
            rcloserecord=1e-6,
        )

        # Create an explicit model solution (EMS) for the MODFLOW 6 prt model
        ems = flopy.mf6.ModflowEms(
            sim,
            pname="ems",
            filename="{}.ems".format(nm_prt),
        )
        sim.register_solution_package(ems, [prt.name])

        return ctx, sim, None, self.eval_mp7_p01

    def eval_mp7_p01(self, ctx, sim):
        print(f"Evaluating results for sim {sim.name}")
        simpath = Path(sim.simpath)
        
        # check budget data
        check_budget_data(
            ctx,
            simpath / f"{sim.name}_prt.lst",
            simpath / f"{sim.name}_prt.cbb")

        # check model-level particle track data
        prt_track_bin_file = simpath / f"{sim.name}_prt.trk"
        prt_track_hdr_file = simpath / f"{sim.name}_prt.trk.hdr"
        prt_track_csv_file = simpath / f"{sim.name}_prt.trk.csv"
        prt_track_bud_file = simpath / f"{sim.name}_prt.cbb"
        assert prt_track_bin_file.exists()
        assert prt_track_hdr_file.exists()
        assert prt_track_csv_file.exists()
        check_track_data(
            track_bin=prt_track_bin_file,
            track_hdr=prt_track_hdr_file,
            track_csv=prt_track_csv_file,
            # kluge: not checking cbb because track data are overwritten
            # for each PRP, so it doesn't contain the full track dataset
            # track_bud=prt_track_bud_file
        )

        # check PRP-level particle track data
        track_csv_file_1a = Path(sim.simpath) / f"{ctx.name}_prt_1a.trk.csv"
        track_csv_file_1b = Path(sim.simpath) / f"{ctx.name}_prt_1b.trk.csv"
        assert track_csv_file_1a.is_file()
        assert track_csv_file_1b.is_file()

        # stack PRP-level track data, then
        # check equal to model-level track data
        # get dtype from ascii header file
        dt = get_track_dtype(prt_track_hdr_file)
        prt_data_bin = np.fromfile(prt_track_bin_file, dtype=dt)
        prt_data_csv = np.genfromtxt(prt_track_csv_file, dtype=dt, delimiter=',', skip_header=True)
        prp_data_csv_1a = np.genfromtxt(track_csv_file_1a, dtype=dt, delimiter=',', skip_header=True)
        prp_data_csv_1b = np.genfromtxt(track_csv_file_1b, dtype=dt, delimiter=',', skip_header=True)
        prp_data_csv = np.hstack((prp_data_csv_1a, prp_data_csv_1b))
        assert prt_data_bin.shape == prt_data_csv.shape == prp_data_csv.shape
        for k in prt_data_csv.dtype.names:
            assert np.allclose(prt_data_csv[k], prp_data_csv[k], equal_nan=True)

    # MODPATH 7 example problem 2
    case_mp7_p02 = Case(
        name="prtmp7p02",
        # discretization
        length_units = "feet",
        time_units = "days",
        tdis_rc = [(1000.0, 1, 1.0)],
        nper = 1,
        Lx = 10000.0,
        Ly = 10500.0,
        nlay = 3,
        nrow = 21,
        ncol = 20,
        delr = 500.0,
        delc = 500.0,
        top = 400,
        botm = [220, 200, 0],
        ncpl = 651,
        nvert = 723,
        # Cell types by layer
        icelltype = [1, 0, 0],
        # Conductivities
        k = [50.0, 0.01, 200.0],
        k33 = [10.0, 0.01, 20.0],
        # Well
        wel_coords = [(4718.45, 5281.25)],
        wel_q = [-150000.0],
        # Recharge
        rch = 0.005,
        rch_iface = 6,
        rch_iflowface = -1,
        # River
        riv_h = 320.0,
        riv_z = 318.0,
        riv_c = 1.0e5,
        riv_iface = 6,
        riv_iflowface = -1,
        # particle tracking
        porosity = 0.1
    )
    cases_mp7_p02 = [
        case_mp7_p02.copy_update(
            name=case_mp7_p02.name + "a",
        ),
        case_mp7_p02.copy_update(
            name=case_mp7_p02.name + "b",
        )
    ]

    @staticmethod
    def get_refined_gridprops_p02(ws: Path, ctx, exe) -> dict:
        ms = flopy.modflow.Modflow()
        dis = flopy.modflow.ModflowDis(
            ms,
            nlay=ctx.nlay,
            nrow=ctx.nrow,
            ncol=ctx.ncol,
            delr=ctx.delr,
            delc=ctx.delc,
            top=ctx.top,
            botm=ctx.botm,
        )

        from flopy.utils.gridgen import Gridgen

        # create Gridgen workspace
        gridgen_ws = ws / "gridgen"
        gridgen_ws.mkdir(parents=True, exist_ok=True)

        # create Gridgen object
        g = Gridgen(ms.modelgrid, model_ws=gridgen_ws, exe_name=exe)

        # add polygon for each refinement level
        outer_polygon = [
            [(3500, 4000), (3500, 6500), (6000, 6500), (6000, 4000), (3500, 4000)]
        ]
        g.add_refinement_features([outer_polygon], "polygon", 1, range(ctx.nlay))
        refshp0 = gridgen_ws / "rf0"

        middle_polygon = [
            [(4000, 4500), (4000, 6000), (5500, 6000), (5500, 4500), (4000, 4500)]
        ]
        g.add_refinement_features([middle_polygon], "polygon", 2, range(ctx.nlay))
        refshp1 = gridgen_ws / "rf1"

        inner_polygon = [
            [(4500, 5000), (4500, 5500), (5000, 5500), (5000, 5000), (4500, 5000)]
        ]
        g.add_refinement_features([inner_polygon], "polygon", 3, range(ctx.nlay))
        refshp2 = gridgen_ws / "rf2"

        g.build(verbose=False)

        return g.get_gridprops_disv()
    
    @parametrize(ctx=cases_mp7_p02, ids=[c.name for c in cases_mp7_p02])
    def case_mp7_p02(self, ctx, function_tmpdir, targets):

        sim_name = "mp7-p02"
        gwf_name = sim_name
        prt_name = sim_name + "_prt"
        mp7_name = sim_name + "_mp7"
        example_name = "ex-prt-" + sim_name

        headfile = "{}.hds".format(gwf_name)
        headfile_bkwd = "{}.hds".format(gwf_name + "_bkwd")
        budgetfile = "{}.cbb".format(gwf_name)
        budgetfile_bkwd = "{}.cbb".format(gwf_name + "_bkwd")
        budgetfile_prt = "{}.cbb".format(prt_name)

        welcells = []
        rivcells = []
        releasepts = []

        tdis_rc = [(1000.0, 1, 1.0)]
        nper = len(tdis_rc)

        length_units = "feet"
        time_units = "days"

        # Cell types by layer
        icelltype = [1, 0, 0]

        # Conductivities
        k = [50.0, 0.01, 200.0]
        k33 = [10.0, 0.01, 20.0]

        # Well
        # issue(flopyex): in the original flopy example, the well
        #   coordinates were at the edge (not the center) of the cell,
        #   but it apparently worked out anyway
        wel_coords = [(4718.45, 5281.25)]
        wel_q = [-150000.0]

        # Recharge
        rch = 0.005
        rch_iface = 6
        rch_iflowface = -1

        # River
        riv_h = 320.0
        riv_z = 318.0
        riv_c = 1.0e5
        riv_iface = 6
        riv_iflowface = -1

        # porosity
        porosity = 0.1
        
        # todo build/refine grid
        disv_props = PrtCases.get_refined_gridprops_p02(function_tmpdir / "gridgen", ctx, targets.gridgen)

        # retrieve GRIDGEN-generated gridprops
        ncpl = disv_props["ncpl"]
        top = disv_props["top"]
        botm = disv_props["botm"]
        nvert = disv_props["nvert"]
        vertices = disv_props["vertices"]
        cell2d = disv_props["cell2d"]

        # determine whether this is part A or B
        part = ctx.name[-1]

        # Instantiate the MODFLOW 6 simulation object
        sim = flopy.mf6.MFSimulation(
            sim_name=gwf_name, exe_name=targets.mf6, version="mf6", sim_ws=function_tmpdir
        )

        # Instantiate the MODFLOW 6 temporal discretization package
        flopy.mf6.ModflowTdis(
            sim, pname="tdis", time_units="DAYS", perioddata=tdis_rc, nper=len(tdis_rc)
        )

        # Instantiate the MODFLOW 6 gwf (groundwater-flow) model
        gwf = flopy.mf6.ModflowGwf(
            sim, modelname=gwf_name, model_nam_file="{}.nam".format(gwf_name)
        )
        gwf.name_file.save_flows = True

        # Instantiate the MODFLOW 6 gwf discretization package
        flopy.mf6.ModflowGwfdisv(
            gwf,
            length_units=length_units,
            **disv_props,
        )

        # GridIntersect object for setting up boundary conditions
        ix = GridIntersect(gwf.modelgrid, method="vertex", rtree=True)

        # Instantiate the MODFLOW 6 gwf initial conditions package
        flopy.mf6.ModflowGwfic(gwf, pname="ic", strt=riv_h)

        # Instantiate the MODFLOW 6 gwf node property flow package
        flopy.mf6.ModflowGwfnpf(
            gwf,
            xt3doptions=[("xt3d")],
            icelltype=icelltype,
            k=k, k33=k33,
            save_saturation=True, save_specific_discharge=True,
        )

        # Instantiate the MODFLOW 6 gwf recharge package
        flopy.mf6.ModflowGwfrcha(
            gwf, recharge=rch,
            auxiliary=["iface", "iflowface"], aux=[rch_iface, rch_iflowface],
        )

        # Instantiate the MODFLOW 6 gwf well package
        welcells = ix.intersects(MultiPoint(wel_coords))
        welcells = [icpl for (icpl,) in welcells]
        welspd = [[(2, icpl), wel_q[idx]] for idx, icpl in enumerate(welcells)]
        flopy.mf6.ModflowGwfwel(
            gwf, print_input=True, stress_period_data=welspd
        )

        # Instantiate the MODFLOW 6 gwf river package
        riverline = [(ctx.Lx - 1.0, ctx.Ly), (ctx.Lx - 1.0, 0.0)]
        rivcells = ix.intersects(LineString(riverline))
        rivcells = [icpl for (icpl,) in rivcells]
        rivspd = [[(0, icpl), riv_h, riv_c, riv_z, riv_iface, riv_iflowface]
                for icpl in rivcells]
        flopy.mf6.ModflowGwfriv(gwf, stress_period_data=rivspd,
            auxiliary=[("iface", "iflowface")])

        # Instantiate the MODFLOW 6 gwf output control package
        headfile = "{}.hds".format(gwf_name)
        head_record = [headfile]
        budgetfile = "{}.cbb".format(gwf_name)
        budget_record = [budgetfile]
        flopy.mf6.ModflowGwfoc(
            gwf,
            pname="oc",
            budget_filerecord=budget_record,
            head_filerecord=head_record,
            headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
            printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        # Create an iterative model solution (IMS) for the MODFLOW 6 gwf model
        ims = flopy.mf6.ModflowIms(
            sim,
            pname="ims",
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
        sim.register_ims_package(ims, [gwf.name])

        # set up PRT release point data
        releasepts = []
        welcellnode = welcells[0]
        xctr = cell2d[welcellnode][1]
        yctr = cell2d[welcellnode][2]
        vert0 = cell2d[welcellnode][4]
        vert1 = cell2d[welcellnode][5]
        x0 = vertices[vert0][1]
        y0 = vertices[vert0][2]
        x1 = vertices[vert1][1]
        y1 = vertices[vert1][2]
        dx = x1 - x0
        dy = y1 - y0
        delcell = sqrt(dx * dx + dy * dy)
        delcellhalf = 0.5 * delcell
        if part == "a":
            zrpt = 0.5 * (botm[1][welcellnode] + botm[2][welcellnode])
            npart_on_side = 4
            delta = delcell / npart_on_side
            baseoffset = 0.5 * (npart_on_side + 1) * delta
            xbase = xctr - baseoffset
            ybase = yctr - baseoffset
            nrpt = -1
            for idiv in range(npart_on_side):
                i = idiv + 1
                xrpt = xbase + i * delta
                yrpt = yctr + delcellhalf
                nrpt += 1
                rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                releasepts.append(rpt)
                yrpt = yctr - delcellhalf
                nrpt += 1
                rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                releasepts.append(rpt)
                yrpt = ybase + i * delta
                xrpt = xctr + delcellhalf
                nrpt += 1
                rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                releasepts.append(rpt)
                xrpt = xctr - delcellhalf
                nrpt += 1
                rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                releasepts.append(rpt)
        else:
            # 4x4 array of particles on top of well cell
            zrpt = botm[1][welcellnode]
            npart_on_side = 4
            delta = delcell / npart_on_side
            baseoffset = 0.5 * (npart_on_side + 1) * delta
            xbase = xctr - baseoffset
            ybase = yctr - baseoffset
            nrpt = -1
            for idivx in range(npart_on_side):
                ix = idivx + 1
                xrpt = xbase + ix * delta
                for idivy in range(npart_on_side):
                    iy = idivy + 1
                    yrpt = ybase + iy * delta
                    nrpt += 1
                    rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                    releasepts.append(rpt)
            # 10x10 arrays of particles on the four sides of well cell
            zwel = 0.5 * (botm[1][welcellnode] + botm[2][welcellnode])
            npart_on_side = 10
            delcellz = botm[1][welcellnode] - botm[2][welcellnode]
            delta = delcell / npart_on_side
            deltaz = delcellz / npart_on_side
            baseoffset = 0.5 * (npart_on_side + 1) * delta
            baseoffsetz = 0.5 * (npart_on_side + 1) * deltaz
            xbase = xctr - baseoffset
            ybase = yctr - baseoffset
            zbase = zwel - baseoffsetz
            for idivz in range(npart_on_side):
                iz = idivz + 1
                zrpt = zbase + iz * deltaz
                for idiv in range(npart_on_side):
                    i = idiv + 1
                    xrpt = xbase + i * delta
                    yrpt = yctr + delcellhalf
                    nrpt += 1
                    rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                    releasepts.append(rpt)
                    yrpt = yctr - delcellhalf
                    nrpt += 1
                    rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                    releasepts.append(rpt)
                    yrpt = ybase + i * delta
                    xrpt = xctr + delcellhalf
                    nrpt += 1
                    rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                    releasepts.append(rpt)
                    xrpt = xctr - delcellhalf
                    nrpt += 1
                    rpt = [nrpt, (2, welcellnode), xrpt, yrpt, zrpt]
                    releasepts.append(rpt)

        nodew = ncpl * 2 + welcells[0]

        # set up MODPATH 7 particle release points
        if part == "a":
            pcoord = np.array(
                [
                    [0.000, 0.125, 0.500],
                    [0.000, 0.375, 0.500],
                    [0.000, 0.625, 0.500],
                    [0.000, 0.875, 0.500],
                    [1.000, 0.125, 0.500],
                    [1.000, 0.375, 0.500],
                    [1.000, 0.625, 0.500],
                    [1.000, 0.875, 0.500],
                    [0.125, 0.000, 0.500],
                    [0.375, 0.000, 0.500],
                    [0.625, 0.000, 0.500],
                    [0.875, 0.000, 0.500],
                    [0.125, 1.000, 0.500],
                    [0.375, 1.000, 0.500],
                    [0.625, 1.000, 0.500],
                    [0.875, 1.000, 0.500],
                ]
            )
            plocs = [nodew for i in range(pcoord.shape[0])]
            pgdata = flopy.modpath.ParticleData(
                plocs,
                structured=False,
                localx=pcoord[:, 0],
                localy=pcoord[:, 1],
                localz=pcoord[:, 2],
                drape=0,
            )
        else:
            facedata = flopy.modpath.FaceDataType(
                drape=0,
                verticaldivisions1=10,
                horizontaldivisions1=10,
                verticaldivisions2=10,
                horizontaldivisions2=10,
                verticaldivisions3=10,
                horizontaldivisions3=10,
                verticaldivisions4=10,
                horizontaldivisions4=10,
                rowdivisions5=0,
                columndivisions5=0,
                rowdivisions6=4,
                columndivisions6=4,
            )
            pgdata = flopy.modpath.NodeParticleData(subdivisiondata=facedata, nodes=nodew)

        # Set particle release point data according to the scenario
        prpname = f"prp2{part}"
        prpfilename = f"{prt_name}_2{part}.prp"
        nreleasepts = len(releasepts)

        # Create a time discretization for backward tracking;
        # since there is only a single period with a single time step,
        # the time discretization is the same backward as forward
        tdis_bkwd = tdis_rc

        # Instantiate the MODFLOW 6 prt model
        prt = flopy.mf6.ModflowPrt(
            sim, modelname=prt_name, model_nam_file="{}.nam".format(prt_name)
        )

        # Instantiate the MODFLOW 6 prt discretization package
        flopy.mf6.ModflowGwfdisv(
            prt,
            length_units=length_units,
            **disv_props,
        )

        # Instantiate the MODFLOW 6 prt model input package
        flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

        # Instantiate the MODFLOW 6 prt particle release point (prp) package
        pd = {0: ["FIRST"], 1: []}
        flopy.mf6.ModflowPrtprp(
            prt, pname=prpname, filename=prpfilename,
            nreleasepts=nreleasepts, packagedata=releasepts,
            perioddata=pd, 
        )

        # Instantiate the MODFLOW 6 prt output control package
        budgetfile_prt = "{}.cbb".format(prt_name)
        budget_record = [budgetfile_prt]
        flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            budget_filerecord=budget_record,
            saverecord=[("BUDGET", "ALL")],
        )

        # Instantiate the MODFLOW 6 prt flow model interface
        # using "time-reversed" budget and head files
        pd = [
            ("GWFHEAD", headfile_bkwd),
            ("GWFBUDGET", budgetfile_bkwd),
        ]
        flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)

        # create exchange
        flopy.mf6.ModflowGwfprt(
            sim, exgtype="GWF6-PRT6",
            exgmnamea=gwf_name, exgmnameb=prt_name,
            filename=f"{prt_name}.gwfprt",
        )

        # Create an explicit model solution (EMS) for the MODFLOW 6 prt model
        ems = flopy.mf6.ModflowEms(
            sim, pname="ems",
            filename="{}.ems".format(prt_name),
        )
        sim.register_solution_package(ems, [prt.name])

        return ctx, sim, None, self.eval_mp7_p02


    def eval_mp7_p02(self, ctx, sim):
        pass

    # MODPATH 7 example problem 3
    cases_mp7_p03 = []

    @pytest.mark.skip(reason="indev: need finer-grained particle release timing")
    @parametrize(ctx=cases_mp7_p03, ids=[c.name for c in cases_mp7_p03])
    def case_mp7_p03(self, ctx, function_tmpdir, targets):
        pass

    def eval_mp7_p03(self, ctx, sim):
        pass

    # MODPATH 7 example problem 4
    cases_mp7_p04 = [
        Case(
            name="prtmp7p04",
            nlay=1,
            nrow=21,
            ncol=26,
            delr=500,
            delc=500,
            porosity=0.1,
            nper=1,
            tdis_rc=[(10000, 1, 1.0)],
            top=100.0,
            botm=np.zeros((1, 21, 26), dtype=np.float32),
        )
    ]

    @parametrize(ctx=cases_mp7_p04, ids=[c.name for c in cases_mp7_p04])
    def case_mp7_p04(self, ctx, function_tmpdir, targets):
        ms = flopy.modflow.Modflow()
        dis = flopy.modflow.ModflowDis(
            ms,
            nlay=ctx.nlay,
            nrow=ctx.nrow,
            ncol=ctx.ncol,
            delr=ctx.delr,
            delc=ctx.delc,
            top=ctx.top,
            botm=ctx.botm,
        )

        gwf_name = f"{ctx.name}_gwf"
        prt_name = f"{ctx.name}_prt"

        headfile = f"{gwf_name}.hds"
        headfile_bkwd = f"{gwf_name + '_bkwd'}.hds"
        budgetfile = f"{gwf_name}.bud"
        budgetfile_bkwd = f"{gwf_name + '_bkwd'}.bud"
        budgetfile_prt = f"{prt_name}.cbb"

        from flopy.utils.gridgen import Gridgen

        # create Gridgen workspace
        gridgen_ws = function_tmpdir / "gridgen"
        gridgen_ws.mkdir(parents=True, exist_ok=True)

        # create Gridgen object
        g = Gridgen(ms.modelgrid, model_ws=gridgen_ws, exe_name=targets.gridgen)

        # add polygon for each refinement level
        outer_polygon = [
            [
                (2500, 6000),
                (2500, 9500),
                (3000, 9500),
                (3000, 10000),
                (6000, 10000),
                (6000, 9500),
                (6500, 9500),
                (6500, 6000),
                (6000, 6000),
                (6000, 5500),
                (3000, 5500),
                (3000, 6000),
                (2500, 6000),
            ]
        ]
        g.add_refinement_features([outer_polygon], "polygon", 1, range(ctx.nlay))
        refshp0 = gridgen_ws / "rf0"

        middle_polygon = [
            [
                (3000, 6500),
                (3000, 9000),
                (3500, 9000),
                (3500, 9500),
                (5500, 9500),
                (5500, 9000),
                (6000, 9000),
                (6000, 6500),
                (5500, 6500),
                (5500, 6000),
                (3500, 6000),
                (3500, 6500),
                (3000, 6500),
            ]
        ]
        g.add_refinement_features([middle_polygon], "polygon", 2, range(ctx.nlay))
        refshp1 = gridgen_ws / "rf1"

        inner_polygon = [
            [
                (3500, 7000),
                (3500, 8500),
                (4000, 8500),
                (4000, 9000),
                (5000, 9000),
                (5000, 8500),
                (5500, 8500),
                (5500, 7000),
                (5000, 7000),
                (5000, 6500),
                (4000, 6500),
                (4000, 7000),
                (3500, 7000),
            ]
        ]
        g.add_refinement_features([inner_polygon], "polygon", 3, range(ctx.nlay))
        refshp2 = gridgen_ws / "rf2"

        g.build(verbose=False)
        grid = flopy.discretization.VertexGrid(**g.get_gridprops_vertexgrid())

        # fmt: off
        idomain = [
            0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,
            0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,
            1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
            0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,
            0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,
            0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0
        ]
        # fmt: on
        disv_props = g.get_gridprops_disv()

        # simulation
        sim = flopy.mf6.MFSimulation(
            sim_name=ctx.name, sim_ws=function_tmpdir, exe_name=targets.mf6, version="mf6"
        )

        # temporal discretization
        tdis = flopy.mf6.ModflowTdis(sim, time_units="days", nper=ctx.nper, perioddata=ctx.tdis_rc)

        # iterative model solver
        ims = flopy.mf6.ModflowIms(
            sim,
            pname="ims",
            complexity="SIMPLE",
            outer_dvclose=1e-4,
            outer_maximum=100,
            inner_dvclose=1e-5,
            under_relaxation_theta=0,
            under_relaxation_kappa=0,
            under_relaxation_gamma=0,
            under_relaxation_momentum=0,
            linear_acceleration="BICGSTAB",
            relaxation_factor=0.99,
            number_orthogonalizations=2,
        )

        # groundwater flow model
        gwf = flopy.mf6.ModflowGwf(
            sim, modelname=gwf_name, model_nam_file=f"{gwf_name}.nam", save_flows=True
        )

        # grid discretization
        
        disv = flopy.mf6.ModflowGwfdisv(gwf, length_units="feet", idomain=idomain, **disv_props)

        # initial conditions
        ic = flopy.mf6.ModflowGwfic(gwf, strt=150.0)

        # wells are tuples (layer, node number, q, iface, iflowface)
        # where iflowface is a PRT parameter corresponding to iface.
        wells = [
            # negative q: discharge
            (0, 861, -30000.0, 0, -1),
            (0, 891, -30000.0, 0, -1),
            # positive q: injection
            (0, 1959, 10000.0, 1, 4),
            (0, 1932, 10000.0, 3, 3),
            (0, 1931, 10000.0, 3, 3),
            (0, 1930, 5000.0, 1, 4),
            (0, 1930, 5000.0, 3, 3),
            (0, 1903, 5000.0, 1, 4),
            (0, 1903, 5000.0, 3, 3),
            (0, 1876, 10000.0, 3, 3),
            (0, 1875, 10000.0, 3, 3),
            (0, 1874, 5000.0, 1, 4),
            (0, 1874, 5000.0, 3, 3),
            (0, 1847, 10000.0, 3, 3),
            (0, 1846, 5000.0, 3, 3),
            (0, 1845, 5000.0, 1, 4),
            (0, 1845, 5000.0, 3, 3),
            (0, 1818, 5000.0, 1, 4),
            (0, 1818, 5000.0, 3, 3),
            (0, 1792, 10000.0, 1, 4),
            (0, 1766, 10000.0, 1, 4),
            (0, 1740, 5000.0, 1, 4),
            (0, 1740, 5000.0, 4, 1),
            (0, 1715, 5000.0, 1, 4),
            (0, 1715, 5000.0, 4, 1),
            (0, 1690, 10000.0, 1, 4),
            (0, 1646, 5000.0, 1, 4),
            (0, 1646, 5000.0, 4, 1),
            (0, 1549, 5000.0, 1, 4),
            (0, 1549, 5000.0, 4, 1),
            (0, 1332, 5000.0, 4, 1),
            (0, 1332, 5000.0, 1, 5),  # why does IFLOWFACE need to be 5 here, not 4?
            (0, 1021, 2500.0, 1, 4),
            (0, 1021, 2500.0, 4, 1),
            (0, 1020, 5000.0, 1, 5),  # "
            (0, 708, 2500.0, 1, 5),  # "
            (0, 708, 2500.0, 4, 1),
            (0, 711, 625.0, 1, 4),
            (0, 711, 625.0, 4, 1),
            (0, 710, 625.0, 1, 4),
            (0, 710, 625.0, 4, 1),
            (0, 409, 1250.0, 1, 4),
            (0, 407, 625.0, 1, 4),
            (0, 407, 625.0, 4, 1),
            (0, 402, 625.0, 1, 4),
            (0, 402, 625.0, 4, 1),
            (0, 413, 1250.0, 1, 4),
            (0, 411, 1250.0, 1, 4),
            (0, 203, 1250.0, 1, 5),  # "
            (0, 202, 1250.0, 1, 4),
            (0, 202, 1250.0, 4, 1),
            (0, 199, 2500.0, 1, 4),
            (0, 197, 1250.0, 1, 4),
            (0, 197, 1250.0, 4, 1),
            (0, 96, 2500.0, 1, 4),
            (0, 97, 1250.0, 1, 4),
            (0, 97, 1250.0, 4, 1),
            (0, 103, 1250.0, 1, 4),
            (0, 103, 1250.0, 4, 1),
            (0, 102, 1250.0, 1, 4),
            (0, 102, 1250.0, 4, 1),
            (0, 43, 2500.0, 1, 4),
            (0, 43, 2500.0, 4, 1),
            (0, 44, 2500.0, 1, 4),
            (0, 44, 2500.0, 4, 1),
            (0, 45, 5000.0, 4, 1),
            (0, 10, 10000.0, 1, 5),  # "
        ]

        flopy.mf6.modflow.mfgwfwel.ModflowGwfwel(
            gwf,
            maxbound=68,
            # auxiliary="IFACE",
            auxiliary=["IFACE", "IFLOWFACE"],
            save_flows=True,
            stress_period_data={0: wells},
        )

        # node property flow
        npf = flopy.mf6.ModflowGwfnpf(
            gwf,
            xt3doptions=True,
            save_flows=True,
            save_specific_discharge=True,
            save_saturation=True,
            icelltype=[0],
            k=[50],
        )

        # constant head boundary (period, node number, head)
        chd_bound = [
            (0, 1327, 150.0),
            (0, 1545, 150.0),
            (0, 1643, 150.0),
            (0, 1687, 150.0),
            (0, 1713, 150.0),
        ]
        chd = flopy.mf6.ModflowGwfchd(
            gwf,
            pname="chd",
            save_flows=True,
            stress_period_data=chd_bound,
            # auxiliary=["IFLOWFACE"]
        )

        # output control
        budget_file = f"{ctx.name}.bud"
        head_file = f"{ctx.name}.hds"
        oc = flopy.mf6.ModflowGwfoc(
            gwf,
            pname="oc",
            budget_filerecord=[budget_file],
            head_filerecord=[head_file],
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        # define particle release points
        particles_prt = [
            (0, (0, 1327), 12500.0, 7062.5, 100.0),
            (1, (0, 1327), 12500.0, 7187.5, 100.0),
            (2, (0, 1327), 12500.0, 7312.5, 100.0),
            (3, (0, 1327), 12500.0, 7437.5, 100.0),
            (4, (0, 1545), 12500.0, 6562.5, 100.0),
            (5, (0, 1545), 12500.0, 6687.5, 100.0),
            (6, (0, 1545), 12500.0, 6812.5, 100.0),
            (7, (0, 1545), 12500.0, 6937.5, 100.0),
            (8, (0, 1643), 12500.0, 6062.5, 100.0),
            (9, (0, 1643), 12500.0, 6187.5, 100.0),
            (10, (0, 1643), 12500.0, 6312.5, 100.0),
            (11, (0, 1643), 12500.0, 6437.5, 100.0),
            (12, (0, 1687), 12500.0, 5562.5, 100.0),
            (13, (0, 1687), 12500.0, 5687.5, 100.0),
            (14, (0, 1687), 12500.0, 5812.5, 100.0),
            (15, (0, 1687), 12500.0, 5937.5, 100.0),
            (16, (0, 1713), 12500.0, 5062.5, 100.0),
            (17, (0, 1713), 12500.0, 5187.5, 100.0),
            (18, (0, 1713), 12500.0, 5312.5, 100.0),
            (19, (0, 1713), 12500.0, 5437.5, 100.0),
            (20, (0, 861), 4500.0, 7820.3125, 100.0),
            (21, (0, 861), 4500.0, 7835.9375, 100.0),
            (22, (0, 861), 4500.0, 7851.5625, 100.0),
            (23, (0, 861), 4500.0, 7867.1875, 100.0),
            (24, (0, 861), 4562.5, 7820.3125, 100.0),
            (25, (0, 861), 4562.5, 7835.9375, 100.0),
            (26, (0, 861), 4562.5, 7851.5625, 100.0),
            (27, (0, 861), 4562.5, 7867.1875, 100.0),
            (28, (0, 861), 4507.8125, 7812.5, 100.0),
            (29, (0, 861), 4523.4375, 7812.5, 100.0),
            (30, (0, 861), 4539.0625, 7812.5, 100.0),
            (31, (0, 861), 4554.6875, 7812.5, 100.0),
            (32, (0, 861), 4507.8125, 7875.0, 100.0),
            (33, (0, 861), 4523.4375, 7875.0, 100.0),
            (34, (0, 861), 4539.0625, 7875.0, 100.0),
            (35, (0, 861), 4554.6875, 7875.0, 100.0),
            (36, (0, 891), 4625.0, 7632.8125, 100.0),
            (37, (0, 891), 4625.0, 7648.4375, 100.0),
            (38, (0, 891), 4625.0, 7664.0625, 100.0),
            (39, (0, 891), 4625.0, 7679.6875, 100.0),
            (40, (0, 891), 4687.5, 7632.8125, 100.0),
            (41, (0, 891), 4687.5, 7648.4375, 100.0),
            (42, (0, 891), 4687.5, 7664.0625, 100.0),
            (43, (0, 891), 4687.5, 7679.6875, 100.0),
            (44, (0, 891), 4632.8125, 7625.0, 100.0),
            (45, (0, 891), 4648.4375, 7625.0, 100.0),
            (46, (0, 891), 4664.0625, 7625.0, 100.0),
            (47, (0, 891), 4679.6875, 7625.0, 100.0),
            (48, (0, 891), 4632.8125, 7687.5, 100.0),
            (49, (0, 891), 4648.4375, 7687.5, 100.0),
            (50, (0, 891), 4664.0625, 7687.5, 100.0),
            (51, (0, 891), 4679.6875, 7687.5, 100.0)
        ]

        # Instantiate the MODFLOW 6 prt model
        prt = flopy.mf6.ModflowPrt(
            sim, modelname=prt_name, model_nam_file="{}.nam".format(prt_name)
        )

        # Instantiate the MODFLOW 6 DISV vertex grid discretization
        disv = flopy.mf6.ModflowGwfdisv(prt, idomain=idomain, **disv_props)

        # Instantiate the MODFLOW 6 prt model input package
        flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=ctx.porosity)

        # Instantiate the MODFLOW 6 prt particle release point (prp) package
        pd = {
            0: ["FIRST"],
        }
        flopy.mf6.ModflowPrtprp(
            prt,
            pname="prp",
            filename="{}_4.prp".format(prt_name),
            nreleasepts=len(particles_prt),
            packagedata=particles_prt,
            perioddata=pd,
        )

        # Instantiate the MODFLOW 6 prt output control package
        budgetfile_prt = "{}.cbb".format(prt_name)
        budget_record = [budgetfile_prt]
        flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            budget_filerecord=budget_record,
            saverecord=[("BUDGET", "ALL")],
        )

        # Instantiate the MODFLOW 6 prt flow model interface
        pd = [
            ("GWFHEAD", headfile_bkwd),
            ("GWFBUDGET", budgetfile_bkwd),
        ]
        flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)

        # create exchange
        flopy.mf6.ModflowGwfprt(
            sim, exgtype="GWF6-PRT6",
            exgmnamea=gwf_name, exgmnameb=prt_name,
            filename=f"{gwf_name}.gwfprt",
        )

        # Create an explicit model solution (EMS) for the MODFLOW 6 prt model
        ems = flopy.mf6.ModflowEms(
            sim,
            pname="ems",
            filename="{}.ems".format(prt_name),
        )
        sim.register_solution_package(ems, [prt.name])

        # todo: add mp7 as comparison simulation
        return ctx, sim, None, self.eval_mp7_p04

    def eval_mp7_p04(self, ctx, sim):
        pass


@parametrize_with_cases("case", cases=PrtCases)
def test_prt_models(case, targets):
    # cases are tuples (context, simulation, optional comparison simulation, evaluation function)
    ctx, sim, cmp, evl = case
    sim.write_simulation()
    if cmp:
        cmp.write_simulation()
    
    test = TestFramework()
    test.run(
        TestSimulation(
            name=ctx.name,
            exe_dict=targets,
            exfunc=lambda s: evl(ctx, s),  # hack the context into the evaluation function for now
            idxsim=0,
            mf6_regression=True,
            require_failure=ctx.xfail,
            make_comparison=False,
        ),
        sim.simulation_data.mfpath.get_sim_path())