"""
Test problem for GWE

One-Dimensional Transport in a Uniform Flow Field.
The purpose of this script is to test the new heat transport model developed
for MODFLOW 6.  To that end, this problem uses the setup of the first MT3DMS
test problem but adapts it for heat. MODFLOW 6 is setup using the new GWE
model with input parameters entered in their native units.

It may be possible to find a 1D heat transport analytical solution in the
future.
"""

# Imports

import os
import numpy as np
import pytest

import flopy

from framework import TestFramework

# Base simulation and model name and workspace

viscosity_on = [False]
cases = ["cnd01"]

# Model units

length_units = "meters"
time_units = "days"

# Table MODFLOW 6 GWE comparison to MT3DMS

nper = 1  # Number of periods
nlay = 1  # Number of layers
ncol = 101  # Number of columns
nrow = 1  # Number of rows
delr = 10.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
top = 0.0  # Top of the model ($m$)
botm = -1.0  # Layer bottom elevations ($m$)
prsity = 0.25  # Porosity
perlen = 2000  # Simulation time ($days$)
k11 = 1.0  # Horizontal hydraulic conductivity ($m/d$)

# Set some static model parameter values

k33 = k11  # Vertical hydraulic conductivity ($m/d$)
laytyp = 1
nstp = 100.0
dt0 = perlen / nstp
Lx = (ncol - 1) * delr
v = 0.24
q = v * prsity
h1 = q * Lx
strt = np.zeros((nlay, nrow, ncol), dtype=float)
strt[0, 0, 0] = h1  # Starting head ($m$)
l = 1000.0  # Needed for plots
icelltype = 1  # Cell conversion type
ibound = np.ones((nlay, nrow, ncol), dtype=int)
ibound[0, 0, 0] = -1
ibound[0, 0, -1] = -1

# Set some static transport related model parameter values

mixelm = 0  # FD
rhob = 1110.0
sp2 = 0.0  # read, but not used in this problem
kd = 1.8168e-4
strt_temp = np.zeros((nlay, nrow, ncol), dtype=float)
dispersivity = 1.0
dmcoef = 3.2519e-7  # Molecular diffusion coefficient

# Set some static heat transport related model parameter values
cpw = 4183.0
rhow = 1000.0
lhv = 2454.0

# Set solver parameter values (and related)
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0
ttsmult = 1.0
dceps = 1.0e-5  # HMOC parameters in case they are invoked
nplane = 1  # HMOC
npl = 0  # HMOC
nph = 4  # HMOC
npmin = 0  # HMOC
npmax = 8  # HMOC
nlsink = nplane  # HMOC
npsink = nph  # HMOC

# Static temporal data used by TDIS file

tdis_rc = []
tdis_rc.append((perlen, nstp, 1.0))

# ### Create MODFLOW 6 GWE MT3DMS Example 1 Boundary Conditions
#
# Constant head cells are specified on both ends of the model

chdspd = [[(0, 0, 0), h1], [(0, 0, ncol - 1), 0.0]]
c0 = 40.0
ctpspd = [[(0, 0, 0), c0]]


def build_models(idx, test):
    # Base MF6 GWE model type
    ws = test.workspace
    name = cases[idx]

    print("Building MF6 model...()".format(name))

    # generate names for each model
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    sim_ws = os.path.join(ws, name)
    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )

    # Instantiating MODFLOW 6 time discretization
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=tdis_rc, time_units=time_units
    )

    # Instantiating MODFLOW 6 groundwater flow model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        model_nam_file="{}.nam".format(gwfname),
    )

    # Instantiating MODFLOW 6 solver for flow model
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwfname),
    )
    sim.register_ims_package(imsgwf, [gwfname])

    # Instantiating MODFLOW 6 discretization package
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
        filename="{}.dis".format(gwfname),
    )

    # Instantiating MODFLOW 6 node-property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=icelltype,
        k=k11,
        k33=k33,
        save_specific_discharge=True,
        filename="{}.npf".format(gwfname),
    )

    # Instantiating MODFLOW 6 initial conditions package for flow model
    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename="{}.ic".format(gwfname))

    # Instantiating VSC
    if viscosity_on[idx]:
        # Instantiate viscosity (VSC) package
        vsc_filerecord = "{}.vsc.bin".format(gwfname)
        vsc_pd = [(0, 0.0, 20.0, gwename, "temperature")]
        flopy.mf6.ModflowGwfvsc(
            gwf,
            viscref=8.904e-4,
            viscosity_filerecord=vsc_filerecord,
            thermal_formulation="nonlinear",
            thermal_a2=10.0,
            thermal_a3=248.37,
            thermal_a4=133.16,
            nviscspecies=len(vsc_pd),
            packagedata=vsc_pd,
            pname="vsc",
            filename="{}.vsc".format(gwfname),
        )

    # Instantiating MODFLOW 6 constant head package
    flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(chdspd),
        stress_period_data=chdspd,
        save_flows=False,
        pname="CHD-1",
        filename="{}.chd".format(gwfname),
    )

    # Instantiating MODFLOW 6 output control package for flow model
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="{}.hds".format(gwfname),
        budget_filerecord="{}.cbc".format(gwfname),
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiating MODFLOW 6 groundwater transport package
    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=gwename,
        model_nam_file="{}.nam".format(gwename),
    )
    gwe.name_file.save_flows = True
    imsgwe = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwename),
    )
    sim.register_ims_package(imsgwe, [gwe.name])

    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename="{}.dis".format(gwename),
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(
        gwe, strt=strt_temp, filename="{}.ic".format(gwename)
    )

    # Instantiating MODFLOW 6 transport advection package
    if mixelm == 0:
        scheme = "UPSTREAM"
    elif mixelm == -1:
        scheme = "TVD"
    else:
        raise Exception()
    flopy.mf6.ModflowGweadv(
        gwe, scheme=scheme, filename="{}.adv".format(gwename)
    )

    # Instantiating MODFLOW 6 transport dispersion package
    if dispersivity != 0:
        flopy.mf6.ModflowGwecnd(
            gwe,
            xt3d_off=True,
            alh=dispersivity,
            ath1=dispersivity,
            ktw=0.5918,
            kts=0.2700,
            filename="{}.cnd".format(gwename),
        )

    # Instantiating MODFLOW 6 transport mass storage package (formerly "reaction" package in MT3DMS)
    flopy.mf6.ModflowGweest(
        gwe,
        save_flows=True,
        porosity=prsity,
        cps=760.0,
        rhos=1500.0,
        packagedata=[cpw, rhow, lhv],
        filename="{}.est".format(gwename),
    )

    # Instantiating MODFLOW 6 transport constant concentration package
    flopy.mf6.ModflowGwectp(
        gwe,
        maxbound=len(ctpspd),
        stress_period_data=ctpspd,
        save_flows=False,
        pname="CTP-1",
        filename="{}.ctp".format(gwename),
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    flopy.mf6.ModflowGwessm(
        gwe, sources=[[]], filename="{}.ssm".format(gwename)
    )

    # Instantiate MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord="{}.cbc".format(gwename),
        temperature_filerecord="{}.ucn".format(gwename),
        temperatureprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiating MODFLOW 6 flow-transport exchange mechanism
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename="{}.gwfgwe".format(name),
    )

    return sim, None


def check_output(idx, test, snapshot):
    print("evaluating results...")
    headfile = flopy.utils.HeadFile(
        test.workspace / f"gwe-{test.name}.ucn", precision="double", text="TEMPERATURE"
    )
    assert snapshot == headfile.get_alldata(), "gwe temperatures don't match expectation"


@pytest.mark.parametrize(
    "idx, name",
    list(enumerate(cases)),
)
def test_mf6model(idx, name, function_tmpdir, targets, array_snapshot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, array_snapshot),
        targets=targets,
    )
    test.run()
