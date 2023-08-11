from pathlib import Path
import shutil
import numpy as np
import flopy
from flopy.mf6 import MFSimulation


sim_name = "regional"
gwf_ws = Path("temp").absolute() / sim_name / "gwf"

# load gwf simulation
sim = MFSimulation.load(sim_name, sim_ws=gwf_ws)
gwf = sim.get_model()

# extract discretization props
delr = gwf.dis.delr.array[0]
delc = gwf.dis.delc.array[0]
nlay = gwf.dis.nlay.data
nrow = gwf.dis.nrow.data
ncol = gwf.dis.ncol.data
idomain = gwf.dis.idomain.array[0].copy()
tdis = sim.tdis

# set grid spatial ref info
gwf.modelgrid.set_coord_info(xoff= 518200.0,yoff=350400.0, epsg=3070)

# starting points: let's initially roll with a particle in each active cell
def get_release_pts():
    sloc_tmpl_size = 1
    release_cells = []
    ndivc = sloc_tmpl_size
    ndivr = sloc_tmpl_size
    deldivc = delc / ndivc
    deldivr = delr / ndivr
    k = 0
    top = gwf.dis.top.array
    nrpt = -1
    for i in range(nrow):
        y0 = (nrow - i - 1) * delc
        for j in range(ncol):
            x0 = j * delr
            if idomain[i,j] != 0:
                zrpt = top[i,j]
                for idiv in range(ndivc):
                    dy = (idiv + 0.5) * deldivc
                    yrpt = y0 + dy
                    for jdiv in range(ndivr):
                        dx = (jdiv + 0.5) * deldivr
                        xrpt = x0 + dx
                        nrpt += 1
                        rpt = [nrpt, k, i, j, xrpt, yrpt, zrpt]
                        release_cells.append(rpt)
    return [[i] + row[1:] for i, row  in enumerate(release_cells)]

# get release points 
release_pts = get_release_pts() 

# create prt workspace
prt_ws = gwf_ws.parent / "prt"
prt_ws.mkdir(exist_ok=True, parents=True)

# create prt simulation (load discretizations from gwf sim)
sim = MFSimulation.load(sim_name, sim_ws=gwf_ws, load_only=["tdis"])
sim.remove_model(sim_name)

# change to prt workspace
sim.set_sim_path(prt_ws)

# add prt model
prt_name = f'{sim_name}_prt'
prt = flopy.mf6.ModflowPrt(
    sim,
    modelname=prt_name,
    model_nam_file=f'{prt_name}.nam'
)

# discretization
flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
    prt,
    nlay=nlay,
    nrow=nrow,
    ncol=ncol,
    length_units="FEET",
    delr=gwf.dis.delr.array,
    delc=gwf.dis.delc.array,
    top=gwf.dis.top.array,
    botm=gwf.dis.botm.array,
    idomain=gwf.dis.idomain.array,
    filename=f'{prt_name}.dis'
)

# model input package
flopy.mf6.ModflowPrtmip(
    prt,
    porosity=0.01
)

# particle release points
flopy.mf6.ModflowPrtprp(
    prt,
    filename=f'{prt_name}.prp',
    nreleasepts=len(release_pts), 
    packagedata=release_pts,
    perioddata={0: ["FIRST"],})

# output control
flopy.mf6.ModflowPrtoc(
    prt,
    budget_filerecord=[f'{prt_name}.cbb'],
    trackcsv_filerecord=f'{prt_name}.csv',
    saverecord=[("BUDGET", "ALL")],
    filename=f'{sim_name}_prt.oc'
)

# flow model interface
flopy.mf6.ModflowPrtfmi(
    prt,
    packagedata=[
        ['GWFBUDGET', gwf_ws / f'{sim_name}.cbc'],
        ['GWFHEAD', gwf_ws / f'{sim_name}.hds']
    ]
)

# explicit model solution
ems = flopy.mf6.ModflowEms(
    sim,
    filename=f'{prt_name}.ems',
)
sim.register_solution_package(ems, [prt.name])

# write prt simulation
sim.write_simulation()

# kludge to fully purge the gwf model from the simulation
# todo: fix this in flopy
# https://github.com/modflowpy/flopy/issues/1916
n = open(prt_ws / 'mfsim.nam', 'r').readlines()
with open(prt_ws / 'mfsim.nam', 'w') as ofp:
    for line in n:
        if 'ims' not in line.lower() and 'gwf' not in line.lower():
            ofp.write(line)

# run simulation
# sim.run_simulation()