"""
Tests a PRT model on the Voronoi grid demonstrated
in Flopy's Voronoi example:

https://flopy.readthedocs.io/en/latest/Notebooks/dis_voronoi_example.html

Two variants are included, first with straight
pathlines then with a well capturing particles.

Particles are released from the x coordinate
of the constant concentration cell from the
transport model, along a range of y coords.

TODO: support parallel adjacent cell faces,
duplicated vertices as flopy.utils.voronoi
can produce via scipy/Qhull (for now flopy
filters these but mf6 probably should too)
"""


from pathlib import Path

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.discretization import VertexGrid
from flopy.utils import GridIntersect
from flopy.utils.triangle import Triangle
from flopy.utils.voronoi import VoronoiGrid
from prt_test_utils import get_model_name
from shapely.geometry import LineString, Point

from framework import TestFramework

simname = "prtfmi08"
cases = [simname, f"{simname}_wel"]
xmin = 0.0
xmax = 2000.0
ymin = 0.0
ymax = 1000.0
top = 1.0
botm = [0.0]
angle_min = 30
area_max = 1000.0
delr = area_max**0.5
nlay = 1
ncol = xmax / delr
nrow = ymax / delr
nodes = ncol * nrow
porosity = 0.1
rpts = [
    [500, 100, 0.5],
    [500, 350, 0.5],
    [500, 450, 0.5],
    [500, 500, 0.5],
    [500, 550, 0.5],
    [500, 650, 0.5],
    [500, 800, 0.5],
]


def get_grid(workspace, targets):
    workspace.mkdir(exist_ok=True, parents=True)
    tri = Triangle(
        maximum_area=area_max,
        angle=angle_min,
        model_ws=workspace,
        exe_name=targets.triangle,
    )
    poly = np.array(((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)))
    tri.add_polygon(poly)
    tri.build(verbose=False)
    return VoronoiGrid(tri)

    # fig = plt.figure(figsize=(10, 10))
    # ax = plt.subplot(1, 1, 1, aspect="equal")
    # pc = tri.plot(ax=ax)


def build_gwf_sim(name, ws, targets):
    ws = Path(ws)
    gwfname = get_model_name(name, "gwf")

    # create grid
    grid = get_grid(ws / "grid", targets)
    vgrid = VertexGrid(**grid.get_gridprops_vertexgrid(), nlay=1)
    ibd = np.zeros(vgrid.ncpl, dtype=int)
    gi = GridIntersect(vgrid)

    # identify cells on left edge
    line = LineString([(xmin, ymin), (xmin, ymax)])
    cells_left = gi.intersect(line)["cellids"]
    cells_left = np.array(list(cells_left))
    ibd[cells_left] = 1

    # identify cells on right edge
    line = LineString([(xmax, ymin), (xmax, ymax)])
    cells_right = gi.intersect(line)["cellids"]
    cells_right = np.array(list(cells_right))
    ibd[cells_right] = 2

    # identify release cell
    point = Point((500, 500))
    cells2 = gi.intersect(point)["cellids"]
    cells2 = np.array(list(cells2))
    # ibd[cells2] = 3

    # identify well cell
    point = Point((1200, 500))
    cell_wel = vgrid.intersect(point.x, point.y)

    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=targets.mf6, sim_ws=ws
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", perioddata=[[1.0, 1, 1.0]]
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="complex",
        outer_dvclose=1.0e-8,
        inner_dvclose=1.0e-8,
    )
    disv = flopy.mf6.ModflowGwfdisv(
        gwf, nlay=nlay, **grid.get_disv_gridprops(), top=top, botm=botm
    )
    if "wel" in name:
        wells = [
            # k, j, q
            (0, cell_wel, -0.5),
        ]
        wel = flopy.mf6.ModflowGwfwel(
            gwf,
            maxbound=len(wells),
            save_flows=True,
            stress_period_data={0: wells},
        )
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=[(True)],
        k=10.0,
        save_saturation=True,
        save_specific_discharge=True,
    )
    ic = flopy.mf6.ModflowGwfic(gwf)

    chdlist = []
    for icpl in cells_left:
        chdlist.append([(0, icpl), 1.0])
    for icpl in cells_right:
        chdlist.append([(0, icpl), 0.0])
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdlist)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return sim


def build_prt_sim(name, gwf_ws, prt_ws, targets):
    prt_ws = Path(prt_ws)
    gwfname = get_model_name(name, "gwf")
    prtname = get_model_name(name, "prt")

    # create grid
    grid = get_grid(prt_ws / "grid", targets)
    gridprops = grid.get_gridprops_vertexgrid()
    vgrid = VertexGrid(**gridprops, nlay=1)
    ibd = np.zeros(vgrid.ncpl, dtype=int)
    gi = GridIntersect(vgrid)

    # identify cells on left edge
    line = LineString([(xmin, ymin), (xmin, ymax)])
    cells0 = gi.intersect(line)["cellids"]
    cells0 = np.array(list(cells0))
    ibd[cells0] = 1

    # identify cells on right edge
    line = LineString([(xmax, ymin), (xmax, ymax)])
    cells1 = gi.intersect(line)["cellids"]
    cells1 = np.array(list(cells1))
    ibd[cells1] = 2

    # identify well cell
    point = Point((800, 500))
    cell_wel = vgrid.intersect(point.x, point.y)

    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=targets.mf6, sim_ws=prt_ws
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", perioddata=[[1.0, 1, 1.0]]
    )
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)
    disv = flopy.mf6.ModflowGwfdisv(
        prt, nlay=nlay, **grid.get_disv_gridprops(), top=top, botm=botm
    )
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    prpdata = [
        # index, (layer, cell index), x, y, z
        (i, (0, vgrid.intersect(p[0], p[1])), p[0], p[1], p[2])
        for i, p in enumerate(rpts)
    ]
    prp_track_file = f"{prtname}.prp.trk"
    prp_track_csv_file = f"{prtname}.prp.trk.csv"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=len(prpdata),
        packagedata=prpdata,
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        boundnames=True,
        stop_at_weak_sink=True,  # currently required for this problem
    )
    prt_track_file = f"{prtname}.trk"
    prt_track_csv_file = f"{prtname}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )
    gwf_budget_file = gwf_ws / f"{gwfname}.bud"
    gwf_head_file = gwf_ws / f"{gwfname}.hds"
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prtname}.ems",
    )
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_models(idx, test):
    gwfsim = build_gwf_sim(test.name, test.workspace, test.targets)
    prtsim = build_prt_sim(
        test.name, test.workspace, test.workspace / "prt", test.targets
    )
    return gwfsim, prtsim


def check_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    gwfname = get_model_name(name, "gwf")
    prtname = get_model_name(name, "prt")
    mp7name = get_model_name(name, "mp7")

    # extract mf6 simulations/models and grid
    gwfsim = test.sims[0]
    prtsim = test.sims[1]

    # get gwf output
    gwf = gwfsim.get_model()
    head = gwf.output.head().get_data()
    bdobj = gwf.output.budget()
    spdis = bdobj.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # get prt output
    prt_track_csv_file = f"{prtname}.prp.trk.csv"
    pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    plot_debug = False
    if plot_debug:
        # plot in 2d with mpl
        fig = plt.figure(figsize=(10, 10))
        ax = plt.subplot(1, 1, 1, aspect="equal")
        pmv = flopy.plot.PlotMapView(model=gwf, ax=ax)
        pmv.plot_grid()
        pmv.plot_array(head, cmap="Blues", alpha=0.25)
        pmv.plot_vector(qx, qy, normalize=True, alpha=0.25)
        if "wel" in name:
            pmv.plot_bc(ftype="WEL")
        mf6_plines = pls.groupby(["iprp", "irpt", "trelease"])
        for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
            pl.plot(
                title=f"MF6 pathlines ({name})",
                kind="line",
                x="x",
                y="y",
                ax=ax,
                legend=False,
                color="black",
            )
        plt.show()

        # plot in 3d with pyvista (via vtk)
        import pyvista as pv
        from flopy.export.vtk import Vtk
        from flopy.plot.plotutil import to_mp7_pathlines

        def get_meshes(model, pathlines):
            vtk = Vtk(model=model, binary=False, smooth=False)
            vtk.add_model(model)
            vtk.add_pathline_points(
                to_mp7_pathlines(pathlines.to_records(index=False))
            )
            grid_mesh, path_mesh = vtk.to_pyvista()
            grid_mesh.rotate_x(-100, point=axes.origin, inplace=True)
            grid_mesh.rotate_z(90, point=axes.origin, inplace=True)
            grid_mesh.rotate_y(120, point=axes.origin, inplace=True)
            path_mesh.rotate_x(-100, point=axes.origin, inplace=True)
            path_mesh.rotate_z(90, point=axes.origin, inplace=True)
            path_mesh.rotate_y(120, point=axes.origin, inplace=True)
            return grid_mesh, path_mesh

        def callback(mesh, value):
            sub = pls[pls.t <= value]
            gm, pm = get_meshes(gwf, sub)
            mesh.shallow_copy(pm)

        pv.set_plot_theme("document")
        axes = pv.Axes(show_actor=True, actor_scale=2.0, line_width=5)
        p = pv.Plotter(notebook=False)
        grid_mesh, path_mesh = get_meshes(gwf, pls)
        p.add_mesh(grid_mesh, scalars=head[0], cmap="Blues", opacity=0.5)
        p.add_mesh(path_mesh, label="Time", style="points", color="black")
        p.camera.zoom(1)
        p.add_slider_widget(lambda v: callback(path_mesh, v), [0, 30202])
        p.show()


@pytest.mark.parametrize("idx, name", enumerate(cases))
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
