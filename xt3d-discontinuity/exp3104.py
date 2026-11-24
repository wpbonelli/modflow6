"""Six-cell 1/4/1 DISU model from #3104, plus a uniform-K control."""

import os
import sys

import flopy

EXE, ROOT = sys.argv[1], sys.argv[2]

vertices = [
    (0, 0, 0),
    (1, 0, 1),
    (2, 1, 1),
    (3, 1, 0),
    (4, 0, 0.5),
    (5, 0.5, 0.5),
    (6, 0.5, 0),
    (7, 1, 0.5),
    (8, 0.5, 1),
]
cell2d = [
    (0, 0.5, 0.5, 4, 0, 1, 2, 3),
    (1, 0.25, 0.25, 4, 0, 4, 5, 6),
    (2, 0.75, 0.25, 4, 6, 5, 7, 3),
    (3, 0.25, 0.75, 4, 4, 1, 8, 5),
    (4, 0.75, 0.75, 4, 5, 8, 2, 7),
    (5, 0.5, 0.5, 4, 0, 1, 2, 3),
]
ja = [
    0,
    1,
    2,
    3,
    4,
    1,
    0,
    2,
    3,
    5,
    2,
    0,
    1,
    4,
    5,
    3,
    0,
    1,
    4,
    5,
    4,
    0,
    2,
    3,
    5,
    5,
    1,
    2,
    3,
    4,
]
angles = [
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    90,
    0,
    0,
    0,
    180,
    90,
    0,
    0,
    0,
    270,
    0,
    0,
    0,
    0,
    270,
    180,
    0,
    0,
    0,
    0,
    0,
    0,
]


def run(name, kmid, xt3d):
    expected = 1 / (1 + 1 / kmid + 1)
    ws = os.path.join(ROOT, name)
    sim = flopy.mf6.MFSimulation(exe_name=EXE, sim_ws=ws)
    sim.simulation_data.float_precision = 17
    sim.simulation_data.float_characters = 25
    flopy.mf6.ModflowTdis(sim, perioddata=[(1, 1, 1)])
    flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="BICGSTAB",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        rcloserecord=1e-12,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdisu(
        gwf,
        nodes=6,
        nja=30,
        nvert=9,
        top=[3, 2, 2, 2, 2, 1],
        bot=[2, 1, 1, 1, 1, 0],
        area=[1, 0.25, 0.25, 0.25, 0.25, 1],
        iac=[5] * 6,
        ja=ja,
        ihc=[0] * 5 + [0, 0, 2, 2, 0] * 4 + [0] * 5,
        cl12=[0, 0.5, 0.5, 0.5, 0.5]
        + [0, 0.5, 0.25, 0.25, 0.5] * 4
        + [0, 0.5, 0.5, 0.5, 0.5],
        hwva=[0, 0.25, 0.25, 0.25, 0.25]
        + [0, 0.25, 0.5, 0.5, 0.25] * 4
        + [0, 0.25, 0.25, 0.25, 0.25],
        angldegx=angles,
        vertices=vertices,
        cell2d=cell2d,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=3.5)
    k = [1, kmid, kmid, kmid, kmid, 1]
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k, k33=k, xt3doptions=xt3d)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[((0,), 4 - 0.5 * expected), ((5,), 3 + 0.5 * expected)]
    )
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord="m.cbc", saverecord=[("BUDGET", "ALL")]
    )
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    assert ok, name
    q = gwf.output.budget().get_data(text="CHD")[0]["q"]
    qd = q[q > 0].sum()
    print(
        f"{name:28s} Q={qd:.9f} expected={expected:.9f} "
        f"err={100 * (qd / expected - 1):+8.2f}%"
    )


for kmid in (0.01, 0.1, 1.0):
    for xt3d in (False, True):
        run(f"k{kmid:g}_{'xt3d' if xt3d else 'std'}", kmid, xt3d)
