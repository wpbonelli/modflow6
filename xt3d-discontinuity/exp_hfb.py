"""
XT3D + HFB experiments for #3208 and #415.

Based on the 20x20 model from #3208 (Huite): CHD h=1 on the left column,
h=0 on the right column, closed top/bottom, barrier between columns 10 and
11 (1-based). Prints, per run:

  hmax/hmin   head range (should stay within [0, 1])
  net         net flow across the barrier line (col 10 | col 11)
  circ        gross CHD inflow at the right boundary, i.e. water entering
              where it should only leave; a measure of spurious circulation

Usage: python exp_hfb.py <mf6 exe> <workspace root>
"""

import os
import sys

import flopy
import numpy as np

EXE, ROOT = sys.argv[1], sys.argv[2]
N = 20


def run(
    name, k=100.0, k22=1.0, angle=45.0, xt3d=True, hydchr=None, lowk=None, width=1.0
):
    """lowk: replace the HFB by an extra column of this K and width."""
    ws = os.path.join(ROOT, name)
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=EXE)
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="BICGSTAB",
        outer_dvclose=1e-12,
        inner_dvclose=1e-12,
        rcloserecord=[1e-10, "STRICT"],
        inner_maximum=2000,
        outer_maximum=200,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    if lowk is None:
        delr = np.ones(N)
    else:
        delr = np.concatenate([np.ones(10), [width], np.ones(10)])
    ncol = len(delr)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=N, ncol=ncol, delr=delr, delc=1.0, top=1.0, botm=0.0
    )
    kk = np.full((1, N, ncol), k)
    kk22 = np.full((1, N, ncol), k22)
    if lowk is not None:
        kk[:, :, 10] = kk22[:, :, 10] = lowk
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=0, k=kk, k22=kk22, angle1=angle, xt3doptions=xt3d
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.5)
    chd = [[(0, r, 0), 1.0] for r in range(N)]
    chd += [[(0, r, ncol - 1), 0.0] for r in range(N)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    if hydchr is not None:
        flopy.mf6.ModflowGwfhfb(
            gwf,
            maxhfb=N,
            stress_period_data=[[(0, r, 9), (0, r, 10), hydchr] for r in range(N)],
        )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    if not ok:
        print(f"{name:34s} FAILED (see {ws}/mfsim.lst)")
        return

    h = gwf.output.head().get_data()[0]
    cbc = gwf.output.budget()
    fja = cbc.get_data(text="FLOW-JA-FACE")[0].ravel()
    chdq = cbc.get_data(text="CHD")[0]
    grb = flopy.mf6.utils.MfGrdFile(os.path.join(ws, "gwf.dis.grb"))
    ia, ja = grb.ia, grb.ja
    net = 0.0
    for r in range(N):
        n = r * ncol + 9
        for p in range(ia[n] + 1, ia[n + 1]):
            if ja[p] == n + 1:
                net -= fja[p]  # flowja > 0 means into n
    q, node = chdq["q"], chdq["node"] - 1
    right = q[(node % ncol) == ncol - 1]
    circ = right[right > 0].sum()
    print(
        f"{name:34s} hmax={h.max():+8.4f} hmin={h.min():+8.4f} "
        f"net={net:+10.3e} circ={circ:10.3e}"
    )


print("== #3208: remove one ingredient at a time")
run("A_aniso45_hfb", hydchr=1e-12)
run("B_aniso45_hfb_noxt3d", hydchr=1e-12, xt3d=False)
run("C_iso_hfb", k22=100.0, hydchr=1e-12)
run("D_aniso0_hfb", angle=0.0, hydchr=1e-12)
run("E_aniso45_nohfb", hydchr=None)

print("== #3208: hydchr sweep (aniso 45, xt3d)")
for hc in [1e3, 1e1, 1e0, 1e-1, 1e-2, 1e-4, 1e-12]:
    run(f"F_hydchr_{hc:g}", hydchr=hc)

print("== #3208: barrier as a column of low-K cells, no HFB")
for w in [1.0, 0.1, 0.01]:
    run(f"G_lowk_w{w:g}", lowk=1e-12 * w, width=w)
    run(f"G_lowk_w{w:g}_noxt3d", lowk=1e-12 * w, width=w, xt3d=False)
run("H_lowk_iso", k22=100.0, lowk=1e-12)

print("== #415: multiplier option (hydchr <= 0), isotropic, xt3d vs std")
print("   (isotropic regular grid: xt3d should match std exactly)")
for hc in [-0.1, -0.5, -1.0, -1.5, 0.0]:
    for x in (False, True):
        tag = "xt3d" if x else "std"
        run(f"I_mult{hc:g}_{tag}", k22=100.0, angle=0.0, hydchr=hc, xt3d=x)
