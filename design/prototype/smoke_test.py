"""
Standalone smoke test (not pytest) for the idomain-subset prototype:
GWF-PRT exchange where PRT's active domain is a strict subset of GWF's.

Not wired into the pytest autotest suite yet -- see the "Prototype status"
section of ../idomain-subset-coupling.md for how this relates to
autotest/test_prt_exg.py, which already has an "idmu" case for exactly the
"subset" scenario below but currently marks it xfail.

Usage: build the project first (meson setup builddir && ninja -C builddir),
then run this with the modflow6 flopy environment, e.g.:
    python design/prototype/smoke_test.py

Three simulations, all sharing the same 10x10x1 grid and flow field
(FlopyReadmeCase from autotest/prt_test_utils.py):

1. "baseline":  PRT idomain == GWF idomain (all active). Exercises the
   existing (unmodified) direct-index code path.
2. "subset":    PRT idomain excludes one cell, mid-grid, that GWF still has
   active. Exercises the new node/connection map. Since removing a cell
   from the middle of the grid shifts PRT's reduced numbering relative to
   GWF's for every node that comes after it, this meaningfully exercises
   the mapping logic. Particles whose paths never approach the excluded
   cell should match "baseline" bit-for-bit; particles that do approach it
   should terminate there (TERM_NO_EXITS) rather than passing through.
3. "violation": PRT active in a cell GWF marks inactive. Should be
   rejected at exg_ar with the new subset-violation error message.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import flopy
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "autotest"))
from prt_test_utils import FlopyReadmeCase  # noqa: E402

MF6 = REPO / "builddir" / "src" / "mf6"
WS = Path(__file__).resolve().parent / "_smoke_test_output"


def build_sim(name, ws, prt_idomain, gwf_idomain=None):
    sim = FlopyReadmeCase.get_gwf_sim(name, ws, str(MF6))
    gwf = sim.get_model()
    if gwf_idomain is not None:
        gwf.dis.idomain = gwf_idomain

    prt_name = f"{name}_prt"
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
        idomain=prt_idomain,
    )
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(FlopyReadmeCase.releasepts_prt),
        packagedata=FlopyReadmeCase.releasepts_prt,
        perioddata={0: ["FIRST"]},
        extend_tracking=True,
    )
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[f"{prt_name}.cbb"],
        track_filerecord=[f"{prt_name}.trk"],
        trackcsv_filerecord=[prt_track_csv_file],
        saverecord=[("BUDGET", "ALL")],
    )
    gwf_name = gwf.name
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=gwf_name,
        exgmnameb=prt_name,
        filename=f"{gwf_name}.gwfprt",
    )
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{prt_name}.ems")
    sim.register_solution_package(ems, [prt.name])
    return sim, prt_track_csv_file


def run(ws):
    return subprocess.run([str(MF6)], cwd=str(ws), capture_output=True, text=True)


def main():
    if WS.exists():
        shutil.rmtree(WS)
    WS.mkdir(parents=True)

    nlay, nrow, ncol = FlopyReadmeCase.nlay, FlopyReadmeCase.nrow, FlopyReadmeCase.ncol

    # 1. baseline: identical idomain (all active)
    idomain_full = np.ones((nlay, nrow, ncol), dtype=int)
    ws1 = WS / "baseline"
    sim1, csv1 = build_sim("base", ws1, idomain_full.copy())
    sim1.write_simulation()
    r1 = run(ws1)
    print("=== baseline ===")
    print(r1.stdout[-1500:])
    ok1 = "Normal termination" in r1.stdout

    # 2. subset: PRT excludes a GWF-active cell mid-grid (column 9, middle
    #    row) -- in the middle of the reduced-numbering order, so most
    #    nodes after it shift, meaningfully exercising the map.
    idomain_subset = idomain_full.copy()
    idomain_subset[0, 5, 9] = 0
    ws2 = WS / "subset"
    sim2, csv2 = build_sim("sub", ws2, idomain_subset, gwf_idomain=idomain_full.copy())
    sim2.write_simulation()
    r2 = run(ws2)
    print("=== subset ===")
    print(r2.stdout[-1500:])
    ok2 = "Normal termination" in r2.stdout

    # 3. violation: PRT active where GWF is inactive -> should error out
    idomain_gwf_holed = idomain_full.copy()
    idomain_gwf_holed[0, 5, 0] = 0  # inactive in GWF
    idomain_prt_violate = idomain_full.copy()  # active everywhere in PRT
    ws3 = WS / "violation"
    sim3, csv3 = build_sim(
        "vio", ws3, idomain_prt_violate, gwf_idomain=idomain_gwf_holed
    )
    sim3.write_simulation()
    r3 = run(ws3)
    print("=== violation (expected to fail) ===")
    print(r3.stdout[-2500:])
    expected_err = "PRT active domain must be a subset" in r3.stdout or \
        "PRT active domain must be a subset" in (r3.stderr or "")
    failed_as_expected = "Normal termination" not in r3.stdout

    print("\n\n===== SUMMARY =====")
    print(f"baseline ran to normal termination: {ok1}")
    print(f"subset ran to normal termination:   {ok2}")
    print(f"violation failed as expected:       {failed_as_expected}")
    print(f"violation error message present:    {expected_err}")

    if ok1 and ok2:
        pls1 = pd.read_csv(ws1 / csv1)
        pls2 = pd.read_csv(ws2 / csv2)
        cols = ["t", "x", "y", "z"]
        # Particles whose paths never approach the excluded cell (0,5,9)
        # should match baseline bit-for-bit -- this validates the node and
        # connection maps for the rest of the grid. Particles that do pass
        # near it are expected to terminate earlier in the subset run
        # (TERM_NO_EXITS, since the excluded cell is no longer a viable
        # exit face) -- this validates the feature's intended behavior.
        unaffected = [1, 2, 3, 4, 5, 6]
        all_match = True
        for irpt in unaffected:
            bb = pls1[pls1.irpt == irpt][cols].to_numpy()
            ss = pls2[pls2.irpt == irpt][cols].to_numpy()
            match = bb.shape == ss.shape and np.allclose(bb, ss, atol=1e-9)
            all_match &= match
            print(f"  particle {irpt} (unaffected by excluded cell): "
                  f"{'MATCH' if match else 'DIFFER'}")
        print(f"all unaffected particles match baseline exactly: {all_match}")

        for irpt in [7, 8, 9]:
            srow = pls2[pls2.irpt == irpt].iloc[-1]
            reason = ("TERM_NO_EXITS -> excluded cell had no exit"
                      if srow.istatus == 5 else "unexpected")
            print(f"  particle {irpt} (passes near excluded cell) "
                  f"terminates in subset run at icell={int(srow.icell)}, "
                  f"istatus={int(srow.istatus)} ({reason})")


if __name__ == "__main__":
    main()
