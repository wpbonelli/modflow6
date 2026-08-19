"""
Standalone smoke test (not pytest) for the idomain-subset prototype:
GWF-GWT exchange where GWT's active domain is a strict subset of GWF's.

Mirrors smoke_test.py (the PRT version). Three simulations, all sharing a
10x10x1 structured grid and a diagonal flow field (CHD at opposite corners,
same setup as prt_test_utils.FlopyReadmeCase):

1. "baseline": GWT idomain == GWF idomain (all active).
2. "subset":   GWT idomain excludes one GWF-active cell mid-grid. Exercises
   the node map (gwfhead/gwfsat/gwfspdis owned+translated once per fmi_ad),
   the connection map (gwfflowja, used by the SSM dry/rewet-adjacent flow
   correction and advection), and the gwfpackages nodelist reverse-map
   (gwf2loc) used by SSM to source concentration from the CHD boundary.
3. "violation": GWT active in a cell GWF marks inactive -> rejected.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import flopy
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
MF6 = REPO / "builddir" / "src" / "mf6"
WS = Path(__file__).resolve().parent / "_smoke_test_gwt_output"

NLAY, NROW, NCOL = 1, 10, 10
PERLEN, NSTP, NPER = 50.0, 25, 1
POROSITY = 0.1


def build_sim(name, ws, gwt_idomain, gwf_idomain=None):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=str(MF6), version="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)])

    gwfname = f"{name}_gwf"
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    imsgwf = flopy.mf6.ModflowIms(sim, filename=f"{gwfname}.ims")
    sim.register_ims_package(imsgwf, [gwf.name])
    idomain_full = np.ones((NLAY, NROW, NCOL), dtype=int)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=NLAY, nrow=NROW, ncol=NCOL,
        idomain=gwf_idomain if gwf_idomain is not None else idomain_full,
    )
    flopy.mf6.ModflowGwfic(gwf)
    flopy.mf6.ModflowGwfnpf(gwf, save_saturation=True, save_specific_discharge=True)
    spd = {0: [[(0, 0, 0), 1.0, 1.0], [(0, NROW - 1, NCOL - 1), 0.0, 0.0]]}
    flopy.mf6.ModflowGwfchd(
        gwf, pname="CHD-1", stress_period_data=spd, auxiliary=["concentration"]
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwtname = f"{name}_gwt"
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, save_flows=True)
    imsgwt = flopy.mf6.ModflowIms(
        sim, linear_acceleration="BICGSTAB", filename=f"{gwtname}.ims"
    )
    sim.register_ims_package(imsgwt, [gwt.name])
    flopy.mf6.ModflowGwtdis(gwt, nlay=NLAY, nrow=NROW, ncol=NCOL, idomain=gwt_idomain)
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="upstream")
    flopy.mf6.ModflowGwtmst(gwt, porosity=POROSITY)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    ucn_file = f"{gwtname}.ucn"
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=ucn_file,
        saverecord=[("CONCENTRATION", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )
    return sim, ucn_file


def run(ws):
    return subprocess.run([str(MF6)], cwd=str(ws), capture_output=True, text=True)


def main():
    if WS.exists():
        shutil.rmtree(WS)
    WS.mkdir(parents=True)

    idomain_full = np.ones((NLAY, NROW, NCOL), dtype=int)

    # 1. baseline
    ws1 = WS / "baseline"
    sim1, ucn1 = build_sim("base", ws1, idomain_full.copy())
    sim1.write_simulation()
    r1 = run(ws1)
    print("=== baseline ===")
    print(r1.stdout[-1200:])
    ok1 = "Normal termination" in r1.stdout

    # 2. subset: exclude a GWF-active cell mid-grid (row 5, col 9), same
    #    spot used in the PRT smoke test.
    idomain_subset = idomain_full.copy()
    idomain_subset[0, 5, 9] = 0
    ws2 = WS / "subset"
    sim2, ucn2 = build_sim("sub", ws2, idomain_subset, gwf_idomain=idomain_full.copy())
    sim2.write_simulation()
    r2 = run(ws2)
    print("=== subset ===")
    print(r2.stdout[-1200:])
    ok2 = "Normal termination" in r2.stdout

    # 3. violation: GWT active where GWF is inactive -> should error out
    idomain_gwf_holed = idomain_full.copy()
    idomain_gwf_holed[0, 5, 0] = 0
    ws3 = WS / "violation"
    sim3, ucn3 = build_sim("vio", ws3, idomain_full.copy(), gwf_idomain=idomain_gwf_holed)
    sim3.write_simulation()
    r3 = run(ws3)
    print("=== violation (expected to fail) ===")
    print(r3.stdout[-2000:])
    expected_err = "GWT active domain must be a subset" in r3.stdout
    failed_as_expected = "Normal termination" not in r3.stdout

    print("\n\n===== SUMMARY =====")
    print(f"baseline ran to normal termination: {ok1}")
    print(f"subset ran to normal termination:   {ok2}")
    print(f"violation failed as expected:       {failed_as_expected}")
    print(f"violation error message present:    {expected_err}")

    if ok1 and ok2:
        c1 = flopy.utils.HeadFile(ws1 / ucn1, precision="double", text="CONCENTRATION")
        c2 = flopy.utils.HeadFile(ws2 / ucn2, precision="double", text="CONCENTRATION")
        conc1 = c1.get_alldata().reshape(-1, NROW, NCOL)
        conc2 = c2.get_alldata().reshape(-1, NROW, NCOL)
        # Cells far from (5,9) should match baseline (this validates the
        # node/connection maps, same as the PRT smoke test). Note: per the
        # KNOWN ISSUE below, the artifact from dropped boundary flux isn't
        # confined to the immediate neighborhood, so this may still show
        # differences beyond the small mask -- see the design doc.
        mask = np.ones((NROW, NCOL), dtype=bool)
        mask[4:7, 7:10] = False  # exclude a small neighborhood around (5,9)
        far = np.allclose(conc1[:, mask], conc2[:, mask], atol=1e-8)
        print(f"concentrations away from excluded cell match baseline exactly: {far}")
        print(f"final concentration at (5,9) baseline={conc1[-1, 5, 9]:.6g} "
              f"subset={conc2[-1, 5, 9]:.6g} (subset should be the inactive-cell "
              f"fill value 1e30: cell is inactive in GWT)")
        # KNOWN ISSUE (see design doc "Prototype status: GWT"): dropping
        # real GWF flux across the excluded boundary without folding it
        # into GWT's transport equation as a boundary term breaks the
        # per-cell water balance GWT's advection scheme assumes, producing
        # local concentration artifacts near the boundary. Demonstrate it
        # rather than assert correctness here.
        near = conc2[:, ~mask]
        near = near[near < 1e29]  # exclude the inactive-cell fill value
        print(f"KNOWN ISSUE: max concentration near boundary in subset run = "
              f"{near.max():.4g} (source concentration is 1.0 -- values "
              f"meaningfully above 1.0 confirm the water-balance artifact)")


if __name__ == "__main__":
    main()
