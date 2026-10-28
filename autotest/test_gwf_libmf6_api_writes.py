"""
Verify that values written through the API persist and reach the solution in
the advanced GWF packages.

IDM integration moved time series resolution for PACKAGEDATA and OUTLETS
values into the input context. The packages copy those values into their own
arrays each time step, which overwrote values written through the API
(they persisted in 6.8.1):

    case            package  written (old name)  overwritten when
    maw             MAW      WELL_HEAD           always (constant-head well)
    sfr / sfr_ts    SFR      HK                  a TS6 file is present
    lak / lak_ts    LAK      OUTWIDTH            a TS6 file is present

Each case writes once and steps through two stress periods; the second
PERIOD block sets a different entry. The written value must survive every
time step, and the solution must match a reference run whose input file sets
the same value. The write happens at one of three points in the first time
step, to show which BMI call overwrites it (the packages copy input values in
their *_ad routines, which prepare_solve calls):

    phase   written after
    init    initialize
    pts     prepare_time_step
    ps      prepare_solve
"""

import flopy
import numpy as np
import pytest
from modflow_devtools.markers import requires_pkg

cases = ["maw", "sfr", "sfr_ts", "lak", "lak_ts"]
phases = ["init", "pts", "ps"]

name = "api_writes"
nper, perlen, nstp = 2, 4.0, 2
top, botm = 10.0, 0.0


def get_sim(ws, case, reference):
    """Build the simulation. For the reference run, the input file sets the
    value the API run writes."""
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(perlen, nstp, 1.0)] * nper)
    flopy.mf6.ModflowIms(
        sim,
        complexity="complex",
        outer_dvclose=1e-9,
        inner_dvclose=1e-9,
        outer_maximum=200,
        inner_maximum=200,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    ncol = 5 if case == "maw" else 3
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=ncol, delr=100.0, delc=100.0, top=top, botm=botm
    )
    k = 10.0 if case == "maw" else 1.0
    strt = 9.0 if case == "maw" else 8.0
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
    ts = case.endswith("_ts")

    if case == "maw":
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data={0: [[(0, 0, 0), strt], [(0, 0, ncol - 1), strt]]},
        )
        # -- well 2 is constant head, with no PERIOD WELL_HEAD in the API run
        pd = [[n, 0.1, botm, strt, "THIEM", 1] for n in range(2)]
        cd = [[n, 0, (0, 0, 1 + 2 * n), top, botm, -999.0, -999.0] for n in range(2)]
        per = [[0, "RATE", -50.0], [1, "STATUS", "CONSTANT"]]
        if reference:
            per.append([1, "WELL_HEAD", 7.0])
        flopy.mf6.ModflowGwfmaw(
            gwf,
            nmawwells=2,
            packagedata=pd,
            connectiondata=cd,
            perioddata={0: per, 1: [[0, "RATE", -50.0]]},
            pname="MAW-1",
        )

    elif case.startswith("sfr"):
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 0), strt]]})
        # -- reach 2 bed conductivity is written; the TS6 file drives reach 1
        rhk = [0.1, 0.0 if reference else 0.1, 0.1]
        pd = [
            [i, (0, 0, i), 100.0, 5.0, 0.001, 9.0, 1.0, rhk[i], 0.03, nc, 1.0, 0]
            for i, nc in enumerate([1, 2, 1])
        ]
        cd = [[0, -1], [1, 0, -2], [2, 1]]
        per = [[0, "INFLOW", "qin" if ts else 10.0]]
        sfr = flopy.mf6.ModflowGwfsfr(
            gwf,
            nreaches=3,
            packagedata=pd,
            connectiondata=cd,
            perioddata={0: per, 1: per},
            pname="SFR-1",
        )
        if ts:
            sfr.ts.initialize(
                filename=f"{name}.sfr.ts",
                timeseries=[(0.0, 10.0), (nper * perlen, 10.0)],
                time_series_namerecord="qin",
                interpolation_methodrecord="linear",
            )

    elif case.startswith("lak"):
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 2), strt]]})
        # -- outlet width is written; the TS6 file drives rainfall
        pd = [[0, 11.0, 1]]
        cd = [[0, 0, (0, 0, 0), "VERTICAL", 0.001, 0.0, 0.0, 0.0, 0.0]]
        width = 10.0 if reference else 2.0
        outlets = [[0, 0, -1, "WEIR", 10.5, width, 0.03, 0.001]]
        per = [[0, "RAINFALL", "rain" if ts else 0.5]]
        lak = flopy.mf6.ModflowGwflak(
            gwf,
            nlakes=1,
            noutlets=1,
            packagedata=pd,
            connectiondata=cd,
            outlets=outlets,
            perioddata={0: per, 1: per},
            pname="LAK-1",
        )
        if ts:
            lak.ts.initialize(
                filename=f"{name}.lak.ts",
                timeseries=[(0.0, 0.5), (nper * perlen, 0.5)],
                time_series_namerecord="rain",
                interpolation_methodrecord="linear",
            )
    return sim


# -- case: (package, written old name, index, value, package state to compare)
writes = {
    "maw": ("MAW-1", "WELL_HEAD", 1, 7.0, "HEAD"),
    "sfr": ("SFR-1", "HK", 1, 0.0, "DSFLOW"),
    "lak": ("LAK-1", "OUTWIDTH", 0, 10.0, "SIMOUTRATE"),
}


def run_sim(lib, ws, case, phase=None):
    """Run through the API, writing once at the given phase of the first time
    step (no write if phase is None). Returns per-step model heads and package
    state, and the written variable after each step."""
    from xmipy import XmiWrapper

    pname, var, idx, value, state = writes[case.removesuffix("_ts")]
    prefix = f"{name.upper()}/{pname}"
    mf6 = XmiWrapper(lib, working_directory=str(ws))
    mf6.initialize()
    kstp = 0

    def write(at):
        if at == phase and kstp == 0:
            mf6.get_value_ptr(mf6.get_var_address(var, name.upper(), pname))[idx] = (
                value
            )

    try:
        write("init")
        heads, states, written = [], [], []
        while mf6.get_current_time() < mf6.get_end_time():
            mf6.prepare_time_step(mf6.get_time_step())
            write("pts")
            nsol = mf6.get_subcomponent_count()
            for sol_id in range(1, nsol + 1):
                mf6.prepare_solve(sol_id)
                if sol_id == nsol:
                    write("ps")
                for _ in range(mf6.get_value("SLN_1/MXITER")[0]):
                    if mf6.solve(sol_id):
                        break
                else:
                    raise RuntimeError(f"no convergence at step {kstp + 1}")
                mf6.finalize_solve(sol_id)
            mf6.finalize_time_step()
            kstp += 1
            heads.append(mf6.get_value(f"{name.upper()}/X"))
            states.append(mf6.get_value(f"{prefix}/{state}"))
            written.append(mf6.get_value(f"{prefix}/{var}")[idx])
    finally:
        mf6.finalize()
    return np.array(heads), np.array(states), np.array(written)


def check(lib, ws, case, phase):
    api_ws, ref_ws = ws / "api", ws / "reference"
    get_sim(api_ws, case, reference=False).write_simulation(silent=True)
    get_sim(ref_ws, case, reference=True).write_simulation(silent=True)

    _, var, _, value, state = writes[case.removesuffix("_ts")]
    heads, states, written = run_sim(lib, api_ws, case, phase)
    ref_heads, ref_states, _ = run_sim(lib, ref_ws, case)

    assert np.allclose(written, value), (
        f"{var} written after {phase} was overwritten: "
        f"{written} after each step (wrote {value})"
    )
    assert np.allclose(states, ref_states), (
        f"{state} differs from the reference run: {states} vs {ref_states}"
    )
    assert np.allclose(heads, ref_heads), (
        f"heads differ from the reference run: {heads} vs {ref_heads}"
    )


@requires_pkg("xmipy")
@pytest.mark.parametrize("phase", phases)
@pytest.mark.parametrize("case", cases)
def test_api_writes(case, phase, function_tmpdir, targets):
    check(targets["libmf6"], function_tmpdir, case, phase)
