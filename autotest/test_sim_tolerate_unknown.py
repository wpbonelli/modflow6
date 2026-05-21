"""
Test that the TOLERATE_UNKNOWN simulation option causes unrecognized
top-level input tags to produce a warning rather than an error.

Two cases:
  1. Without the option: an unknown tag in a package OPTIONS block is an error.
  2. With the option:    the same unknown tag produces a warning and the
                         simulation completes normally.
"""

import subprocess

import flopy


def run_mf6(exe, ws):
    proc = subprocess.Popen(
        [exe], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ws
    )
    result, _ = proc.communicate()
    buff = result.decode("utf-8").splitlines() if result else []
    return proc.returncode, buff


def build_sim(ws, exe, tolerate_unknown=False):
    """Build and write a minimal GWF simulation."""
    name = "tolunknown"
    simkwargs = {"tolerate_unknown": True} if tolerate_unknown else {}
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=exe, sim_ws=ws, **simkwargs
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowIms(sim, print_option="SUMMARY")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, top=0.0, botm=[-1.0])
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(0, 0, 0), 0.0], [(0, 2, 2), 1.0]]}
    )
    sim.write_simulation()
    return name


def inject_bogus_option(ws, name):
    """Prepend an unrecognized keyword to the NPF OPTIONS block."""
    npf_file = ws / f"{name}.npf"
    with open(npf_file, "w") as f:
        f.write("BEGIN options\n")
        f.write("  BOGUS_OPTION\n")
        f.write("END options\n\n")
        f.write("BEGIN griddata\n")
        f.write("  ICELLTYPE\n")
        f.write("    CONSTANT 0\n")
        f.write("  K\n")
        f.write("    CONSTANT 1.0\n")
        f.write("END griddata\n")


def test_tolerate_unknown_error(function_tmpdir, targets):
    """Without TOLERATE_UNKNOWN, an unrecognized option tag is a fatal error."""
    mf6 = targets["mf6"]
    name = build_sim(function_tmpdir, mf6, tolerate_unknown=False)
    inject_bogus_option(function_tmpdir, name)

    returncode, _ = run_mf6(mf6, function_tmpdir)
    assert returncode != 0, "mf6 should have failed with an unknown input tag"


def test_tolerate_unknown_warning(function_tmpdir, targets):
    """With TOLERATE_UNKNOWN, an unrecognized option tag is a warning and the
    simulation completes normally."""
    mf6 = targets["mf6"]
    name = build_sim(function_tmpdir, mf6, tolerate_unknown=True)
    inject_bogus_option(function_tmpdir, name)

    returncode, buff = run_mf6(mf6, function_tmpdir)
    assert returncode == 0, "mf6 failed unexpectedly:\n" + "\n".join(buff)

    lst = (function_tmpdir / "mfsim.lst").read_text()
    assert "BOGUS_OPTION" in lst, "expected unknown tag name in listing file warning"
    assert "ignored" in lst.lower(), "expected 'ignored' in listing file warning"
    assert "Normal termination" in lst, "expected normal termination in listing file"
