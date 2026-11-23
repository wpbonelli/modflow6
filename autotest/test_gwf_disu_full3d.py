"""
Test the FULL_3D option of the DISU Package, which supports fully
three-dimensional unstructured grids (for example, 3D Voronoi grids) in which
cells have no top and bottom elevation.

A structured box (3 layers, 4 rows, 5 columns, irregular cell sizes) is
simulated twice, with a pair of constant-head faces and two wells:

  1. as a classic multi-layer DISU grid (TOP, BOT, AREA, IHC = 0 or 1), and
  2. as a FULL_3D DISU grid (VOLUME, XCELL, YCELL, ZCELL, IHC = 3, HWVA is
     the face area, CL12 is the distance from the cell center to the face).

Heads must be identical for steady-state and transient conditions, with and
without anisotropy (K11, K22, K33 aligned with x, y and z).  The remaining
cases check that unsupported combinations are rejected with an error.
"""

import subprocess
from pathlib import Path

import flopy
import numpy as np
import pytest

nlay, nrow, ncol = 3, 4, 5
N = nlay * nrow * ncol


def nid(l, r, c):
    return (l * nrow + r) * ncol + c


def _array(name, values, fmt="%.15g"):
    values = np.asarray(values).ravel()
    text = f"  {name}\n    INTERNAL\n"
    for i in range(0, values.size, 6):
        text += "      " + " ".join(fmt % v for v in values[i : i + 6]) + "\n"
    return text


def build_grid():
    rng = np.random.default_rng(1)
    dx = rng.uniform(0.8, 1.6, ncol)
    dy = rng.uniform(0.7, 1.4, nrow)
    dz = rng.uniform(0.5, 1.5, nlay)
    xe = np.r_[0, np.cumsum(dx)]
    ye = np.r_[0, np.cumsum(dy)]
    ze = np.r_[0, -np.cumsum(dz)]
    g = {
        k: np.zeros(N)
        for k in ("top", "bot", "area", "xc", "yc", "zc", "vol", "dx", "dy")
    }
    nb = [[] for _ in range(N)]
    for l in range(nlay):
        for r in range(nrow):
            for c in range(ncol):
                n = nid(l, r, c)
                g["top"][n], g["bot"][n] = ze[l], ze[l + 1]
                g["area"][n] = dx[c] * dy[r]
                g["dx"][n], g["dy"][n] = dx[c], dy[r]
                g["xc"][n] = 0.5 * (xe[c] + xe[c + 1])
                g["yc"][n] = 0.5 * (ye[r] + ye[r + 1])
                g["zc"][n] = 0.5 * (ze[l] + ze[l + 1])
                g["vol"][n] = g["area"][n] * dz[l]
                for dl, dr, dc in [
                    (-1, 0, 0),
                    (0, -1, 0),
                    (0, 0, -1),
                    (0, 0, 1),
                    (0, 1, 0),
                    (1, 0, 0),
                ]:
                    L, R, C = l + dl, r + dr, c + dc
                    if not (0 <= L < nlay and 0 <= R < nrow and 0 <= C < ncol):
                        continue
                    # m, ihc, cl12, hwva (classic), angldegx, face area
                    if dl:
                        t = (
                            nid(L, R, C),
                            0,
                            0.5 * dz[l],
                            dx[c] * dy[r],
                            0.0,
                            dx[c] * dy[r],
                        )
                    elif dc:
                        t = (
                            nid(L, R, C),
                            1,
                            0.5 * dx[c],
                            dy[r],
                            0.0 if dc > 0 else 180.0,
                            dy[r] * dz[l],
                        )
                    else:
                        t = (
                            nid(L, R, C),
                            1,
                            0.5 * dy[r],
                            dx[c],
                            90.0 if dr > 0 else 270.0,
                            dx[c] * dz[l],
                        )
                    nb[n].append(t)
    nb = [sorted(x) for x in nb]
    return g, nb


def write_disu(path, g, nb, full3d, bad_ihc=False):
    iac = [len(x) + 1 for x in nb]
    ja, ihc, cl12, hwva, ang = [], [], [], [], []
    for n, x in enumerate(nb):
        (
            ja.append(n + 1),
            ihc.append(0),
            cl12.append(0.0),
            hwva.append(0.0),
            ang.append(0.0),
        )
        for k, (m, ih, cl, hw, an, fa) in enumerate(x):
            ja.append(m + 1)
            ihc.append(3 if full3d else ih)
            cl12.append(cl)
            hwva.append(fa if full3d else hw)
            ang.append(an)
    if bad_ihc:
        ihc[1] = 0
    text = "BEGIN OPTIONS\n" + ("  FULL_3D\n" if full3d else "") + "END OPTIONS\n\n"
    text += f"BEGIN DIMENSIONS\n  NODES {N}\n  NJA {len(ja)}\nEND DIMENSIONS\n\n"
    text += "BEGIN GRIDDATA\n"
    if full3d:
        text += _array("VOLUME", g["vol"]) + _array("XCELL", g["xc"])
        text += _array("YCELL", g["yc"]) + _array("ZCELL", g["zc"])
    else:
        text += (
            _array("TOP", g["top"])
            + _array("BOT", g["bot"])
            + _array("AREA", g["area"])
        )
    text += "END GRIDDATA\n\nBEGIN CONNECTIONDATA\n"
    text += _array("IAC", iac, "%d") + _array("JA", ja, "%d") + _array("IHC", ihc, "%d")
    text += _array("CL12", cl12) + _array("HWVA", hwva)
    if not full3d:
        text += _array("ANGLDEGX", ang)
    text += "END CONNECTIONDATA\n"
    Path(path).write_text(text)


def write_model(
    ws,
    g,
    nb,
    full3d,
    aniso,
    transient,
    npf_extra="",
    sto_extra="",
    bad_ihc=False,
    accel="CG",
    chd=None,
    wel=None,
    npf_grid="",
):
    ws = Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    write_disu(ws / "m.disu", g, nb, full3d, bad_ihc)
    (ws / "mfsim.nam").write_text(
        "BEGIN TIMING\n  TDIS6 m.tdis\nEND TIMING\n"
        "BEGIN MODELS\n  GWF6 m.nam m\nEND MODELS\n"
        "BEGIN EXCHANGES\nEND EXCHANGES\n"
        "BEGIN SOLUTIONGROUP 1\n  IMS6 m.ims m\nEND SOLUTIONGROUP\n"
    )
    nstp = 6 if transient else 1
    (ws / "m.tdis").write_text(
        "BEGIN DIMENSIONS\n  NPER 1\nEND DIMENSIONS\n"
        f"BEGIN PERIODDATA\n  5.0 {nstp} 1.0\nEND PERIODDATA\n"
    )
    (ws / "m.ims").write_text(
        "BEGIN NONLINEAR\n  OUTER_DVCLOSE 1e-9\nEND NONLINEAR\n"
        "BEGIN LINEAR\n  INNER_MAXIMUM 500\n  INNER_DVCLOSE 1e-11\n"
        "  INNER_RCLOSE 1e-9\n"
        f"  LINEAR_ACCELERATION {accel}\nEND LINEAR\n"
    )
    (ws / "m.nam").write_text(
        "BEGIN PACKAGES\n  DISU6 m.disu\n  NPF6 m.npf\n  IC6 m.ic\n  OC6 m.oc\n"
        "  STO6 m.sto\n"
        "  CHD6 m.chd\n  WEL6 m.wel\nEND PACKAGES\n"
    )
    npf = f"BEGIN OPTIONS\n{npf_extra}END OPTIONS\nBEGIN GRIDDATA\n"
    npf += _array("ICELLTYPE", np.zeros(N), "%d") + _array("K", np.full(N, 3.0))
    if aniso:
        npf += _array("K22", np.full(N, 1.0)) + _array("K33", np.full(N, 0.4))
        npf += _array("ANGLE1", np.zeros(N))
    (ws / "m.npf").write_text(npf + npf_grid + "END GRIDDATA\n")
    (ws / "m.ic").write_text(
        "BEGIN GRIDDATA\n" + _array("STRT", np.zeros(N)) + "END GRIDDATA\n"
    )
    (ws / "m.sto").write_text(
        f"BEGIN OPTIONS\n{sto_extra}END OPTIONS\nBEGIN GRIDDATA\n"
        + _array("ICONVERT", np.zeros(N), "%d")
        + _array("SS", np.full(N, 2e-3 if transient else 0.0))
        + "END GRIDDATA\n"
        + f"BEGIN PERIOD 1\n  {'TRANSIENT' if transient else 'STEADY-STATE'}\n"
        + "END PERIOD\n"
    )
    (ws / "m.oc").write_text(
        "BEGIN OPTIONS\n  HEAD FILEOUT m.hds\nEND OPTIONS\n"
        "BEGIN PERIOD 1\n  SAVE HEAD ALL\nEND PERIOD\n"
    )
    if chd is None:
        chd = [(nid(l, r, 0), 10.0) for l in range(nlay) for r in range(nrow)]
        chd += [(nid(l, r, ncol - 1), 2.0) for l in range(nlay) for r in range(nrow)]
    (ws / "m.chd").write_text(
        f"BEGIN DIMENSIONS\n  MAXBOUND {len(chd)}\nEND DIMENSIONS\nBEGIN PERIOD 1\n"
        + "".join(f"  {c + 1} {h}\n" for c, h in chd)
        + "END PERIOD\n"
    )
    if wel is None:
        wel = [(nid(1, 2, 2), -0.7), (nid(2, 1, 3), 0.4)]
    (ws / "m.wel").write_text(
        f"BEGIN DIMENSIONS\n  MAXBOUND {len(wel)}\nEND DIMENSIONS\nBEGIN PERIOD 1\n"
        + "".join(f"  {c + 1} {q}\n" for c, q in wel)
        + "END PERIOD\n"
    )


def run_mf6(mf6, ws):
    return subprocess.run([str(mf6)], cwd=ws, capture_output=True, text=True)


def heads(ws):
    hf = flopy.utils.HeadFile(Path(ws) / "m.hds", precision="double")
    return np.array([hf.get_data(totim=t).ravel() for t in hf.get_times()])


@pytest.mark.parametrize("aniso", [False, True], ids=["iso", "aniso"])
@pytest.mark.parametrize("transient", [False, True], ids=["steady", "transient"])
def test_full3d_matches_classic_disu(function_tmpdir, targets, aniso, transient):
    g, nb = build_grid()
    out = {}
    for full3d in (False, True):
        ws = function_tmpdir / ("full3d" if full3d else "classic")
        write_model(ws, g, nb, full3d, aniso, transient)
        r = run_mf6(targets["mf6"], ws)
        assert r.returncode == 0, r.stdout[-2000:]
        out[full3d] = heads(ws)
    assert out[True].shape == out[False].shape
    np.testing.assert_allclose(out[True], out[False], rtol=0, atol=1e-8)
    # the solution is not trivial
    assert out[True].max() - out[True].min() > 5.0


@pytest.mark.parametrize(
    "case, kwargs, message",
    [
        (
            "stocoef",
            dict(sto_extra="  STORAGECOEFFICIENT\n"),
            "STORAGECOEFFICIENT option is not supported",
        ),
        ("badihc", dict(bad_ihc=True), "IHC must be 3 for all connections"),
    ],
)
def test_full3d_unsupported(function_tmpdir, targets, case, kwargs, message):
    g, nb = build_grid()
    write_model(function_tmpdir, g, nb, True, False, False, **kwargs)
    r = run_mf6(targets["mf6"], function_tmpdir)
    assert r.returncode != 0
    assert message in " ".join(r.stdout.split())


def test_ihc3_requires_full3d(function_tmpdir, targets):
    g, nb = build_grid()
    write_model(function_tmpdir, g, nb, False, False, False)
    text = (function_tmpdir / "m.disu").read_text()
    # change the first connection of cell 1 to IHC = 3
    lines = text.split("\n")
    i = lines.index("  IHC") + 2
    first = lines[i].split()
    first[1] = "3"
    lines[i] = "      " + " ".join(first)
    (function_tmpdir / "m.disu").write_text("\n".join(lines))
    r = run_mf6(targets["mf6"], function_tmpdir)
    assert r.returncode != 0
    assert "IHC = 3 cannot be specified otherwise" in " ".join(r.stdout.split())


def rotation(a1, a2, a3):
    """Rotation matrix of the NPF K tensor axes (angles in radians)."""
    s1, c1, s2, c2, s3, c3 = (
        np.sin(a1),
        np.cos(a1),
        np.sin(a2),
        np.cos(a2),
        np.sin(a3),
        np.cos(a3),
    )
    return np.array(
        [
            [c1 * c2, c1 * s2 * s3 - s1 * c3, -c1 * s2 * c3 - s1 * s3],
            [s1 * c2, s1 * s2 * s3 + c1 * c3, -s1 * s2 * c3 + c1 * s3],
            [s2, -c2 * s3, c2 * c3],
        ]
    )


@pytest.mark.parametrize(
    "angles", [None, (25.0, 15.0, 10.0)], ids=["aligned", "rotated"]
)
def test_full3d_xt3d_and_specific_discharge(function_tmpdir, targets, angles):
    """
    Heads are fixed to the linear field h = a . x on all boundary cells.  For a
    homogeneous medium the XT3D solution is exactly h = a . x in the interior
    (also for a rotated, anisotropic K tensor) and the specific discharge is
    -K a in every interior cell.
    """
    g, nb = build_grid()
    a = np.array([3.0, -2.0, 5.0])
    hlin = g["xc"] * a[0] + g["yc"] * a[1] + g["zc"] * a[2]
    boundary = [
        nid(l, r, c)
        for l in range(nlay)
        for r in range(nrow)
        for c in range(ncol)
        if l in (0, nlay - 1) or r in (0, nrow - 1) or c in (0, ncol - 1)
    ]
    interior = np.setdiff1d(np.arange(N), boundary)
    k11, k22, k33 = 4.0, 1.5, 0.7
    npf_grid = _array("K22", np.full(N, k22)) + _array("K33", np.full(N, k33))
    if angles is None:
        kt = np.diag([k11, k22, k33])
    else:
        npf_grid += _array("ANGLE1", np.full(N, angles[0]))
        npf_grid += _array("ANGLE2", np.full(N, angles[1]))
        npf_grid += _array("ANGLE3", np.full(N, angles[2]))
        r_ = rotation(*np.radians(angles))
        kt = r_ @ np.diag([k11, k22, k33]) @ r_.T
    write_model(
        function_tmpdir,
        g,
        nb,
        True,
        False,
        False,
        npf_extra="  XT3D\n  SAVE_SPECIFIC_DISCHARGE\n  SAVE_FLOWS\n",
        accel="BICGSTAB",
        chd=[(n, hlin[n]) for n in boundary],
        wel=[(int(interior[0]), 0.0)],
        npf_grid=npf_grid.replace("K22", "K22"),
    )
    # K was written by write_model (3.0): replace it by k11
    npf = (function_tmpdir / "m.npf").read_text()
    npf = npf.replace(_array("K", np.full(N, 3.0)), _array("K", np.full(N, k11)))
    (function_tmpdir / "m.npf").write_text(npf)
    # save budget
    (function_tmpdir / "m.oc").write_text(
        "BEGIN OPTIONS\n  HEAD FILEOUT m.hds\n  BUDGET FILEOUT m.cbc\nEND OPTIONS\n"
        "BEGIN PERIOD 1\n  SAVE HEAD ALL\n  SAVE BUDGET ALL\nEND PERIOD\n"
    )
    r = run_mf6(targets["mf6"], function_tmpdir)
    assert r.returncode == 0, r.stdout[-2000:]
    h = heads(function_tmpdir)[-1]
    np.testing.assert_allclose(h[interior], hlin[interior], rtol=0, atol=1e-7)
    cbc = flopy.utils.CellBudgetFile(function_tmpdir / "m.cbc", precision="double")
    rec = cbc.get_data(text="DATA-SPDIS")[-1]
    q = np.zeros((N, 3))
    for e in rec:
        q[e["node"] - 1] = [e["qx"], e["qy"], e["qz"]]
    np.testing.assert_allclose(
        q[interior], np.tile(-(kt @ a), (interior.size, 1)), rtol=1e-6, atol=1e-6
    )
