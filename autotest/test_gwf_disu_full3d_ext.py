"""
Test the FULL_3D option of the DISU Package together with transport (GWT and
GWE), the GWF-GWF exchange and the RCH, EVT, SFR, MAW and LAK Packages.

Every test simulates the structured box of test_gwf_disu_full3d.py twice, as a
classic multi-layer DISU grid and as a FULL_3D DISU grid, and the results must
be identical.  The GWF-GWF exchange test splits the FULL_3D box in two models
connected through IHC = 3 exchanges and compares with the single model.
"""

import shutil
from pathlib import Path

import flopy
import numpy as np
import pytest
from test_gwf_disu_full3d import (
    N,
    _array,
    build_grid,
    ncol,
    nid,
    nlay,
    nrow,
    run_mf6,
    write_model,
)

SPDIS = "  SAVE_FLOWS\n  SAVE_SPECIFIC_DISCHARGE\n  SAVE_SATURATION\n"
TOP_CELLS = [nid(0, r, c) for r in range(nrow) for c in range(ncol)]


def _edit(path, old, new):
    p = Path(path)
    text = p.read_text()
    assert old in text
    p.write_text(text.replace(old, new, 1))


def add_vertices(ws, g):
    """Add VERTICES and CELL2D (required for specific discharge on classic DISU)."""
    verts, cells = [], []
    for n in range(N):
        x0, x1 = g["xc"][n] - 0.5 * g["dx"][n], g["xc"][n] + 0.5 * g["dx"][n]
        y0, y1 = g["yc"][n] - 0.5 * g["dy"][n], g["yc"][n] + 0.5 * g["dy"][n]
        base = len(verts)
        verts += [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        cells.append((n + 1, g["xc"][n], g["yc"][n], [base + i for i in (1, 2, 3, 4)]))
    text = (
        (Path(ws) / "m.disu")
        .read_text()
        .replace("NJA", f"NVERT {len(verts)}\n  NJA", 1)
    )
    text += "\nBEGIN VERTICES\n"
    text += "".join(f"  {i + 1} {x:.15g} {y:.15g}\n" for i, (x, y) in enumerate(verts))
    text += "END VERTICES\n\nBEGIN CELL2D\n"
    text += "".join(
        f"  {n} {x:.15g} {y:.15g} 4 {' '.join(map(str, v))}\n" for n, x, y, v in cells
    )
    (Path(ws) / "m.disu").write_text(text + "END CELL2D\n")


def add_package(ws, ftype, name, text):
    _edit(Path(ws) / "m.nam", "END PACKAGES", f"  {ftype}6 m.{name}\nEND PACKAGES")
    (Path(ws) / f"m.{name}").write_text(text)


def save_budget(ws):
    (Path(ws) / "m.oc").write_text(
        "BEGIN OPTIONS\n  HEAD FILEOUT m.hds\n  BUDGET FILEOUT m.cbc\nEND OPTIONS\n"
        "BEGIN PERIOD 1\n  SAVE HEAD ALL\n  SAVE BUDGET ALL\nEND PERIOD\n"
    )
    (Path(ws) / "m.nam").write_text(
        "BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n" + (Path(ws) / "m.nam").read_text()
    )


def hds(ws, name="m.hds"):
    hf = flopy.utils.HeadFile(Path(ws) / name, precision="double")
    return hf.get_data(totim=hf.get_times()[-1]).ravel()


def run_pair(function_tmpdir, targets, tag, build, **kwargs):
    """Run classic and FULL_3D versions; build(ws, full3d, g) adds packages."""
    g, nb = build_grid()
    out = {}
    for full3d in (False, True):
        ws = function_tmpdir / f"{tag}_{'full3d' if full3d else 'classic'}"
        write_model(ws, g, nb, full3d, False, False, **kwargs)
        build(ws, full3d, g)
        r = run_mf6(targets["mf6"], ws)
        assert r.returncode == 0, r.stdout[-2500:]
        out[full3d] = hds(ws)
    return out, g, nb


def add_exposed_area(ws, g):
    area = np.where(np.arange(N) < nrow * ncol, g["area"], 0.0)
    _edit(Path(ws) / "m.disu", "END GRIDDATA", _array("AREA", area) + "END GRIDDATA")


@pytest.mark.parametrize("adv", ["UPSTREAM", "CENTRAL", "TVD"])
@pytest.mark.parametrize("model", ["GWT", "GWE"])
def test_full3d_transport_matches_classic(function_tmpdir, targets, model, adv):
    g, nb = build_grid()
    out = {}
    for full3d in (False, True):
        ws = function_tmpdir / ("full3d" if full3d else "classic")
        write_model(ws, g, nb, full3d, False, False, npf_extra=SPDIS)
        save_budget(ws)
        if not full3d:
            add_vertices(ws, g)
        assert run_mf6(targets["mf6"], ws).returncode == 0
        tws = ws / "t"
        tws.mkdir()
        shutil.copy(ws / "m.disu", tws / "m.disu")
        (tws / "mfsim.nam").write_text(
            "BEGIN TIMING\n  TDIS6 m.tdis\nEND TIMING\n"
            f"BEGIN MODELS\n  {model}6 m.nam t\nEND MODELS\n"
            "BEGIN EXCHANGES\nEND EXCHANGES\n"
            "BEGIN SOLUTIONGROUP 1\n  IMS6 m.ims t\nEND SOLUTIONGROUP\n"
        )
        (tws / "m.tdis").write_text(
            "BEGIN DIMENSIONS\n  NPER 1\nEND DIMENSIONS\n"
            "BEGIN PERIODDATA\n  0.4 8 1.0\nEND PERIODDATA\n"
        )
        (tws / "m.ims").write_text(
            "BEGIN NONLINEAR\n  OUTER_DVCLOSE 1e-8\n  OUTER_MAXIMUM 20\nEND NONLINEAR\n"
            "BEGIN LINEAR\n  INNER_MAXIMUM 500\n  INNER_DVCLOSE 1e-10\n"
            "  INNER_RCLOSE 1e-9\n  LINEAR_ACCELERATION BICGSTAB\nEND LINEAR\n"
        )
        gwt = model == "GWT"
        bc = "CNC" if gwt else "CTP"
        var = "CONCENTRATION" if gwt else "TEMPERATURE"
        pk = "  MST6 m.mst\n  DSP6 m.dsp\n" if gwt else "  EST6 m.mst\n  CND6 m.dsp\n"
        (tws / "m.nam").write_text(
            "BEGIN PACKAGES\n  DISU6 m.disu\n  IC6 m.ic\n  FMI6 m.fmi\n"
            f"  SSM6 m.ssm\n  OC6 m.oc\n  {bc}6 m.bc\n  ADV6 m.adv\n{pk}END PACKAGES\n"
        )
        (tws / "m.ic").write_text(
            "BEGIN GRIDDATA\n" + _array("STRT", np.zeros(N)) + "END GRIDDATA\n"
        )
        (tws / "m.fmi").write_text(
            "BEGIN PACKAGEDATA\n  GWFHEAD FILEIN ../m.hds\n"
            "  GWFBUDGET FILEIN ../m.cbc\nEND PACKAGEDATA\n"
        )
        (tws / "m.ssm").write_text("BEGIN SOURCES\nEND SOURCES\n")
        por = _array("POROSITY", np.full(N, 0.3))
        if gwt:
            (tws / "m.mst").write_text(f"BEGIN GRIDDATA\n{por}END GRIDDATA\n")
            disp = (
                _array("DIFFC", np.full(N, 1e-3))
                + _array("ALH", np.full(N, 0.5))
                + _array("ATH1", np.full(N, 0.1))
                + _array("ATV", np.full(N, 0.05))
            )
        else:
            (tws / "m.mst").write_text(
                f"BEGIN GRIDDATA\n{por}"
                + _array("HEAT_CAPACITY_SOLID", np.full(N, 800.0))
                + _array("DENSITY_SOLID", np.full(N, 2650.0))
                + "END GRIDDATA\n"
            )
            disp = (
                _array("KTS", np.full(N, 2.0))
                + _array("KTW", np.full(N, 0.6))
                + _array("ALH", np.full(N, 0.5))
                + _array("ATH1", np.full(N, 0.1))
                + _array("ATV", np.full(N, 0.05))
            )
        (tws / "m.dsp").write_text(f"BEGIN GRIDDATA\n{disp}END GRIDDATA\n")
        (tws / "m.adv").write_text(f"BEGIN OPTIONS\n  SCHEME {adv}\nEND OPTIONS\n")
        (tws / "m.oc").write_text(
            f"BEGIN OPTIONS\n  {var} FILEOUT m.ucn\nEND OPTIONS\n"
            f"BEGIN PERIOD 1\n  SAVE {var} ALL\nEND PERIOD\n"
        )
        left = [nid(l, r, 0) for l in range(nlay) for r in range(nrow)]
        val = 1.0 if gwt else 30.0
        (tws / "m.bc").write_text(
            f"BEGIN DIMENSIONS\n  MAXBOUND {len(left)}\nEND DIMENSIONS\n"
            "BEGIN PERIOD 1\n"
            + "".join(f"  {c + 1} {val}\n" for c in left)
            + "END PERIOD\n"
        )
        r = run_mf6(targets["mf6"], tws)
        assert r.returncode == 0, r.stdout[-2500:]
        cf = flopy.utils.HeadFile(tws / "m.ucn", precision="double", text=var.lower())
        out[full3d] = cf.get_data(totim=cf.get_times()[-1]).ravel()
    np.testing.assert_allclose(out[True], out[False], rtol=0, atol=1e-6)
    # the transport solution is not trivial
    assert out[True].max() - out[True].min() > 0.2 * (1.0 if model == "GWT" else 30.0)


def test_full3d_rch_evt(function_tmpdir, targets):
    def build(ws, full3d, g):
        if full3d:
            add_exposed_area(ws, g)
        add_package(
            ws,
            "RCH",
            "rch",
            "BEGIN OPTIONS\n  FIXED_CELL\nEND OPTIONS\n"
            f"BEGIN DIMENSIONS\n  MAXBOUND {len(TOP_CELLS)}\nEND DIMENSIONS\n"
            "BEGIN PERIOD 1\n"
            + "".join(f"  {c + 1} {0.3 + 0.02 * i}\n" for i, c in enumerate(TOP_CELLS))
            + "END PERIOD\n",
        )
        add_package(
            ws,
            "EVT",
            "evt",
            "BEGIN OPTIONS\n  FIXED_CELL\nEND OPTIONS\n"
            f"BEGIN DIMENSIONS\n  MAXBOUND {len(TOP_CELLS)}\n  NSEG 1\nEND DIMENSIONS\n"
            "BEGIN PERIOD 1\n"
            + "".join(f"  {c + 1} 9.5 0.25 1.5\n" for c in TOP_CELLS)
            + "END PERIOD\n",
        )

    chd = [(nid(l, r, 0), 10.0) for l in range(nlay) for r in range(nrow)]
    out, _, _ = run_pair(function_tmpdir, targets, "rch", build, chd=chd)
    np.testing.assert_allclose(out[True], out[False], rtol=0, atol=1e-9)
    assert out[True].max() - out[True].min() > 0.3


def test_full3d_rch_requires_area(function_tmpdir, targets):
    g, nb = build_grid()
    write_model(function_tmpdir, g, nb, True, False, False)
    add_package(
        function_tmpdir,
        "RCH",
        "rch",
        "BEGIN OPTIONS\n  FIXED_CELL\nEND OPTIONS\n"
        "BEGIN DIMENSIONS\n  MAXBOUND 1\nEND DIMENSIONS\n"
        "BEGIN PERIOD 1\n  1 0.1\nEND PERIOD\n",
    )
    r = run_mf6(targets["mf6"], function_tmpdir)
    assert r.returncode != 0
    assert "requires the exposed surface area" in r.stdout


@pytest.mark.parametrize("eqn", ["SPECIFIED", "MEAN"])
def test_full3d_maw(function_tmpdir, targets, eqn):
    cells = [nid(l, 2, 2) for l in range(nlay)]

    def build(ws, full3d, g):
        ze = [g["top"][c] for c in cells] + [g["bot"][cells[-1]]]
        s = "BEGIN DIMENSIONS\n  NMAWWELLS 1\nEND DIMENSIONS\nBEGIN PACKAGEDATA\n"
        s += f"  1 0.15 {ze[-1] - 1.0} 5.0 {eqn} {nlay}\nEND PACKAGEDATA\n"
        s += "BEGIN CONNECTIONDATA\n"
        for l, c in enumerate(cells):
            hk, rs = (4.0 + l, 0.3) if eqn == "MEAN" else (2.0 + l, 0.0)
            s += f"  1 {l + 1} {c + 1} {ze[l]} {ze[l + 1]} {hk} {rs}\n"
        s += "END CONNECTIONDATA\nBEGIN PERIOD 1\n  1 RATE -3.0\nEND PERIOD\n"
        add_package(ws, "MAW", "maw", s)

    out, _, _ = run_pair(function_tmpdir, targets, "maw", build)
    np.testing.assert_allclose(out[True], out[False], rtol=0, atol=1e-9)


def test_full3d_maw_thiem_not_supported(function_tmpdir, targets):
    g, nb = build_grid()
    write_model(function_tmpdir, g, nb, True, False, False)
    s = "BEGIN DIMENSIONS\n  NMAWWELLS 1\nEND DIMENSIONS\n"
    s += "BEGIN PACKAGEDATA\n  1 0.15 -3.0 0.0 THIEM 1\nEND PACKAGEDATA\n"
    s += "BEGIN CONNECTIONDATA\n  1 1 1 0.0 -1.0 0.0 0.0\nEND CONNECTIONDATA\n"
    s += "BEGIN PERIOD 1\n  1 RATE -1.0\nEND PERIOD\n"
    add_package(function_tmpdir, "MAW", "maw", s)
    r = run_mf6(targets["mf6"], function_tmpdir)
    assert r.returncode != 0
    assert "must be SPECIFIED or MEAN" in r.stdout


def test_full3d_sfr(function_tmpdir, targets):
    reach = [nid(0, r, 2) for r in range(nrow)]

    def build(ws, full3d, g):
        n = len(reach)
        s = f"BEGIN DIMENSIONS\n  NREACHES {n}\nEND DIMENSIONS\nBEGIN PACKAGEDATA\n"
        for i, c in enumerate(reach):
            ncon = 1 if i in (0, n - 1) else 2
            s += f"  {i + 1} {c + 1} 1.0 2.0 0.01 -0.3 0.2 0.8 0.04 {ncon} 1.0 0\n"
        s += "END PACKAGEDATA\nBEGIN CONNECTIONDATA\n"
        for i in range(n):
            cs = ([i] if i > 0 else []) + ([-(i + 2)] if i < n - 1 else [])
            s += f"  {i + 1} " + " ".join(map(str, cs)) + "\n"
        s += "END CONNECTIONDATA\nBEGIN PERIOD 1\n  1 INFLOW 1.5\nEND PERIOD\n"
        add_package(ws, "SFR", "sfr", "BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n" + s)

    out, _, _ = run_pair(function_tmpdir, targets, "sfr", build)
    np.testing.assert_allclose(out[True], out[False], rtol=0, atol=1e-9)


@pytest.mark.parametrize("ictype", ["HORIZONTAL", "VERTICAL"])
def test_full3d_lak(function_tmpdir, targets, ictype):
    cells = [nid(l, 2, 2) for l in range(nlay)]

    def lak(ws, g):
        ze = [g["top"][c] for c in cells] + [g["bot"][cells[-1]]]
        s = "BEGIN OPTIONS\n  SAVE_FLOWS\nEND OPTIONS\n"
        s += "BEGIN DIMENSIONS\n  NLAKES 1\n  NOUTLETS 0\n  NTABLES 1\nEND DIMENSIONS\n"
        s += f"BEGIN PACKAGEDATA\n  1 5.0 {nlay}\nEND PACKAGEDATA\n"
        s += "BEGIN CONNECTIONDATA\n"
        for l, c in enumerate(cells):
            bot, top = ze[l + 1] + 0.05, ze[l] - 0.05
            s += f"  1 {l + 1} {c + 1} {ictype} 0.5 {bot} {top} 0.6 0.8\n"
        s += "END CONNECTIONDATA\nBEGIN TABLES\n  1 TAB6 FILEIN m.laktab\nEND TABLES\n"
        s += "BEGIN PERIOD 1\n  1 STATUS ACTIVE\n  1 RAINFALL 0.0\nEND PERIOD\n"
        add_package(ws, "LAK", "lak", s)
        (Path(ws) / "m.laktab").write_text(
            "BEGIN DIMENSIONS\n  NROW 3\n  NCOL 3\nEND DIMENSIONS\nBEGIN TABLE\n"
            f"  {ze[-1] - 0.5} 0.0 0.0\n  5.0 50.0 10.0\n  12.0 150.0 10.0\nEND TABLE\n"
        )

    if ictype == "VERTICAL":
        g, nb = build_grid()
        write_model(function_tmpdir, g, nb, True, False, False)
        lak(function_tmpdir, g)
        r = run_mf6(targets["mf6"], function_tmpdir)
        assert r.returncode != 0
        assert "must be HORIZONTAL" in r.stdout
        return
    out, _, _ = run_pair(function_tmpdir, targets, "lak", lambda ws, f, g: lak(ws, g))
    np.testing.assert_allclose(out[True], out[False], rtol=0, atol=1e-9)


def test_full3d_uzf_not_supported(function_tmpdir, targets):
    g, nb = build_grid()
    write_model(function_tmpdir, g, nb, True, False, False)
    add_package(
        function_tmpdir,
        "UZF",
        "uzf",
        "BEGIN DIMENSIONS\n  NUZFCELLS 1\n  NTRAILWAVES 7\n  NWAVESETS 40\n"
        "END DIMENSIONS\nBEGIN PACKAGEDATA\n"
        "  1 1 1 1 0.1 0.01 0.0 0.3 0.1 3.5\nEND PACKAGEDATA\n"
        "BEGIN PERIOD 1\n  1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\nEND PERIOD\n",
    )
    r = run_mf6(targets["mf6"], function_tmpdir)
    assert r.returncode != 0
    assert "UZF Package is not supported" in r.stdout


def test_full3d_gwfgwf_exchange(function_tmpdir, targets):
    g, nb = build_grid()
    cut = 3
    in_a = np.array([(n % ncol) < cut for n in range(N)])
    idx = {"a": np.where(in_a)[0], "b": np.where(~in_a)[0]}
    loc = {}
    for tag, ids in idx.items():
        for i, n in enumerate(ids):
            loc[n] = (tag, i)
    kvals = np.array([3.0 + 0.2 * (n % 7) for n in range(N)])
    chd_all = [(nid(l, r, 0), 10.0) for l in range(nlay) for r in range(nrow)]
    chd_all += [(nid(l, r, ncol - 1), 2.0) for l in range(nlay) for r in range(nrow)]
    wel_all = [(nid(1, 2, 2), -0.7), (nid(2, 1, 3), 0.4)]

    # reference: single FULL_3D model with the same K field
    single = function_tmpdir / "single"
    write_model(single, g, nb, True, False, False, chd=chd_all, wel=wel_all)
    _edit(single / "m.npf", _array("K", np.full(N, 3.0)), _array("K", kvals))
    assert run_mf6(targets["mf6"], single).returncode == 0
    href = hds(single)

    split = function_tmpdir / "split"
    split.mkdir()
    for tag in ("a", "b"):
        ids = idx[tag]
        sel = {n: i for i, n in enumerate(ids)}
        sub = split / tag
        sub_nb = []
        for n in ids:
            sub_nb.append([(sel[t[0]],) + t[1:] for t in nb[n] if t[0] in sel])
        sg = {k: v[ids] for k, v in g.items()}
        chd = [(loc[c][1], h) for c, h in chd_all if loc[c][0] == tag]
        wel = [(loc[c][1], q) for c, q in wel_all if loc[c][0] == tag]
        sub.mkdir()
        # write_disu / write_model use the global node count, so write the
        # sub-model files with a local node count
        import test_gwf_disu_full3d as base

        old = base.N
        base.N = len(ids)
        try:
            write_model(sub, sg, sub_nb, True, False, False, chd=chd, wel=wel)
        finally:
            base.N = old
        _edit(
            sub / "m.npf", _array("K", np.full(len(ids), 3.0)), _array("K", kvals[ids])
        )
        for fn in ("m.nam", "m.oc"):
            (sub / fn).write_text((sub / fn).read_text().replace(" m.", f" {tag}/m."))
    shutil.copy(split / "a" / "m.tdis", split / "m.tdis")
    shutil.copy(split / "a" / "m.ims", split / "m.ims")
    exg = []
    for n in np.where(in_a)[0]:
        for m, ih, cl, hw, an, fa in nb[n]:
            if not in_a[m]:
                cl2 = next(t[2] for t in nb[m] if t[0] == n)
                exg.append((loc[n][1] + 1, loc[m][1] + 1, 3, cl, cl2, fa))
    (split / "mfsim.nam").write_text(
        "BEGIN TIMING\n  TDIS6 m.tdis\nEND TIMING\n"
        "BEGIN MODELS\n  GWF6 a/m.nam a\n  GWF6 b/m.nam b\nEND MODELS\n"
        "BEGIN EXCHANGES\n  GWF6-GWF6 m.exg a b\nEND EXCHANGES\n"
        "BEGIN SOLUTIONGROUP 1\n  IMS6 m.ims a b\nEND SOLUTIONGROUP\n"
    )
    (split / "m.exg").write_text(
        f"BEGIN DIMENSIONS\n  NEXG {len(exg)}\nEND DIMENSIONS\nBEGIN EXCHANGEDATA\n"
        + "".join(f"  {a} {b} {i} {c1} {c2} {h}\n" for a, b, i, c1, c2, h in exg)
        + "END EXCHANGEDATA\n"
    )
    r = run_mf6(targets["mf6"], split)
    assert r.returncode == 0, r.stdout[-2500:]
    hs = np.empty(N)
    hs[idx["a"]] = hds(split / "a")
    hs[idx["b"]] = hds(split / "b")
    np.testing.assert_allclose(hs, href, rtol=0, atol=1e-9)
    assert href.max() - href.min() > 5.0
