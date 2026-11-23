"""
Test the VTU output of the DISU Package with the FULL_3D option.

The options VTU_GEOMETRY (polyhedra of the cells, free-format file) and
VTU_FILEOUT (root of the output files) make the model write, whenever the
dependent variable is saved by the Output Control Package, one VTK XML
unstructured grid file per time step (HEAD and SPECIFIC_DISCHARGE for
GWF, CONCENTRATION or TEMPERATURE for GWT and GWE, plus CELL_VOLUME) and a
ParaView collection (.pvd).  The box grid of test_gwf_disu_full3d is used,
and the cells are written as hexahedra.
"""

import xml.etree.ElementTree as ET

import numpy as np
import pytest
from test_gwf_disu_full3d import (
    N,
    build_grid,
    heads,
    run_mf6,
    write_model,
)


def write_v3d(path, g):
    """Cells of the box as hexahedra: per cell nf=6 faces with 4 vertices."""
    dz = g["vol"] / (g["dx"] * g["dy"])
    pts, cells = [], []
    for n in range(N):
        x0, x1 = g["xc"][n] - g["dx"][n] / 2, g["xc"][n] + g["dx"][n] / 2
        y0, y1 = g["yc"][n] - g["dy"][n] / 2, g["yc"][n] + g["dy"][n] / 2
        z0, z1 = g["zc"][n] - dz[n] / 2, g["zc"][n] + dz[n] / 2
        c = [(x, y, z) for z in (z0, z1) for y in (y0, y1) for x in (x0, x1)]
        b = len(pts)
        pts += c
        # vertex order: index = ix + 2*iy + 4*iz (1-based below)
        f = [
            (0, 2, 3, 1),
            (4, 5, 7, 6),
            (0, 1, 5, 4),
            (2, 6, 7, 3),
            (0, 4, 6, 2),
            (1, 3, 7, 5),
        ]
        cells.append([[b + v + 1 for v in face] for face in f])
    text = f"NCELLS {N}\nNPOINTS {len(pts)}\n"
    text += "".join(f"{x:.15g} {y:.15g} {z:.15g}\n" for x, y, z in pts)
    for fc in cells:
        text += f"{len(fc)}\n" + "".join(f"4 {a} {b} {c} {d}\n" for a, b, c, d in fc)
    path.write_text(text)


def read_vtu(path):
    root = ET.parse(path).getroot()
    piece = root.find("UnstructuredGrid/Piece")
    out = {"ncells": int(piece.get("NumberOfCells"))}
    for da in piece.iter("DataArray"):
        name = da.get("Name")
        if name:
            out[name] = np.array(da.text.split(), dtype=float)
    return out


def add_vtu_options(ws, geom="m.v3d", fileout="out"):
    p = ws / "m.disu"
    opts = ""
    if geom:
        opts += f"  VTU_GEOMETRY {geom}\n"
    if fileout:
        opts += f"  VTU_FILEOUT {fileout}\n"
    p.write_text(p.read_text().replace("  FULL_3D\n", "  FULL_3D\n" + opts, 1))


def test_full3d_vtu_flow(function_tmpdir, targets):
    g, nb = build_grid()
    ws = function_tmpdir
    write_model(ws, g, nb, True, False, True, npf_extra="  SAVE_SPECIFIC_DISCHARGE\n")
    write_v3d(ws / "m.v3d", g)
    add_vtu_options(ws)
    r = run_mf6(targets["mf6"], ws)
    assert r.returncode == 0, r.stdout[-2000:]
    h = heads(ws)
    nstp = h.shape[0]
    files = sorted(f for f in ws.glob("out_*.vtu") if "vectors" not in f.name)
    assert len(files) == nstp == 6
    pvd = (ws / "out.pvd").read_text()
    assert all(f.name in pvd for f in files)
    for i, f in enumerate(files):
        d = read_vtu(f)
        assert d["ncells"] == N
        np.testing.assert_allclose(d["HEAD"], h[i], rtol=0, atol=1e-9)
        np.testing.assert_allclose(d["CELL_VOLUME"], g["vol"], rtol=1e-12)
        assert d["SPECIFIC_DISCHARGE"].size == 3 * N
    # the discharge is along the x axis (constant heads on x faces): mean qx > 0
    qx = read_vtu(files[-1])["SPECIFIC_DISCHARGE"].reshape(-1, 3)[:, 0]
    assert qx.mean() > 0.0
    # vectors at the cell centers: point cloud with the same discharge
    vfiles = sorted(ws.glob("out_vectors_*.vtu"))
    assert len(vfiles) == nstp
    assert all(f.name in (ws / "out_vectors.pvd").read_text() for f in vfiles)
    root = ET.parse(vfiles[-1]).getroot()
    assert int(root.find(".//Piece").get("NumberOfPoints")) == N
    names = [a.get("Name") for a in root.iter("DataArray")]
    assert "SPECIFIC_DISCHARGE" in names and "HEAD" in names
    # the binary grid file is written also for a fully 3D grid
    assert len(list(ws.glob("*.grb"))) == 1


@pytest.mark.parametrize(
    "geom, fileout, message",
    [
        ("m.v3d", None, "VTU_GEOMETRY and VTU_FILEOUT"),
        (None, "out", "VTU_GEOMETRY and VTU_FILEOUT"),
        ("missing.v3d", "out", "Could not open the VTU geometry file"),
    ],
)
def test_full3d_vtu_errors(function_tmpdir, targets, geom, fileout, message):
    g, nb = build_grid()
    ws = function_tmpdir
    write_model(ws, g, nb, True, False, False)
    write_v3d(ws / "m.v3d", g)
    add_vtu_options(ws, geom, fileout)
    r = run_mf6(targets["mf6"], ws)
    assert r.returncode != 0
    assert message in r.stdout


def test_vtu_requires_full3d(function_tmpdir, targets):
    g, nb = build_grid()
    ws = function_tmpdir
    write_model(ws, g, nb, False, False, False)
    write_v3d(ws / "m.v3d", g)
    p = ws / "m.disu"
    p.write_text(
        p.read_text().replace(
            "BEGIN OPTIONS\n",
            "BEGIN OPTIONS\n  VTU_GEOMETRY m.v3d\n  VTU_FILEOUT out\n",
            1,
        )
    )
    r = run_mf6(targets["mf6"], ws)
    assert r.returncode != 0
    assert "FULL_3D" in r.stdout
