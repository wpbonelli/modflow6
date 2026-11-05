"""Run mf6 binaries against one fixed input set; report the test's metrics.

Args: input_dir out_dir name=path/to/mf6 [name=path ...]
Metrics match test_mmr_to_mf6_dfw_regression. Appends to out_dir/results.jsonl.
"""

import hashlib
import json
import pathlib as pl
import re
import shutil
import subprocess
import sys

import flopy
import numpy as np

ANS = {"stage_all": 1.0372047024253908, "flow_all": 44.70210250910929}
input_dir, out_dir = map(pl.Path, sys.argv[1:3])
out_dir.mkdir(parents=True, exist_ok=True)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()[:12]


def metrics(ws):
    sobj = flopy.utils.HeadFile(ws / "drb_dfw.stage", text="STAGE")
    grb = flopy.mf6.utils.MfGrdFile(ws / "drb_dfw.disv1d.grb", verbose=False)
    stage = sobj.get_alldata().squeeze()
    bottom = grb.bot if hasattr(grb, "bot") else grb._datadict["BOTTOM"]
    stage = np.maximum(stage - bottom, 0)
    bud = flopy.utils.binaryfile.CellBudgetFile(ws / "drb_dfw.bud")
    flowja = bud.get_data(text="FLOW-JA-FACE")
    ia, ja = grb.ia, grb.ja
    flow = np.full_like(stage, np.nan, dtype=float)
    for tt in range(flow.shape[0]):
        out = np.zeros(ia.shape[0] - 1)
        q = flowja[tt].flatten()
        for n in range(grb.nodes):
            seg = q[ia[n] + 1 : ia[n + 1]]
            out[n] = -seg[seg < 0].sum()
        flow[tt, :] = out
    return stage.mean(), flow.mean()


for spec in sys.argv[3:]:
    name, exe = spec.split("=", 1)
    ws = out_dir / name
    if ws.exists():
        shutil.rmtree(ws)
    shutil.copytree(input_dir, ws)
    r = subprocess.run(
        [str(pl.Path(exe).resolve())], cwd=ws, capture_output=True, text=True
    )
    ok = "Normal termination" in r.stdout
    rec = {"name": name, "ok": ok}
    if ok:
        s, f = metrics(ws)
        lst = (ws / "mfsim.lst").read_text(errors="replace")
        m = re.findall(r"(\d+)\s+TOTAL ITERATIONS|TOTAL ITER\D+(\d+)", lst)
        rec.update(
            stage=float(s),
            flow=float(f),
            d_stage=float(s - ANS["stage_all"]),
            d_flow=float(f - ANS["flow_all"]),
            sha_stage=sha(ws / "drb_dfw.stage"),
            sha_bud=sha(ws / "drb_dfw.bud"),
            n_total_iter=len(re.findall(r"^\s*Solving:", r.stdout, re.M)),
        )
    else:
        rec["tail"] = r.stdout[-400:]
    print(json.dumps(rec), flush=True)
    with open(out_dir / "results.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    # keep only listing, drop bulky outputs but keep stage for direct compare
    for p in ws.iterdir():
        if p.name not in ("drb_dfw.stage", "mfsim.lst", "drb_dfw.lst"):
            if p.is_file():
                p.unlink()
