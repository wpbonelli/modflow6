"""Run each nightly binary against each input set and print a summary.

Args: asset(linux|win64|macarm) input_dir [input_dir ...]
"""

import glob
import json
import os
import pathlib as pl
import subprocess
import sys
import zipfile
from concurrent.futures import ThreadPoolExecutor

TAGS = """20260923 20261008 20261009 20261010 20261011 20261012 20261013
20261014 20261021 20261022 20261023 20261024 20261025 20261026 20261027
20261028 20261029 20261031 20261103 20261104""".split()
REPO = "MODFLOW-ORG/modflow6-nightly-build"
asset, inputs = sys.argv[1], sys.argv[2:]
here = pl.Path(__file__).parent
bins = pl.Path("bins")


def fetch(tag):
    d = bins / tag
    d.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["gh", "release", "download", tag, "-R", REPO, "-p", f"{asset}.zip"],
        cwd=d,
        check=True,
    )
    with zipfile.ZipFile(d / f"{asset}.zip") as z:
        z.extractall(d)
    exe = [p for p in glob.glob(f"{d}/**/bin/mf6*", recursive=True)
           if pl.Path(p).name in ("mf6", "mf6.exe")][0]
    os.chmod(exe, 0o755)
    ver = subprocess.run([exe, "--version"], capture_output=True, text=True)
    print("BIN", tag, ver.stdout.strip(), flush=True)
    return tag, exe


def run(job):
    inp, tag, exe = job
    out = pl.Path("out") / pl.Path(inp).name / tag
    r = subprocess.run(
        [sys.executable, str(here / "run_matrix.py"), inp, str(out),
         f"{tag}={exe}"],
        capture_output=True,
        text=True,
    )
    if r.returncode:
        print("ERR", inp, tag, r.stderr[-800:], flush=True)
        return inp, {"name": tag, "ok": False, "tail": r.stderr[-300:]}
    rec = json.loads((out / "results.jsonl").read_text().splitlines()[-1])
    return inp, rec


with ThreadPoolExecutor(4) as ex:
    exes = list(ex.map(fetch, TAGS))
    jobs = [(i, t, e) for i in inputs for t, e in exes]
    results = list(ex.map(run, jobs))

for inp in inputs:
    print(f"\n=== input {inp}")
    for ii, rec in results:
        if ii != inp:
            continue
        if rec["ok"]:
            print(
                f"RES {pl.Path(inp).name} {rec['name']} "
                f"dflow={rec['d_flow']:+.6e} dstage={rec['d_stage']:+.6e} "
                f"flow={rec['flow']:.12f} stage={rec['stage']:.12f} "
                f"{rec['sha_stage']} {rec['sha_bud']}"
            )
        else:
            print(f"RES {pl.Path(inp).name} {rec['name']} FAILED "
                  f"{rec.get('tail', '')[-200:]!r}")

pl.Path("summary.json").write_text(json.dumps(results, indent=1))
