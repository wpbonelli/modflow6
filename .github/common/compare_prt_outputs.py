"""
Helpers for the PRT compiler-comparison workflow (.github/workflows/prt-compilers.yml).

Subcommands
-----------
config     Resolve a (os, toolchain, fp) matrix entry to a variant name, a
           setup-fortran compiler/version, FFLAGS and extra meson arguments,
           and append them to $GITHUB_ENV.
buildinfo  Record the compile command meson used for a PRT source file.
rerun      Re-run every PRT-only simulation under a directory with a given
           mf6 binary, after deleting all particle track output. This is
           used to rerun PRT with one build on flow model output produced
           by another build, to separate GWF differences from PRT ones.
compare    Compare track (and head) output across variants and write a
           markdown report plus a JSON summary.

Track output is compared from the binary track files where present (exact),
otherwise from the CSV track files. CSV values are printed with the G0 edit
descriptor, whose digit count is compiler-dependent, so CSV comparisons are
only exact to printed precision.
"""

import argparse
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import OrderedDict
from pathlib import Path

# ---------------------------------------------------------------------------
# variant configuration
# ---------------------------------------------------------------------------

REFERENCE = "ifort-default"
GFORTRAN_REFERENCE = "gfortran-default"

# fp name -> (meson args, linux intel FFLAGS, windows intel FFLAGS)
INTEL_FP = {
    "default": ("", "", ""),
    "o0-strict": ("-Doptimization=0", "-fp-model strict", "/fp:strict"),
    "precise": ("", "-fp-model precise", "/fp:precise"),
    "precise-nofma": (
        "",
        "-fp-model precise -no-fma -fimf-arch-consistency=true",
        "/fp:precise /Qfma- /Qimf-arch-consistency:true",
    ),
}

# fp name -> (meson args, FFLAGS)
GNU_FP = {
    "default": ("", ""),
    "o0": ("-Doptimization=0", ""),
    # flush out reads of uninitialized variables: signaling NaN reals
    # trap on first use, integers get an implausible value
    "uninit": (
        "-Doptimization=0",
        "-finit-real=snan -finit-integer=-2147483647 -finit-logical=true "
        "-ffpe-trap=invalid,zero,overflow -fcheck=bounds",
    ),
}


def os_short(os_name):
    os_name = os_name.lower()
    if os_name.startswith("ubuntu") or os_name.startswith("linux"):
        return "linux"
    if os_name.startswith("windows"):
        return "windows"
    if os_name.startswith("macos"):
        return "macos"
    raise ValueError(f"unknown os: {os_name}")


def resolve_variant(os_name, toolchain, fp, versions):
    oss = os_short(os_name)
    if toolchain in ("ifort", "ifx"):
        meson_args, lin, win = INTEL_FP[fp]
        fflags = win if oss == "windows" else lin
        if oss == "windows" and toolchain == "ifx":
            # ifx 2025.2 on Windows hits an internal compiler error (C0000005)
            # with /Qimf-arch-consistency:true, even on blas1_d.f90
            fflags = fflags.replace(" /Qimf-arch-consistency:true", "")
        compiler = "intel-classic" if toolchain == "ifort" else "intel"
    elif toolchain == "gfortran":
        meson_args, fflags = GNU_FP[fp]
        compiler = "gcc"
    else:
        raise ValueError(f"unknown toolchain: {toolchain}")
    return {
        "VARIANT": f"{oss}-{toolchain}-{fp}",
        "REF_VARIANT": f"{oss}-{REFERENCE}",
        "FC_COMPILER": compiler,
        "FC_VERSION": versions[toolchain],
        "FFLAGS": fflags,
        "MESON_ARGS": meson_args,
    }


def cmd_config(args):
    versions = {
        "ifort": args.ifort_version,
        "ifx": args.ifx_version,
        "gfortran": args.gcc_version,
    }
    cfg = resolve_variant(args.os, args.toolchain, args.fp, versions)
    lines = [f"{k}={v}" for k, v in cfg.items()]
    print("\n".join(lines))
    env_file = os.environ.get("GITHUB_ENV")
    if env_file:
        with open(env_file, "a") as f:
            f.write("\n".join(lines) + "\n")


def cmd_buildinfo(args):
    cc = Path(args.builddir) / "compile_commands.json"
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not cc.is_file():
        out.write_text("compile_commands.json not found\n")
        return
    entries = json.loads(cc.read_text())
    pick = None
    for e in entries:
        if "ParticleTracker" in e.get("file", "") or "ParticleTracker" in e.get(
            "command", ""
        ):
            pick = e
            break
    pick = pick or (entries[0] if entries else None)
    out.write_text((pick or {}).get("command", "no entries") + "\n")


# ---------------------------------------------------------------------------
# reading output
# ---------------------------------------------------------------------------

TRACK_FIELDS = (
    "kper kstp imdl iprp irpt ilay icell izone istatus ireason trelease t x y z name"
).split()
TRACK_STRUCT = struct.Struct("<10i5d40s")
DISCRETE = ("kper", "kstp", "ilay", "icell", "izone", "istatus", "ireason")
CONTINUOUS = ("t", "x", "y", "z")
SHOW = ("kper", "kstp", "ilay", "icell", "izone", "istatus", "ireason") + CONTINUOUS
IDX = {k: i for i, k in enumerate(TRACK_FIELDS)}


def is_track_csv(path):
    try:
        with open(path) as f:
            head = f.readline()
    except (OSError, UnicodeDecodeError):
        return False
    return head.replace(" ", "").lower().startswith("kper,kstp,imdl,iprp,irpt")


def find_track_files(root):
    """Map relative path -> kind ('bin' or 'csv'). Binary track files are
    recognized by their '.hdr' sibling. A CSV is only used if there is no
    binary track file in the same directory."""
    root = Path(root)
    found = {}
    bin_dirs = set()
    for hdr in root.rglob("*.hdr"):
        f = hdr.with_suffix("")
        if f.is_file():
            found[f.relative_to(root).as_posix()] = "bin"
            bin_dirs.add(f.parent)
    for csv in root.rglob("*.csv"):
        if csv.parent in bin_dirs:
            continue
        if is_track_csv(csv):
            found[csv.relative_to(root).as_posix()] = "csv"
    return found


def read_track(path, kind):
    """Return OrderedDict particle key -> list of records (tuples)."""
    rows = []
    if kind == "bin":
        data = Path(path).read_bytes()
        if len(data) % TRACK_STRUCT.size:
            raise ValueError(
                f"{path}: size {len(data)} not a multiple of {TRACK_STRUCT.size}"
            )
        for rec in TRACK_STRUCT.iter_unpack(data):
            rows.append(rec[:15])
    else:
        with open(path) as f:
            f.readline()
            for line in f:
                p = line.strip().split(",")
                if len(p) < 15:
                    continue
                rows.append(
                    tuple(int(v) for v in p[:10]) + tuple(float(v) for v in p[10:15])
                )
    particles = OrderedDict()
    for r in rows:
        key = (r[IDX["imdl"]], r[IDX["iprp"]], r[IDX["irpt"]], r[IDX["trelease"]])
        particles.setdefault(key, []).append(r)
    return particles


def read_heads(path):
    """Read an MF6 binary head file into a list of numpy arrays."""
    import numpy as np

    data = Path(path).read_bytes()
    hdr = struct.Struct("<2i2d16s3i")
    out = []
    pos = 0
    while pos + hdr.size <= len(data):
        kstp, kper, pertim, totim, text, n1, n2, ilay = hdr.unpack_from(data, pos)
        pos += hdr.size
        text = text.decode(errors="replace").strip().upper()
        n = n2 - n1 + 1 if text.endswith("HEADU") else n1 * n2
        if n <= 0 or pos + 8 * n > len(data):
            break
        out.append(np.frombuffer(data, dtype="<f8", count=n, offset=pos))
        pos += 8 * n
    return out


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------


def close(a, b, rtol, atol):
    return abs(a - b) <= atol + rtol * abs(b)


def compare_particle(ref, cand, rtol, atol):
    """Return (status, info). status: exact | close | diverged."""
    n = min(len(ref), len(cand))
    exact = True
    maxd = {k: 0.0 for k in CONTINUOUS}
    for i in range(n):
        r, c = ref[i], cand[i]
        for k in DISCRETE:
            if r[IDX[k]] != c[IDX[k]]:
                return "diverged", {"index": i, "reason": k, "maxd": maxd}
        for k in CONTINUOUS:
            a, b = c[IDX[k]], r[IDX[k]]
            if a != b:
                exact = False
                d = abs(a - b)
                maxd[k] = max(maxd[k], d)
                if not close(a, b, rtol, atol):
                    return "diverged", {"index": i, "reason": k, "maxd": maxd}
    if len(ref) != len(cand):
        return "diverged", {"index": n, "reason": "length", "maxd": maxd}
    return ("exact" if exact else "close"), {"maxd": maxd}


def fate(recs):
    last = recs[-1]
    return (last[IDX["istatus"]], last[IDX["ilay"]], last[IDX["icell"]])


def fmt_rec(r):
    if r is None:
        return "—"
    vals = []
    for k in SHOW:
        v = r[IDX[k]]
        vals.append(f"{v:.17g}" if isinstance(v, float) else str(v))
    return " ".join(vals)


def compare_track_file(ref_path, cand_path, kind, rtol, atol, max_examples):
    ref = read_track(ref_path, kind)
    cand = read_track(cand_path, kind)
    res = {
        "kind": kind,
        "particles": len(ref),
        "exact": 0,
        "close": 0,
        "diverged": 0,
        "missing": 0,
        "extra": 0,
        "fate_changed": 0,
        "reasons": {},
        "maxd": {k: 0.0 for k in CONTINUOUS},
        "examples": [],
    }
    for key, rrecs in ref.items():
        crecs = cand.get(key)
        if crecs is None:
            res["missing"] += 1
            res["fate_changed"] += 1
            continue
        status, info = compare_particle(rrecs, crecs, rtol, atol)
        res[status] += 1
        for k in CONTINUOUS:
            res["maxd"][k] = max(res["maxd"][k], info["maxd"][k])
        if fate(rrecs) != fate(crecs):
            res["fate_changed"] += 1
        if status == "diverged":
            reason = info["reason"]
            res["reasons"][reason] = res["reasons"].get(reason, 0) + 1
            if len(res["examples"]) < max_examples:
                i = info["index"]
                res["examples"].append(
                    {
                        "particle": f"imdl={key[0]} iprp={key[1]} irpt={key[2]} "
                        f"trelease={key[3]:.17g}",
                        "index": i,
                        "reason": reason,
                        "nrec": [len(rrecs), len(crecs)],
                        "prev": fmt_rec(rrecs[i - 1]) if i > 0 else None,
                        "ref": fmt_rec(rrecs[i] if i < len(rrecs) else None),
                        "cand": fmt_rec(crecs[i] if i < len(crecs) else None),
                    }
                )
    res["extra"] = sum(1 for k in cand if k not in ref)
    if res["diverged"] or res["missing"] or res["extra"]:
        res["status"] = "diverged"
    elif res["close"]:
        res["status"] = "close"
    else:
        res["status"] = "exact"
    return res


def compare_heads(ref_path, cand_path):
    import numpy as np

    r = read_heads(ref_path)
    c = read_heads(cand_path)
    if len(r) != len(c):
        return {"status": "diverged", "maxd": math.inf, "dry_mismatch": 0}
    maxd = 0.0
    dry = 0
    for a, b in zip(r, c):
        if a.shape != b.shape:
            return {"status": "diverged", "maxd": math.inf, "dry_mismatch": 0}
        da, db = np.abs(a) >= 1e29, np.abs(b) >= 1e29
        dry += int(np.count_nonzero(da != db))
        m = ~(da | db)
        if m.any():
            maxd = max(maxd, float(np.max(np.abs(a[m] - b[m]))))
    exact = maxd == 0.0 and dry == 0 and all(np.array_equal(a, b) for a, b in zip(r, c))
    return {
        "status": "exact" if exact else "differs",
        "maxd": maxd,
        "dry_mismatch": dry,
    }


def compare_dirs(ref_root, cand_root, rtol, atol, max_examples, heads=True):
    ref_tracks = find_track_files(ref_root)
    cand_tracks = find_track_files(cand_root)
    files = {}
    for rel, kind in sorted(ref_tracks.items()):
        if rel not in cand_tracks:
            files[rel] = {"status": "missing", "kind": kind}
            continue
        try:
            files[rel] = compare_track_file(
                Path(ref_root) / rel,
                Path(cand_root) / rel,
                kind,
                rtol,
                atol,
                max_examples,
            )
        except Exception as e:
            files[rel] = {"status": "error", "kind": kind, "error": str(e)}
    hds = {}
    if heads:
        for h in sorted(Path(ref_root).rglob("*.hds")):
            rel = h.relative_to(ref_root).as_posix()
            ch = Path(cand_root) / rel
            if not ch.is_file():
                hds[rel] = {"status": "missing"}
                continue
            try:
                hds[rel] = compare_heads(h, ch)
            except Exception as e:
                hds[rel] = {"status": "error", "error": str(e)}
    return {"tracks": files, "heads": hds}


def summarize(cmp):
    t = cmp["tracks"].values()
    s = {
        "files": len(cmp["tracks"]),
        "exact": sum(1 for f in t if f["status"] == "exact"),
        "close": sum(1 for f in t if f["status"] == "close"),
        "diverged": sum(1 for f in t if f["status"] == "diverged"),
        "missing": sum(1 for f in t if f["status"] in ("missing", "error")),
        "particles": sum(f.get("particles", 0) for f in t),
        "pdiverged": sum(
            f.get("diverged", 0) + f.get("missing", 0) + f.get("extra", 0) for f in t
        ),
        "fate": sum(f.get("fate_changed", 0) for f in t),
        "maxdxyz": max(
            [max(f["maxd"][k] for k in ("x", "y", "z")) for f in t if "maxd" in f]
            or [0.0]
        ),
        "maxdt": max([f["maxd"]["t"] for f in t if "maxd" in f] or [0.0]),
    }
    h = cmp["heads"].values()
    s["hfiles"] = len(cmp["heads"])
    s["hdiff"] = sum(1 for x in h if x["status"] != "exact")
    s["maxdh"] = max([x.get("maxd", 0.0) for x in h] or [0.0])
    return s


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------


def g(x):
    if x == 0:
        return "0"
    if x == math.inf:
        return "inf"
    return f"{x:.2e}"


def junit_counts(path):
    try:
        root = ET.parse(path).getroot()
    except Exception:
        return None
    suites = [root] if root.tag == "testsuite" else list(root.iter("testsuite"))
    c = {"tests": 0, "failures": 0, "errors": 0, "skipped": 0, "failed": []}
    for s in suites:
        for k in ("tests", "failures", "errors", "skipped"):
            c[k] += int(s.get(k, 0))
        for tc in s.iter("testcase"):
            if tc.find("failure") is not None or tc.find("error") is not None:
                c["failed"].append(f"{tc.get('classname')}::{tc.get('name')}")
    return c


def read_text(p, n=None):
    p = Path(p)
    if not p.is_file():
        return ""
    lines = p.read_text(errors="replace").splitlines()
    return "\n".join(lines[:n] if n else lines)


SUMMARY_HEADER = (
    "| comparison | track files | exact | close | diverged | missing "
    "| particles diverged | fates changed | max Δxyz (pre-divergence) | max Δt "
    "| head files differing | max Δh |\n"
    "|---|---|---|---|---|---|---|---|---|---|---|---|\n"
)


def summary_row(label, s):
    return (
        f"| {label} | {s['files']} | {s['exact']} | {s['close']} | {s['diverged']} "
        f"| {s['missing']} | {s['pdiverged']} / {s['particles']} | {s['fate']} "
        f"| {g(s['maxdxyz'])} | {g(s['maxdt'])} | {s['hdiff']} / {s['hfiles']} "
        f"| {g(s['maxdh'])} |\n"
    )


def details(label, cmp, with_examples, max_files):
    bad = [(k, v) for k, v in cmp["tracks"].items() if v["status"] != "exact"]
    hbad = [(k, v) for k, v in cmp["heads"].items() if v["status"] != "exact"]
    if not bad and not hbad:
        return ""
    out = [f"<details><summary>{label}</summary>\n"]
    order = {"diverged": 0, "error": 1, "missing": 2, "close": 3}
    bad.sort(key=lambda kv: (order.get(kv[1]["status"], 9), kv[0]))
    for rel, f in bad[:max_files]:
        st = f["status"]
        if st in ("missing", "error"):
            out.append(f"- `{rel}` ({f['kind']}): **{st}** {f.get('error', '')}")
            continue
        reasons = ", ".join(f"{k}: {v}" for k, v in sorted(f["reasons"].items()))
        out.append(
            f"- `{rel}` ({f['kind']}): **{st}**, {f['diverged'] + f['missing']}"
            f"/{f['particles']} particles diverged, {f['fate_changed']} fates "
            f"changed, max Δxyz {g(max(f['maxd'][k] for k in 'xyz'))}"
            + (f", first divergence in: {reasons}" if reasons else "")
        )
        if with_examples and f["examples"]:
            out.append("\n  ```")
            out.append("  columns: " + " ".join(SHOW))
            for ex in f["examples"]:
                out.append(
                    f"  {ex['particle']}: record {ex['index']} ({ex['reason']}), "
                    f"nrec ref/cand {ex['nrec'][0]}/{ex['nrec'][1]}"
                )
                out.append(f"    prev  {ex['prev'] or '—'}")
                out.append(f"    ref   {ex['ref']}")
                out.append(f"    cand  {ex['cand']}")
            out.append("  ```\n")
    if len(bad) > max_files:
        out.append(f"- … {len(bad) - max_files} more (see report.json)")
    if hbad:
        out.append("\nHead files differing:\n")
        for rel, h in hbad[:max_files]:
            out.append(
                f"- `{rel}`: {h['status']}, max Δh {g(h.get('maxd', 0))}, "
                f"dry mismatches {h.get('dry_mismatch', 0)}"
            )
    out.append("\n</details>\n")
    return "\n".join(out) + "\n"


def cmd_compare(args):
    results = Path(args.results)
    variants = {}
    cross = {}
    for d in sorted(p for p in results.iterdir() if p.is_dir()):
        name = d.name.removeprefix("out-")
        if name.startswith("cross-"):
            cross[name.removeprefix("cross-")] = d
        else:
            variants[name] = d

    def split(v):
        oss, rest = v.split("-", 1)
        tc, fp = rest.split("-", 1)
        return oss, tc, fp

    report = {"tolerance": {"rtol": args.rtol, "atol": args.atol}, "sections": {}}
    md = [
        "# PRT compiler comparison\n",
        f"Tolerance for *close*: |a-b| <= {args.atol} + {args.rtol}·|b|. "
        "*Exact* means bitwise equal (binary track files) or equal at printed "
        "precision (CSV, when no binary track file exists). A particle "
        "*diverges* at the first record where a discrete field (period, step, "
        "layer, cell, zone, status, reason) differs, a coordinate or time "
        "differs beyond tolerance, or the record count differs. *Fate* is the "
        "final (status, layer, cell). Max Δxyz is taken over records before "
        "any divergence.\n",
    ]

    # builds
    md.append("## Builds\n")
    md.append(
        "| variant | compiler | FFLAGS | meson args | pytest (fail/err/total) "
        "| option warnings |\n|---|---|---|---|---|---|\n"
    )
    builds = {}
    for v, d in variants.items():
        env = {}
        for line in read_text(d / "env.txt").splitlines():
            if "=" in line:
                k, val = line.split("=", 1)
                env[k] = val
        ver = next(
            (ln for ln in read_text(d / "compiler-version.txt").splitlines() if ln),
            "?",
        )
        warn = [
            ln.strip()
            for ln in read_text(d / "build.log").splitlines()
            if "unknown option" in ln.lower()
            or "ignoring option" in ln.lower()
            or "#10006" in ln
            or "#10148" in ln
        ]
        warn = sorted(set(warn))
        ju = junit_counts(d / "junit.xml")
        builds[v] = {
            "env": env,
            "compiler": ver,
            "warnings": warn,
            "junit": ju,
            "compile_command": read_text(d / "compile-command.txt"),
        }
        pt = f"{ju['failures']}/{ju['errors']}/{ju['tests']}" if ju else "no junit"
        md.append(
            f"| {v} | {ver[:80]} | `{env.get('FFLAGS', '')}` "
            f"| `{env.get('MESON_ARGS', '')}` | {pt} "
            f"| {'<br>'.join(w[:100] for w in warn[:3]) or '—'} |\n"
        )
    report["builds"] = builds
    failing = [(v, b["junit"]["failed"]) for v, b in builds.items() if b["junit"]]
    failing = [(v, f) for v, f in failing if f]
    if failing:
        md.append("\n<details><summary>Failing tests by variant</summary>\n\n")
        for v, f in failing:
            md.append(f"**{v}**\n\n" + "\n".join(f"- `{t}`" for t in f) + "\n\n")
        md.append("</details>\n")
    md.append("\n<details><summary>Compile commands (one PRT source)</summary>\n\n")
    for v, b in builds.items():
        md.append(f"**{v}**\n```\n{b['compile_command'].strip()}\n```\n")
    md.append("</details>\n\n")

    def section(key, title, blurb, pairs, with_examples, heads=True):
        md.append(f"## {title}\n\n{blurb}\n\n")
        if not pairs:
            md.append("_no data_\n\n")
            return
        md.append(SUMMARY_HEADER)
        sec = {}
        dets = []
        for label, rdir, cdir in pairs:
            cmp = compare_dirs(rdir, cdir, args.rtol, args.atol, 3, heads=heads)
            s = summarize(cmp)
            sec[label] = {"summary": s, **cmp}
            md.append(summary_row(label, s))
            dets.append(details(label, cmp, with_examples, args.max_files))
        md.append("\n")
        md.extend(d for d in dets if d)
        report["sections"][key] = sec

    oses = sorted({split(v)[0] for v in variants})

    # A: ifx vs ifort with the same flags
    pairs = []
    for oss in oses:
        for v in variants:
            o, tc, fp = split(v)
            if o == oss and tc == "ifx" and f"{oss}-ifort-{fp}" in variants:
                pairs.append(
                    (
                        f"{oss} {fp}: ifx vs ifort",
                        variants[f"{oss}-ifort-{fp}"],
                        variants[v],
                    )
                )
    section(
        "A",
        "A. ifx vs ifort, same flags",
        "If `o0-strict` agrees but `default` does not, the difference comes "
        "from optimization or the FP model. If `o0-strict` also differs, "
        "suspect a real bug (uninitialized data, exact real comparisons) or "
        "a difference in language semantics.",
        pairs,
        True,
    )

    # B: everything vs the release baseline (ifort default)
    pairs = []
    for oss in oses:
        ref = f"{oss}-{REFERENCE}"
        if ref not in variants:
            continue
        for v in variants:
            if v != ref and split(v)[0] == oss:
                pairs.append((f"{v} vs {ref}", variants[ref], variants[v]))
    section(
        "B",
        "B. All builds vs the release baseline (ifort, default flags)",
        "Which ifx flag set (if any) reproduces current release results?",
        pairs,
        False,
    )

    # C: PRT-only reruns on baseline GWF output
    pairs = []
    for v, d in cross.items():
        oss = split(v)[0]
        ref = f"{oss}-{REFERENCE}"
        if ref in variants:
            pairs.append((f"PRT {v} on GWF {ref}", variants[ref], d))
    section(
        "C",
        "C. PRT rerun on baseline flow output",
        "Each build reran only the PRT simulations, reading heads and budgets "
        "written by the baseline build. Divergence here comes from PRT itself; "
        "divergence in B that disappears here comes from the flow solution. "
        "Simulations that combine GWF and PRT can't be split this way and show "
        "as missing.",
        pairs,
        True,
        heads=False,
    )

    # D: vs gfortran
    pairs = []
    for oss in oses:
        ref = f"{oss}-{GFORTRAN_REFERENCE}"
        if ref not in variants:
            continue
        for v in variants:
            if v != ref and split(v)[0] == oss:
                pairs.append((f"{v} vs {ref}", variants[ref], variants[v]))
    section(
        "D",
        "D. Builds vs gfortran (default flags)",
        "A third opinion: if ifx agrees with gfortran where ifort doesn't, "
        "ifort shouldn't be assumed correct.",
        pairs,
        False,
    )

    Path(args.out_md).write_text("".join(md))
    Path(args.out_json).write_text(json.dumps(report, indent=1, default=str))
    print(f"wrote {args.out_md} and {args.out_json}")


# ---------------------------------------------------------------------------
# rerun
# ---------------------------------------------------------------------------


def cmd_rerun(args):
    root = Path(args.root)
    mf6 = str(Path(args.mf6).resolve())
    # delete all track output so nothing stale survives a failed rerun
    for rel, kind in find_track_files(root).items():
        p = root / rel
        p.unlink()
        if kind == "bin":
            hdr = Path(str(p) + ".hdr")
            if hdr.exists():
                hdr.unlink()
    log = []
    for nam in sorted(root.rglob("mfsim.nam")):
        text = nam.read_text(errors="replace").lower()
        if "prt6" not in text:
            continue
        if any(m in text for m in ("gwf6", "gwt6", "gwe6")):
            log.append(f"SKIP (combined) {nam.parent.relative_to(root)}")
            continue
        try:
            p = subprocess.run(
                [mf6],
                cwd=nam.parent,
                capture_output=True,
                text=True,
                timeout=args.timeout,
            )
            ok = p.returncode == 0 and "normal termination" in p.stdout.lower()
        except subprocess.TimeoutExpired:
            ok = False
        log.append(f"{'OK  ' if ok else 'FAIL'} {nam.parent.relative_to(root)}")
        if not ok:
            # a partial track file would be misleading
            for rel, kind in find_track_files(nam.parent).items():
                (nam.parent / rel).unlink()
    Path(args.log).write_text("\n".join(log) + "\n")
    print("\n".join(log))
    if args.collect:
        # copy track output (with binary headers) to a separate tree for upload
        dest = Path(args.collect)
        for rel, kind in find_track_files(root).items():
            srcs = [root / rel]
            if kind == "bin":
                srcs.append(root / (rel + ".hdr"))
            for s in srcs:
                t = dest / s.relative_to(root)
                t.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(s, t)
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.log, dest / Path(args.log).name)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("config")
    c.add_argument("--os", required=True)
    c.add_argument("--toolchain", required=True)
    c.add_argument("--fp", required=True)
    c.add_argument("--ifort-version", default="2021.6")
    c.add_argument("--ifx-version", default="2025.2")
    c.add_argument("--gcc-version", default="13")
    c.set_defaults(func=cmd_config)

    b = sub.add_parser("buildinfo")
    b.add_argument("--builddir", required=True)
    b.add_argument("--out", required=True)
    b.set_defaults(func=cmd_buildinfo)

    r = sub.add_parser("rerun")
    r.add_argument("--root", required=True)
    r.add_argument("--mf6", required=True)
    r.add_argument("--log", required=True)
    r.add_argument("--timeout", type=int, default=900)
    r.add_argument("--collect", help="copy track output here after rerunning")
    r.set_defaults(func=cmd_rerun)

    m = sub.add_parser("compare")
    m.add_argument("--results", required=True)
    m.add_argument("--out-md", default="report.md")
    m.add_argument("--out-json", default="report.json")
    m.add_argument("--rtol", type=float, default=1e-8)
    m.add_argument("--atol", type=float, default=1e-8)
    m.add_argument("--max-files", type=int, default=25)
    m.set_defaults(func=cmd_compare)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
