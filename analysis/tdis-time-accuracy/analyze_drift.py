"""Measure TDIS time drift from instrumented MODFLOW 6 runs.

Each run directory holds ``drift_<mode>.csv`` files written by a build of
MODFLOW 6 with TDIS drift logging (``MF6_TDIS_DRIFT``), one per first-step
method (``MF6_TDIS_FIRSTSTEP`` = ``int``, ``pow``, ``hybrid``). Each row is one
attempt at a time step; a retry row replaces the earlier attempt.

References are computed in high precision from the logged double values:

- ``ideal``: the time the TDIS inputs describe, ``T0 + P (r^k - 1) / (r^n - 1)``
  for fixed-stepping periods (``T0`` is the exact sum of earlier ``PERLEN``).
  For ATS periods only the period end is defined.
- ``sum``: the exact sum of the step lengths actually used. ``totim - sum`` is
  the rounding accumulated by ``totim = totim + delt`` alone.

Usage: analyze_drift.py <out_dir> <run_dir>...
"""

import csv
import json
import math
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

getcontext().prec = 80
MODES = ("int", "pow", "hybrid")
EPS = 2.0**-52
DSAME = 100 * EPS
TOLS = {
    "DSAME*t": ("rel", DSAME),
    "BMI 1e-10*t": ("rel", 1e-10),
    "PRT 2.2e-5": ("abs", DSAME * 1e9),
}


def read_steps(path):
    """Accepted steps in order: a retry row replaces the previous attempt."""
    steps, retries = [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            r = {k: v.strip() for k, v in row.items()}
            rec = {
                "kper": int(r["kper"]),
                "kstp": int(r["kstp"]),
                "nstp": int(r["nstp"]),
                "ats": r["ats"] == "1",
                "eop": r["eop"] == "1",
                **{
                    k: float(r[k])
                    for k in (
                        "perlen",
                        "tsmult",
                        "delt",
                        "totim",
                        "pertim",
                        "totimc",
                        "totalsimtime",
                    )
                },
            }
            if r["retry"] == "1":
                retries.append(rec)
                steps[-1] = rec
            else:
                steps.append(rec)
    return steps, retries


def ulp_err(x, ref):
    """Signed error of double x against an exact reference, in ulps of the reference."""
    u = math.ulp(float(ref)) if ref != 0 else math.ulp(0.0)
    return float((Decimal(x) - Decimal(ref)) / Decimal(u))


def analyze(steps):
    """Per-step errors against the ideal and exact-sum references."""
    out = []
    acc = Fraction(0)
    t0 = Decimal(0)  # exact start of current period
    cur = None
    cur_perlen = Decimal(0)
    for s in steps:
        if s["kper"] != cur:
            t0 += cur_perlen
            cur = s["kper"]
            p, n, r = Decimal(s["perlen"]), s["nstp"], Decimal(s["tsmult"])
            rn = r**n if r != 1 else None
            rk = Decimal(1)
            cur_perlen = p
        acc += Fraction(s["delt"])
        rec = {"kper": s["kper"], "kstp": s["kstp"], "eop": s["eop"], "ats": s["ats"]}
        rec["acc_ulp"] = ulp_err(
            s["totim"], Decimal(acc.numerator) / Decimal(acc.denominator)
        )
        ideal = dt_ideal = None
        if not s["ats"]:
            k = s["kstp"]
            if rn is None:
                ideal = t0 + p * k / n
                dt_ideal = p / n
            else:
                rkm1 = rk
                rk = rk * r
                ideal = t0 + p * (rk - 1) / (rn - 1)
                dt_ideal = p * (rk - rkm1) / (rn - 1)
        elif s["eop"]:
            ideal = t0 + p
        if ideal is not None:
            rec["t"] = float(ideal)
            rec["ulp"] = ulp_err(s["totim"], ideal)
            rec["abs"] = float(Decimal(s["totim"]) - ideal)
            rec["rel"] = rec["abs"] / float(ideal) if ideal else 0.0
        if dt_ideal is not None:
            rec["dt_rel"] = float((Decimal(s["delt"]) - dt_ideal) / dt_ideal)
        # end of step as some packages compute it (totimc + delt), vs totim
        rec["end_ulp"] = ulp_err(s["totimc"] + s["delt"], Decimal(s["totim"]))
        out.append(rec)
    return out


def summarize(name, mode, steps, retries, recs):
    withref = [r for r in recs if "ulp" in r]
    eops = [r for r in withref if r["eop"]]
    last = [r for r in recs if r["eop"] and "dt_rel" in r]
    s = {
        "model": name,
        "mode": mode,
        "nper": steps[-1]["kper"],
        "nsteps": len(steps),
        "T": steps[-1]["totalsimtime"],
        "retries": len(retries),
        "retries_eop": sum(r["eop"] for r in retries),
        "retry_end_of_sim": bool(
            retries and retries[-1]["eop"] and retries[-1]["kper"] == steps[-1]["kper"]
        ),
        "max_ulp": max((abs(r["ulp"]) for r in withref), default=0.0),
        "final_ulp": withref[-1]["ulp"] if withref else 0.0,
        "max_acc_ulp": max(abs(r["acc_ulp"]) for r in recs),
        "max_eop_ulp": max((abs(r["ulp"]) for r in eops), default=0.0),
        "eop_mismatch": sum(abs(r["ulp"]) > 0.5 for r in eops),
        "n_eop": len(eops),
        "max_dt_rel": max(
            (abs(r["dt_rel"]) for r in recs if "dt_rel" in r), default=0.0
        ),
        "max_last_dt_rel": max((abs(r["dt_rel"]) for r in last), default=0.0),
        "max_end_ulp": max(abs(r["end_ulp"]) for r in recs),
    }
    for tol, (kind, v) in TOLS.items():
        hit = next(
            (
                r
                for r in withref
                if abs(r["abs"]) > (v * abs(r["t"]) if kind == "rel" else v)
            ),
            None,
        )
        s[f"cross {tol}"] = f"{hit['kper']}/{hit['kstp']}" if hit else ""
    return s


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    rows, series = [], {}
    for d in map(Path, sys.argv[2:]):
        name = (
            d.name if d.name not in ("gwe", "mf6gwt") else f"{d.parent.name}/{d.name}"
        )
        for mode in MODES:
            f = d / f"drift_{mode}.csv"
            if not f.exists() or f.stat().st_size < 200:
                print(f"skip {name} {mode}", file=sys.stderr)
                continue
            steps, retries = read_steps(f)
            recs = analyze(steps)
            rows.append(summarize(name, mode, steps, retries, recs))
            series[f"{name}|{mode}"] = [
                [r.get("t"), r.get("ulp"), r.get("dt_rel"), r["eop"], r["acc_ulp"]]
                for r in recs
            ]
            print(f"done {name} {mode}", file=sys.stderr)
    with open(out / "drift_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(out / "drift_series.json", "w") as f:
        json.dump(series, f)


if __name__ == "__main__":
    main()
