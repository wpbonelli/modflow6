"""Figures comparing ways of computing the TDIS first time step length.

Methods, for a geometric period (tsmult != 1):

  current  integer exponent, tsmult**nstp, by repeated squaring (gfortran's
           __powidf2), as in tdis_set_delt today
  #3204    real exponent, tsmult**real(nstp, DP), i.e. libm pow
  expm1    perlen * (tsmult - 1) / expm1(nstp * log1p(tsmult - 1))
  hybrid   expm1 form for |nstp * log1p(tsmult - 1)| < 1, #3204 form otherwise

Errors are against exact rational arithmetic for the same double inputs, in ulps
of the exact result. Writes three PNGs into figures/ next to this script:

  first_step_map.png     first-step error over tsmult and nstp, per method
  first_step_slice.png   the same at nstp = 10 and 100
  period_errors.png      max totim error within the period and last-step delt
                         error, replaying tdis_set_delt for each method

Needs numpy and matplotlib. Run with: python3 -I plot_errors.py
Takes a few seconds.
"""

import math
from fractions import Fraction as F
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

OUT = Path(__file__).resolve().parent / "figures"
METHODS = ["current", "#3204", "expm1", "hybrid"]
TITLES = {
    "current": "current (integer **)",
    "#3204": "#3204 (real **)",
    "expm1": "expm1/log1p only",
    "hybrid": "hybrid (both)",
}
COLORS = {
    "current": "#888888",
    "#3204": "#d1495b",
    "expm1": "#2e86ab",
    "hybrid": "#1b1b1b",
}
FLOOR = 0.1  # plot errors below this (including exact results) at the floor


def powi(r, n):
    """r**n by repeated squaring, rounding at each multiply (__powidf2)."""
    res, b, m = 1.0, r, n
    while m:
        if m & 1:
            res *= b
        m >>= 1
        if m:
            b *= b
    return res


def first_steps(p, n, r):
    y = n * math.log1p(r - 1)
    f = {
        "current": p * (1 - r) / (1 - powi(r, n)),
        "#3204": p * (1 - r) / (1 - r ** float(n)),
        "expm1": p * (r - 1) / math.expm1(y),
    }
    f["hybrid"] = f["expm1"] if abs(y) < 1 else f["#3204"]
    return f


def ulps(x, ex):
    return max(FLOOR, float(abs(F(x) - ex) / F(math.ulp(float(ex)))))


def first_step_errors(p, n, r):
    rn = F(r) ** n
    ex = F(p) * (1 - F(r)) / (1 - rn)
    return {k: ulps(v, ex) for k, v in first_steps(p, n, r).items()}


def period_errors(p, n, r, d):
    """Replay one period from t = 0 with first step d; return max totim error
    (ulp) over the period and the last-step delt error (relative)."""
    fr = F(r)
    d1 = F(p) * (1 - fr) / (1 - fr**n)
    t, max_t = 0.0, FLOOR
    ex_d, ex_t = d1, F(0)
    for k in range(1, n + 1):
        if k > 1:
            d = p - t if k == n else r * d
            ex_d *= fr
        t = t + d
        ex_t += ex_d
        if k < n:
            max_t = max(max_t, ulps(t, ex_t))
    rel = max(1e-17, float(abs(F(d) - ex_d) / ex_d))
    return max_t, rel


def sides(h):
    """tsmult values above and below 1 for offsets h."""
    return {">1": [1 + x for x in h], "<1": [1 - x for x in h]}


def plot_map():
    h = np.logspace(-8, math.log10(0.5), 70)
    nstp = np.unique(np.round(np.logspace(math.log10(2), 3, 45)).astype(int))
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=True, sharey=True)
    norm = LogNorm(vmin=FLOOR, vmax=1e7)
    for row, (side, rs) in enumerate(sides(h).items()):
        err = {m: np.zeros((len(nstp), len(h))) for m in METHODS}
        ys = np.zeros((len(nstp), len(h)))
        for i, n in enumerate(nstp):
            for j, r in enumerate(rs):
                e = first_step_errors(1.0, int(n), r)
                for m in METHODS:
                    err[m][i, j] = e[m]
                ys[i, j] = abs(n * math.log1p(r - 1))
        for col, m in enumerate(METHODS):
            ax = axes[row, col]
            pc = ax.pcolormesh(h, nstp, err[m], norm=norm, cmap="magma_r")
            ax.contour(h, nstp, ys, levels=[1.0], colors="#2e86ab", linewidths=1)
            ax.set_xscale("log")
            ax.set_yscale("log")
            worst = err[m].max()
            ax.set_title(
                f"{TITLES[m]}\ntsmult {side}, worst {worst:.3g} ulp", fontsize=9
            )
            if row == 1:
                ax.set_xlabel("|tsmult - 1|")
            if col == 0:
                ax.set_ylabel("nstp")
    cb = fig.colorbar(pc, ax=axes, shrink=0.8)
    cb.set_label("first-step delt error (ulp)")
    fig.suptitle(
        "First step length error, perlen = 1. Blue line: |nstp * log1p(tsmult - 1)| = 1"
        " (hybrid switches there)",
        fontsize=10,
    )
    fig.savefig(OUT / "first_step_map.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_slice():
    h = np.logspace(-8, math.log10(0.5), 250)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True)
    for row, n in enumerate([10, 100]):
        for col, (side, rs) in enumerate(sides(h).items()):
            ax = axes[row, col]
            errs = [first_step_errors(1.0, n, r) for r in rs]
            for m in METHODS:
                ax.plot(
                    h,
                    [e[m] for e in errs],
                    ".",
                    ms=3,
                    color=COLORS[m],
                    label=TITLES[m],
                )
            hs = math.expm1(1 / n) if side == ">1" else -math.expm1(-1 / n)
            ax.axvline(hs, color="#2e86ab", lw=0.8, ls="--")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_title(f"nstp = {n}, tsmult {side}", fontsize=9)
            if row == 1:
                ax.set_xlabel("|tsmult - 1|")
            if col == 0:
                ax.set_ylabel("first-step delt error (ulp)")
    axes[0, 0].legend(fontsize=8, markerscale=3)
    fig.suptitle(
        f"First step length error, perlen = 1 (errors below {FLOOR} ulp, including"
        " exact, drawn at the floor). Dashed: hybrid switch",
        fontsize=10,
    )
    fig.savefig(OUT / "first_step_slice.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_period():
    n, p = 100, 10.0
    h = np.logspace(-8, math.log10(0.2), 60)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for col, (side, rs) in enumerate(sides(h).items()):
        res = {m: [] for m in METHODS}
        for r in rs:
            f = first_steps(p, n, r)
            for m in METHODS:
                res[m].append(period_errors(p, n, r, f[m]))
        for m in METHODS:
            kw = dict(color=COLORS[m], label=TITLES[m], lw=1, marker=".", ms=3)
            axes[0, col].plot(h, [x[0] for x in res[m]], **kw)
            axes[1, col].plot(h, [x[1] for x in res[m]], **kw)
        axes[0, col].set_title(f"tsmult {side}", fontsize=9)
        axes[1, col].set_xlabel("|tsmult - 1|")
        for ax in axes[:, col]:
            ax.set_xscale("log")
            ax.set_yscale("log")
    axes[0, 0].set_ylabel("max totim error in period (ulp)")
    axes[1, 0].set_ylabel("last-step delt error (relative)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        f"Within one period (perlen = {p:g}, nstp = {n}), replaying tdis_set_delt"
        " with each first-step method",
        fontsize=10,
    )
    fig.savefig(OUT / "period_errors.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    plot_map()
    plot_slice()
    plot_period()
    print(f"wrote figures to {OUT}")
