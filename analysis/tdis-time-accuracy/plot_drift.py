"""Plot TDIS drift measured by analyze_drift.py.

Usage: plot_drift.py <analysis_out_dir> <figure_dir>
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLORS = {"int": "#888888", "pow": "#d1495b", "hybrid": "#1b1b1b"}
LABELS = {
    "int": "current (integer **)",
    "pow": "#3204 (real **)",
    "hybrid": "hybrid (expm1/log1p)",
}
MARKERS = {"int": "o", "pow": "s", "hybrid": "D"}
FLOOR = 0.1

SERIES = [
    "uniform_1e5_tenth",
    "test1200_gwtbuy-goswami",
    "periods_2000x0.1",
    "near1_r1.00001_n10000",
    "ex-gwe-prt/gwe",
    "test205_gwtbuy-henrytidal",
]


def load(d):
    with open(d / "drift_summary.csv") as f:
        rows = list(csv.DictReader(f))
    with open(d / "drift_series.json") as f:
        series = json.load(f)
    return rows, series


def plot_summary(rows, fig_dir):
    models = list(dict.fromkeys(r["model"] for r in rows))
    by = {(r["model"], r["mode"]): r for r in rows}
    fig, axes = plt.subplots(1, 3, figsize=(13, 0.32 * len(models) + 1.5), sharey=True)
    y = np.arange(len(models))
    panels = [
        ("max_ulp", "max |totim error| vs inputs (ulp)"),
        ("max_eop_ulp", "max |totim error| at period ends (ulp)"),
        ("max_last_dt_rel", "max last-step delt relative error (exact at 1e-17)"),
    ]
    for ax, (key, title) in zip(axes, panels):
        for i, mode in enumerate(("int", "pow", "hybrid")):
            x = [
                np.nan
                if (m, mode) not in by
                or (key == "max_last_dt_rel" and m.startswith("ats_"))
                else max(float(by[(m, mode)][key]), FLOOR if "ulp" in key else 1e-17)
                for m in models
            ]
            ax.scatter(
                x,
                y + (i - 1) * 0.22,
                s=22,
                marker=MARKERS[mode],
                color=COLORS[mode],
                label=LABELS[mode],
                zorder=3,
            )
        ax.set_xscale("log")
        ax.set_title(title, fontsize=10)
        ax.grid(axis="x", color="#e5e5e5", lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].axvspan(100, 200, color="#2a78d6", alpha=0.12, lw=0)
    axes[1].axvspan(100, 200, color="#2a78d6", alpha=0.12, lw=0)
    axes[2].axvline(1e-10, color="#2a78d6", lw=0.8, ls="--")
    axes[2].text(1e-10, -0.8, " 1e-10", fontsize=8, color="#184f95")
    axes[0].set_yticks(y, models, fontsize=8)
    axes[0].invert_yaxis()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        fontsize=9,
        frameon=False,
        bbox_to_anchor=(0.5, 0.975),
    )
    fig.suptitle(
        f"TDIS drift by model (values below {FLOOR} ulp shown at {FLOOR}; "
        "shaded: DSAME·t)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(fig_dir / "drift_summary.png", dpi=150)
    plt.close(fig)


def plot_series(series, fig_dir):
    names = [n for n in SERIES if f"{n}|int" in series]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for ax, name in zip(axes.flat, names):
        for mode in ("int", "pow", "hybrid"):
            pts = [p for p in series.get(f"{name}|{mode}", []) if p[1] is not None]
            if not pts:
                continue
            t = np.array([p[0] for p in pts])
            e = np.abs(np.array([p[1] for p in pts]))
            ax.plot(
                t,
                np.maximum(e, FLOOR),
                color=COLORS[mode],
                lw=1.2 if mode != "int" else 2.0,
                label=LABELS[mode],
            )
        ax.set_yscale("log")
        ax.set_ylim(FLOOR / 2, None)
        ax.axhspan(100, 200, color="#2a78d6", alpha=0.12, lw=0)
        ax.set_title(name, fontsize=10)
        ax.set_xlabel("simulation time", fontsize=9)
        ax.grid(color="#eeeeee", lw=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("|totim error| (ulp)")
    for ax in axes.flat[len(names) :]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=8, frameon=False)
    fig.suptitle(
        "totim error against the TDIS inputs, per step "
        "(shaded: DSAME·t; where methods agree, lines coincide)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(fig_dir / "drift_series.png", dpi=150)
    plt.close(fig)


def main():
    d, fig_dir = Path(sys.argv[1]), Path(sys.argv[2])
    fig_dir.mkdir(parents=True, exist_ok=True)
    rows, series = load(d)
    plot_summary(rows, fig_dir)
    plot_series(series, fig_dir)


if __name__ == "__main__":
    main()
