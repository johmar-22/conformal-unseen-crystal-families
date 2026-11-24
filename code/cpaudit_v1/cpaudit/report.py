"""Aggregate Phase 1 runs into tables, a pre-registered endpoint check, and one figure."""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

from .config import atomic_write_csv, atomic_write_json, atomic_write_text, get_paths

LATTICE_ORDER = ["C2/m", "Pna21", "R-3c", "P63/mmc", "Ia-3", "Fd-3m"]
METHOD_ORDER = ["cvplus_lattice", "cvplus_random", "split_global", "split_mondrian",
                "split_normalized_knn", "cqr"]
METHOD_LABEL = {"cvplus_lattice": "CV+ (lattice folds, submitted)", "cvplus_random": "CV+ (random folds)",
                "split_global": "Split (global)", "split_mondrian": "Split (Mondrian)",
                "split_normalized_knn": "Split (normalized kNN)", "cqr": "CQR"}


def load_runs(paths) -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(paths.runs, "r*_*.csv")))
    if not files:
        raise RuntimeError(f"No Phase 1 results in {paths.runs}")
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    keys = ["target", "tag", "method", "alpha", "group"]
    rows = []
    for k, g in df.groupby(keys, sort=False):
        a = float(k[3])
        c = g["coverage"]
        rows.append(dict(zip(keys, k), **dict(
            n_repeats=len(g), cov_mean=c.mean(), cov_median=c.median(), cov_p10=c.quantile(0.1),
            cov_p90=c.quantile(0.9), frac_below_nominal=float((c < 1 - a - 1e-12).mean()),
            frac_below_expected_if_valid=g["p_below_if_valid"].mean(),
            width_median=g["width_median"].median(),
            width_mean=g["width_mean"].replace(np.inf, np.nan).mean(),
            interval_score=g["interval_score"].replace(np.inf, np.nan).mean(), mae=g["mae"].mean(),
            spearman_width_err=g["spearman_width_abserr"].mean(), n_test_mean=g["n"].mean(),
            n_test_groups_mean=g["n_test_comp_groups"].mean())))
    return pd.DataFrame(rows)


def endpoint(summary: pd.DataFrame, alpha: float) -> dict:
    """Pre-registered primary endpoint and decision (see PREREGISTRATION.md)."""
    res = {}
    s = summary[(summary.tag == "hpo_once") & (summary.alpha.round(4) == round(alpha, 4))]
    for tgt in ("bandgap", "formation"):
        row = s[(s.target == tgt) & (s.method == "cvplus_random") & (s.group == "Fd-3m")]
        if row.empty:
            continue
        v = float(row.cov_mean.iloc[0])
        decision = "H1a supported" if v < 0.85 else ("H1a rejected (gate G1)" if v >= 0.88 else "inconclusive")
        mon = s[(s.target == tgt) & (s.method == "split_mondrian") & (s.group.isin(LATTICE_ORDER))]
        h1b = bool((mon.cov_mean >= 0.88).all() and
                   ((mon.frac_below_nominal - mon.frac_below_expected_if_valid) <= 0.15).all()) if not mon.empty else None
        res[tgt] = dict(primary_endpoint_fd3m_cov_mean=v, decision_H1a=decision, H1b_supported=h1b)
    return res


def hpo_sensitivity(df: pd.DataFrame, alpha: float) -> pd.DataFrame:
    d = df[(df.alpha.round(4) == round(alpha, 4)) & (df.group.isin(LATTICE_ORDER + ["ALL"]))]
    own = d[d.tag == "hpo_own"]
    if own.empty:
        return pd.DataFrame()
    key = ["repeat", "target", "method", "group"]
    m = own.merge(d[d.tag == "hpo_once"], on=key, suffixes=("_own", "_once"))
    m["delta_coverage"] = m["coverage_own"] - m["coverage_once"]
    m["delta_width_median"] = m["width_median_own"] - m["width_median_once"]
    return m[key + ["coverage_once", "coverage_own", "delta_coverage", "width_median_once",
                    "width_median_own", "delta_width_median"]]


def figure(df: pd.DataFrame, paths, alpha: float):
    """Per-lattice coverage across repeats for four calibrations, both targets."""
    sys.path.insert(0, paths.drive_root)
    try:
        import figstyle as fs   # the thesis figure standard already in the Drive folder
        fig, axes = fs.fig_double(0.75, 2, 1)
        style = lambda i: fs.series(i, line=False)
        panel = fs.panel
        save = lambda f, name: fs.save_fig(f, name, outdir=paths.report)
    except Exception as e:  # fall back to plain matplotlib if figstyle is unavailable
        print(f"  figstyle not usable ({e}); using matplotlib defaults")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(5.9, 5.3), constrained_layout=True)
        style = lambda i: dict(marker="osD^"[i % 4], linestyle="none")
        panel = lambda ax, l: ax.text(0, 1.02, f"({l})", transform=ax.transAxes)
        save = lambda f, name: f.savefig(os.path.join(paths.report, name + ".png"), dpi=600)
    methods = ["cvplus_lattice", "cvplus_random", "split_global", "split_mondrian"]
    d = df[(df.tag == "hpo_once") & (df.alpha.round(4) == round(alpha, 4))]
    for ax, tgt, letter in zip(axes, ("bandgap", "formation"), "ab"):
        x0 = np.arange(len(LATTICE_ORDER))
        for i, m in enumerate(methods):
            g = d[(d.target == tgt) & (d.method == m)]
            off = (i - 1.5) * 0.18
            for j, lat in enumerate(LATTICE_ORDER):
                c = g[g.group == lat].coverage.to_numpy()
                if c.size == 0:
                    continue
                med, lo, hi = np.median(c), np.quantile(c, 0.1), np.quantile(c, 0.9)
                st = style(i)
                ax.errorbar(x0[j] + off, med, yerr=[[med - lo], [hi - med]], label=METHOD_LABEL[m] if j == 0 else None,
                            capsize=1.5, **st)
        ax.axhline(1 - alpha, color="0.4", lw=0.7, ls="--")
        ax.set_xticks(x0, LATTICE_ORDER)
        ax.set_ylabel(f"{'Band gap' if tgt == 'bandgap' else 'Formation energy'} coverage")
        ax.set_ylim(min(0.3, ax.get_ylim()[0]), 1.02)
        panel(ax, letter)
    placer = getattr(sys.modules.get("figstyle"), "place_legend", None)
    if placer is not None:
        placer(axes[0], prefer="lower left", ncol=2)   # moves the legend off the data if needed
    else:
        axes[0].legend(ncol=2, loc="lower left")
    fig.text(0.5, -0.01, "Markers: median over splits; bars: 10th-90th percentile; dashed: nominal",
             ha="center", va="top", fontsize=7)
    save(fig, "fig_phase1_lattice_coverage")


def run_report(cfg):
    paths = get_paths(cfg)
    df = load_runs(paths)
    a = float(cfg["alpha_primary"])
    summary = aggregate(df)
    atomic_write_csv(os.path.join(paths.report, "phase1_all_runs.csv"), df)
    atomic_write_csv(os.path.join(paths.report, "phase1_summary.csv"), summary)
    ep = endpoint(summary, a)
    atomic_write_json(os.path.join(paths.report, "phase1_endpoint.json"), ep)
    sens = hpo_sensitivity(df, a)
    if not sens.empty:
        atomic_write_csv(os.path.join(paths.report, "phase1_hpo_sensitivity.csv"), sens)
    try:
        figure(df, paths, a)
    except Exception as e:
        print(f"  figure skipped: {e}")

    # readable digest
    lines = [f"# Phase 1 summary (alpha = {a})", "",
             f"Repeats found: {df.repeat.nunique()}  |  units: {len(df.groupby(['repeat','target','tag']))}", ""]
    for tgt in ("bandgap", "formation"):
        s = summary[(summary.target == tgt) & (summary.tag == "hpo_once") & (summary.alpha.round(4) == round(a, 4))
                    & (summary.group.isin(LATTICE_ORDER + ["ALL"]))]
        if s.empty:
            continue
        piv = s.pivot(index="group", columns="method", values="cov_mean").reindex(["ALL"] + LATTICE_ORDER)
        piv = piv[[m for m in METHOD_ORDER if m in piv.columns]]
        wid = s.pivot(index="group", columns="method", values="width_median").reindex(["ALL"] + LATTICE_ORDER)
        wid = wid[[m for m in METHOD_ORDER if m in wid.columns]]
        lines += [f"## {tgt}: mean coverage over repeats", "", piv.round(3).to_markdown(), "",
                  f"## {tgt}: median interval width", "", wid.round(3).to_markdown(), ""]
    lines += ["## Pre-registered endpoint", "", "```", json.dumps(ep, indent=2), "```", ""]
    if not sens.empty:
        lines += ["## Hyperparameter sensitivity (own vs split-0 tuning)", "",
                  sens.groupby(["target", "method"])["delta_coverage"].agg(["mean", "min", "max"]).round(3).to_markdown(), ""]
    atomic_write_text(os.path.join(paths.report, "PHASE1_SUMMARY.md"), "\n".join(lines))
    print("\n".join(lines))
    print(f"\nReport written to {paths.report}")
