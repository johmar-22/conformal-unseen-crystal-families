"""Paper figures from the saved Phase 1-4 outputs (read only; nothing is refitted).

Style: figures4papers house style (see cpaudit/pubstyle.py for attribution and licence).
Every figure is drawn only from files already written by the pipeline:

  Fig. 1  phase01_v1/report/phase1_all_runs.csv                     in-distribution coverage
  Fig. 2  phase02_v1/post1/cluster_mondrian_{summary,by_run}.csv    group-level Mondrian
  Fig. 3  phase02_v1/report/phase2_{loso,l2o}_all.csv + loso/*.npz  unseen-lattice coverage
  Fig. 4  phase02_v1/loso/*.npz + report csv                        mechanism (offset, transfer)
  Fig. 5  phase02_v1/overlap/*                                      no covariate overlap
  Fig. 6  phase04_v1/report/*                                       failure detection (Phase 4)
  Fig. S1 calibration-tail composition; Fig. S2 formation-energy transfer; Fig. S3 detection radar

A figure whose inputs are missing is skipped with a message, so the command can be run at any
stage. Output: <drive_root>/<figures_dir>/ as PNG (600 dpi), PDF and SVG, plus a manifest with
the SHA-256 of every input file used.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from .config import sha256_file
from .pubstyle import (DOUBLE_COL_IN, LATTICE_COLORS, LATTICE_ORDER, LATTICE_TEX, PALETTE,
                       apply_publication_style, coverage_cmap, coverage_norm, finalize_figure,
                       new_figure, panel_label, sequential_red)

import matplotlib.pyplot as plt  # noqa: E402  (pubstyle selects the Agg backend first)
from matplotlib.colors import LogNorm, TwoSlopeNorm, LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

SG2LAT = {12: "C2/m", 33: "Pna21", 167: "R-3c", 194: "P63/mmc", 206: "Ia-3", 227: "Fd-3m"}
LAT2SG = {v: k for k, v in SG2LAT.items()}
TARGETS = {"bandgap": ("Band gap", "eV"), "formation": ("Formation energy", "eV/atom")}
ALPHA = 0.10
NOMINAL = 1 - ALPHA

P1_METHODS = [  # (key, label, colour)
    ("cvplus_random", "CV+, random folds", PALETTE["red_strong"]),
    ("split_global", "Split, global", PALETTE["red_2"]),
    ("split_normalized_knn", "Split, normalized (kNN)", PALETTE["magenta"]),
    ("cqr", "CQR", PALETTE["teal"]),
    ("split_mondrian", "Mondrian (per lattice)", PALETTE["blue_main"]),
    ("cvplus_lattice", "CV+, lattice folds (submitted)", PALETTE["grey_dark"]),
]
DETECT_SCORES = [  # (key, label, colour)
    ("width", "CV+ interval width", PALETTE["red_strong"]),
    ("knn", "kNN distance", PALETTE["blue_main"]),
    ("kde", "KDE dissimilarity", PALETTE["teal"]),
    ("committee", "Committee SD (lattice folds)", PALETTE["violet"]),
    ("subsample", "Subsample-ensemble SD", PALETTE["magenta"]),
    ("classifier", "Domain-classifier probability", PALETTE["grey_mid"]),
]


# --------------------------------------------------------------------------------------- helpers
class Inputs:
    """Resolves input files and records the SHA-256 of each one that a figure reads."""

    def __init__(self, phase1: str, phase2: str, phase4: str):
        self.dirs = dict(phase1=phase1, phase2=phase2, phase4=phase4)
        self.used: dict[str, str] = {}

    def path(self, phase: str, *parts) -> str:
        return os.path.join(self.dirs[phase], *parts)

    def exists(self, phase: str, *parts) -> bool:
        return os.path.exists(self.path(phase, *parts))

    def csv(self, phase: str, *parts) -> pd.DataFrame:
        p = self.path(phase, *parts)
        self.used[p] = sha256_file(p)
        return pd.read_csv(p)

    def npz(self, phase: str, *parts):
        p = self.path(phase, *parts)
        self.used[p] = sha256_file(p)
        return np.load(p, allow_pickle=False)

    def json(self, phase: str, *parts):
        p = self.path(phase, *parts)
        self.used[p] = sha256_file(p)
        with open(p) as f:
            return json.load(f)


def _missing(inp: Inputs, needs) -> list[str]:
    return [inp.path(ph, *pt) for ph, *pt in needs if not inp.exists(ph, *pt)]


def _text_color(rgba) -> str:
    r, g, b = rgba[:3]
    return "white" if (0.299 * r + 0.587 * g + 0.114 * b) < 0.5 else PALETTE["ink"]


def _heatmap(ax, values, row_labels, col_labels, cmap, norm, texts=None, na_text="n/a", fontsize=7.0):
    """Annotated heatmap with white cell borders (figures4papers composition-heatmap style)."""
    nr, nc = values.shape
    for i in range(nr):
        for j in range(nc):
            v = values[i, j]
            if np.isnan(v):
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor="#F2F2F2", edgecolor="white", lw=1.6,
                                       hatch="////", zorder=1))
                ax.text(j + 0.5, i + 0.5, na_text, ha="center", va="center", fontsize=fontsize - 0.5,
                        color=PALETTE["grey_mid"], zorder=3)
                continue
            rgba = cmap(norm(v))
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=rgba, edgecolor="white", lw=1.6, zorder=1))
            if texts is not None and texts[i][j]:
                ax.text(j + 0.5, i + 0.5, texts[i][j], ha="center", va="center", fontsize=fontsize,
                        color=_text_color(rgba), zorder=3, linespacing=1.15)
    ax.set_xlim(0, nc)
    ax.set_ylim(nr, 0)
    ax.set_xticks(np.arange(nc) + 0.5)
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(nr) + 0.5)
    ax.set_yticklabels(row_labels)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)


def _lat_labels(lats):
    return [LATTICE_TEX[l] for l in lats]


def _loso_npz(inp: Inputs, target: str, lat: str, design: str = "loso"):
    name = f"{target}_{design}{LAT2SG[lat]}.npz"
    return inp.npz("phase2", "loso", name) if inp.exists("phase2", "loso", name) else None


# --------------------------------------------------------------------------------------- Fig. 1
def fig1_indistribution(inp: Inputs, out: str):
    miss = _missing(inp, [("phase1", "report", "phase1_all_runs.csv")])
    if miss:
        return None, f"missing {miss}"
    d = inp.csv("phase1", "report", "phase1_all_runs.csv")
    d = d[(d.tag == "hpo_once") & np.isclose(d.alpha, ALPHA)]
    groups = ["ALL"] + LATTICE_ORDER
    glabels = ["All"] + _lat_labels(LATTICE_ORDER)
    strip = [m for m in P1_METHODS if m[0] in ("cvplus_random", "split_mondrian", "cvplus_lattice")]
    rng = np.random.default_rng(7)

    fig = new_figure(DOUBLE_COL_IN, 5.4, wspace=0.06, hspace=0.08)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1.0])
    for r, (target, (tname, unit)) in enumerate(TARGETS.items()):
        ax = fig.add_subplot(gs[r, 0])
        D = d[d.target == target]
        fdx = len(groups) - 1
        ax.axvspan(fdx - 0.46, fdx + 0.46, color=PALETTE["red_1"], alpha=0.45, lw=0, zorder=0)
        ax.axhline(NOMINAL, color=PALETTE["ink"], lw=1.0, ls=(0, (4, 3)), zorder=1)
        lo_all = 1.0
        for k, (m, lab, col) in enumerate(strip):
            off = (k - 1) * 0.27
            for gi, g in enumerate(groups):
                v = D[(D.method == m) & (D.group == g)].coverage.to_numpy()
                if v.size == 0:
                    continue
                lo_all = min(lo_all, v.min())
                x = gi + off + rng.uniform(-0.06, 0.06, v.size)
                ax.scatter(x, v, s=5, color=col, alpha=0.35, lw=0, zorder=2)
                q = np.percentile(v, [2.5, 97.5])
                ax.plot([gi + off] * 2, q, color=col, lw=1.2, zorder=3, solid_capstyle="round")
                ax.scatter([gi + off], [v.mean()], s=24, marker="D", facecolor="white", edgecolor=col, lw=1.3,
                           zorder=4)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(glabels)
        ax.set_xlim(-0.5, len(groups) - 0.5)
        ax.set_ylim(max(0.0, np.floor((lo_all - 0.03) * 20) / 20), 1.015)
        ax.set_ylabel("Coverage (per split)")
        ax.text(fdx, ax.get_ylim()[0] + 0.012, "failing\nfamily", ha="center", va="bottom", fontsize=6.5,
                color=PALETTE["red_strong"])
        panel_label(ax, "ac"[r], f"{tname}: coverage by lattice, 25 splits")

        ax2 = fig.add_subplot(gs[r, 1])
        ax2.axhline(NOMINAL, color=PALETTE["ink"], lw=1.0, ls=(0, (4, 3)), zorder=1)
        fd_cov = []
        for m, lab, col in P1_METHODS:
            M = D[D.method == m]
            if M.empty:
                continue
            per = M.groupby("group").agg(cov=("coverage", "mean"), w=("width_median", "mean"))
            oth = per.loc[[g for g in LATTICE_ORDER if g != "Fd-3m"]]
            x0, y0 = oth.w.mean(), oth["cov"].mean()
            x1, y1 = per.loc["Fd-3m", "w"], per.loc["Fd-3m", "cov"]
            fd_cov.append(y1)
            ax2.annotate("", xy=(x1, y1), xytext=(x0, y0),
                         arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3, alpha=0.85,
                                         shrinkA=3.5, shrinkB=3.5, mutation_scale=8))
            ax2.scatter([x0], [y0], s=34, facecolor="white", edgecolor=col, lw=1.6, zorder=4)
            ax2.scatter([x1], [y1], s=40, color=col, edgecolor="white", lw=0.6, zorder=5)
        ax2.set_xlabel(f"Median interval width ({unit})")
        ax2.set_ylabel("Mean coverage")
        ax2.set_ylim(max(0.0, np.floor((min(fd_cov) - 0.04) * 20) / 20), 1.015)
        panel_label(ax2, "bd"[r], "Other five lattices → " + LATTICE_TEX["Fd-3m"])

    h = [Line2D([], [], color=c, marker="o", lw=0, ms=5, label=l) for _, l, c in P1_METHODS]
    h += [Line2D([], [], color=PALETTE["grey_dark"], marker="D", mfc="white", lw=1.3, ms=4,
                 label="mean and 95 % range over splits"),
          Line2D([], [], color=PALETTE["grey_dark"], marker="o", mfc="white", lw=0, ms=5,
                 label="mean of the other five lattices"),
          Line2D([], [], color=PALETTE["grey_dark"], marker="o", lw=0, ms=5, label=LATTICE_TEX["Fd-3m"]),
          Line2D([], [], color=PALETTE["ink"], lw=1.0, ls=(0, (4, 3)), label=r"nominal $1-\alpha$ = 0.90")]
    fig.legend(handles=h, loc="outside lower center", ncol=4, columnspacing=1.2, handletextpad=0.4, fontsize=7)
    return finalize_figure(fig, os.path.join(out, "fig1_indistribution_coverage")), None


# --------------------------------------------------------------------------------------- Fig. 2
def fig2_group_mondrian(inp: Inputs, out: str):
    miss = _missing(inp, [("phase2", "post1", "cluster_mondrian_summary.csv"),
                          ("phase2", "post1", "cluster_mondrian_by_run.csv")])
    if miss:
        return None, f"missing {miss}"
    S = inp.csv("phase2", "post1", "cluster_mondrian_summary.csv")
    R = inp.csv("phase2", "post1", "cluster_mondrian_by_run.csv")
    S = S[np.isclose(S.alpha, ALPHA)]
    R = R[np.isclose(R.alpha, ALPHA)]
    lats = LATTICE_ORDER
    ypos = {l: i for i, l in enumerate(lats)}

    fig = new_figure(DOUBLE_COL_IN, 3.05, wspace=0.03)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.8])
    for c, (target, (tname, _)) in enumerate(TARGETS.items()):
        ax = fig.add_subplot(gs[0, c])
        T = S[S.target == target].set_index("lattice")
        for l in lats:
            y = ypos[l]
            n = T.loc[l, "n_cal_groups"]
            ax.add_patch(Rectangle((NOMINAL, y - 0.32), 1 / (n + 1), 0.64, facecolor=PALETTE["green_1"],
                                   edgecolor="none", zorder=0))
            pts = [T.loc[l, "pooled_cov_point"], T.loc[l, "pooled_cov_group"], T.loc[l, "single_sub_cov_group"]]
            ax.plot([min(pts), max(pts)], [y, y], color=PALETTE["neutral"], lw=3.2, solid_capstyle="round",
                    zorder=1)
            ax.scatter(pts[0], y, s=34, facecolor="white", edgecolor=PALETTE["red_2"], lw=1.5, zorder=3)
            ax.scatter(pts[1], y, s=36, color=PALETTE["red_strong"], zorder=3)
            ax.scatter(pts[2], y, s=44, marker="D", color=PALETTE["blue_main"], zorder=4)
            fi = T.loc[l, "single_sub_frac_infinite"]
            if fi > 0:
                ax.text(max(pts) + 0.005, y - 0.18, f"∞ in {fi:.0%}", va="bottom", fontsize=6.3,
                        color=PALETTE["blue_main"])
        ax.axvline(NOMINAL, color=PALETTE["ink"], lw=1.0, ls=(0, (4, 3)), zorder=2)
        ax.set_yticks(range(len(lats)))
        ax.set_yticklabels(_lat_labels(lats) if c == 0 else [])
        ax.set_ylim(len(lats) - 0.5, -0.5)
        ax.set_xlim(0.84, 0.985)
        ax.set_xlabel("Coverage (mean over 25 splits)")
        panel_label(ax, "ab"[c], f"{tname}: Mondrian at " + r"$\alpha$ = 0.1")

    ax = fig.add_subplot(gs[0, 2])
    for k, (target, (tname, _)) in enumerate(TARGETS.items()):
        col = [PALETTE["blue_main"], PALETTE["teal"]][k]
        for l in lats:
            G = R[(R.target == target) & (R.lattice == l) & np.isfinite(R.single_sub_width_median)]
            if G.empty:
                continue
            ratio = (G.single_sub_width_median / G.pooled_width).to_numpy()
            q1, med, q3 = np.percentile(ratio, [25, 50, 75])
            y = ypos[l] + (k - 0.5) * 0.3
            ax.plot([q1, q3], [y, y], color=col, lw=1.4, alpha=0.7, zorder=2, solid_capstyle="round")
            ax.scatter(med, y, s=30, marker="os"[k], color=col, zorder=3)
    ax.axvline(1, color=PALETTE["ink"], lw=1.0)
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.8, 1, 1.25, 1.6, 2])
    ax.set_xticklabels(["0.8", "1", "1.25", "1.6", "2"])
    ax.minorticks_off()
    ax.set_yticks(range(len(lats)))
    ax.set_yticklabels([])
    ax.set_ylim(len(lats) - 0.5, -0.5)
    ax.set_xlabel("Width ratio, one-per-group / pooled")
    panel_label(ax, "c", "Price of validity (finite draws)")

    h = [Line2D([], [], marker="o", mfc="white", mec=PALETTE["red_2"], lw=0, ms=5,
                label="pooled calibration, per-structure coverage"),
         Line2D([], [], marker="o", color=PALETTE["red_strong"], lw=0, ms=5,
                label="pooled calibration, group-weighted coverage"),
         Line2D([], [], marker="D", color=PALETTE["blue_main"], lw=0, ms=5,
                label="one point per composition group (group-weighted)"),
         Rectangle((0, 0), 1, 1, facecolor=PALETTE["green_1"],
                   label=r"exchangeable range $[1-\alpha,\ 1-\alpha+1/(n+1)]$"),
         Line2D([], [], marker="o", color=PALETTE["blue_main"], lw=1.4, ms=5, label="band gap: median, IQR"),
         Line2D([], [], marker="s", color=PALETTE["teal"], lw=1.4, ms=5, label="formation energy: median, IQR")]
    fig.legend(handles=h, loc="outside lower center", ncol=3, fontsize=6.8, columnspacing=1.0, handletextpad=0.4)
    return finalize_figure(fig, os.path.join(out, "fig2_group_level_mondrian")), None


# --------------------------------------------------------------------------------------- Fig. 3
P2_COLS = [  # (source, method, label)
    ("loso", "cvplus_random_in", "CV+ random\nfolds"),
    ("loso", "cvplus_lattice", "CV+ lattice\nfolds"),
    ("loso", "ncvplus_knn", "Normalized\n(kNN)"),
    ("loso", "cqr_cvplus", "CQR"),
    ("loso", "group_level", "Group-level\n(K = 5)"),
    ("loso", "cvplus_lattice_noFd3m_cal", "Calibration\nwithout " + LATTICE_TEX["Fd-3m"]),
    ("l2o", "cvplus_lattice", "Refit without\n" + LATTICE_TEX["Fd-3m"]),
]


def _lattice_sd(inp: Inputs, target: str) -> dict:
    sd = {}
    for l in LATTICE_ORDER:
        z = _loso_npz(inp, target, l)
        if z is not None:
            sd[l] = float(np.std(z["y_test"], ddof=1))
    return sd


def fig3_unseen_lattice(inp: Inputs, out: str):
    miss = _missing(inp, [("phase2", "report", "phase2_loso_all.csv"), ("phase2", "report", "phase2_l2o_all.csv")])
    if miss:
        return None, f"missing {miss}"
    A = inp.csv("phase2", "report", "phase2_loso_all.csv")
    B = inp.csv("phase2", "report", "phase2_l2o_all.csv")
    cmap, norm = coverage_cmap(), coverage_norm(NOMINAL)
    fig = new_figure(DOUBLE_COL_IN, 6.0, hspace=0.03)
    axes = fig.subplots(2, 1)
    for r, (target, (tname, unit)) in enumerate(TARGETS.items()):
        sd = _lattice_sd(inp, target)
        V = np.full((len(LATTICE_ORDER), len(P2_COLS)), np.nan)
        txt = [[""] * len(P2_COLS) for _ in LATTICE_ORDER]
        for i, l in enumerate(LATTICE_ORDER):
            for j, (src, m, _) in enumerate(P2_COLS):
                D = A if src == "loso" else B
                row = D[(D.target == target) & (D.heldout == l) & (D.method == m) & (D.group == "heldout")
                        & np.isclose(D.alpha, ALPHA)]
                if row.empty:
                    continue
                cov, w = float(row.coverage.iloc[0]), float(row.width_median.iloc[0])
                V[i, j] = cov
                if not np.isfinite(w):
                    wt = "width ∞"
                elif l in sd:
                    wt = f"{w / sd[l]:.1f} SD"
                else:
                    wt = f"{w:.2f} {unit}"
                txt[i][j] = f"{cov:.2f}\n{wt}"
        ax = axes[r]
        _heatmap(ax, V, _lat_labels(LATTICE_ORDER), [c[2] for c in P2_COLS], cmap, norm, txt, fontsize=6.8)
        ax.xaxis.tick_top()
        ax.tick_params(axis="x", labelsize=7)
        ax.set_ylabel("Held-out lattice")
        ax.axvline(1, color=PALETTE["ink"], lw=1.2)
        ax.axvline(5, color=PALETTE["ink"], lw=1.2)
        if r == 1:
            ax.set_xticklabels([])
        panel_label(ax, "ab"[r], tname)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=axes, location="bottom", fraction=0.04, aspect=50, pad=0.02)
    cb.set_ticks([0, 0.25, 0.5, 0.75, NOMINAL, 1.0])
    cb.set_label("Coverage of the held-out lattice at nominal 0.90   (cell: coverage, median width in lattice SD)",
                 fontsize=7.5)
    cb.outline.set_visible(False)
    return finalize_figure(fig, os.path.join(out, "fig3_unseen_lattice_coverage")), None


# --------------------------------------------------------------------------------------- Fig. 4
def _transfer_matrix(inp: Inputs, target: str):
    """MAE on the held-out lattice of each lattice-fold model (column = training lattice it left out),
    plus the leave-two-out ensemble; and the five-model LOSO ensemble MAE per row."""
    cols = LATTICE_ORDER + ["L2O"]
    M = np.full((len(LATTICE_ORDER), len(cols)), np.nan)
    ens = {}
    for i, l in enumerate(LATTICE_ORDER):
        z = _loso_npz(inp, target, l)
        if z is None:
            continue
        mu, y, lats = z["mu_te_l"], z["y_test"], [SG2LAT[int(s)] for s in z["lattices"]]
        for k, lo in enumerate(lats):
            M[i, cols.index(lo)] = np.mean(np.abs(mu[k] - y))
        ens[l] = float(np.mean(np.abs(mu.mean(0) - y)))
        z2 = _loso_npz(inp, target, l, "l2o") if l != "Fd-3m" else None
        if z2 is not None:
            M[i, -1] = float(np.mean(np.abs(z2["mu_te"].mean(0) - z2["y_test"])))
    return M, cols, ens


def _transfer_panel(ax, M, cols, ens, cax=None, fmt="{:.2f}"):
    R = np.log2(M / np.array([ens.get(l, np.nan) for l in LATTICE_ORDER])[:, None])
    cmap = LinearSegmentedColormap.from_list("transfer", [PALETTE["blue_main"], "#B9CDE5", "#FFFFFF",
                                                          PALETTE["red_2"], PALETTE["red_strong"]])
    norm = TwoSlopeNorm(vmin=-1.6, vcenter=0.0, vmax=1.6)
    txt = [[("" if np.isnan(v) else fmt.format(v)) for v in row] for row in M]
    labels = [LATTICE_TEX[c] for c in LATTICE_ORDER] + ["refit w/o\nboth"]
    _heatmap(ax, R, _lat_labels(LATTICE_ORDER), labels, cmap, norm, txt, na_text="", fontsize=6.4)
    plt.setp(ax.get_xticklabels(), rotation=40, ha="right", rotation_mode="anchor")
    j = cols.index("Fd-3m")
    ax.add_patch(Rectangle((j, 0), 1, len(LATTICE_ORDER) - 1, fill=False, edgecolor=PALETTE["ink"], lw=1.4,
                           zorder=4))
    ax.tick_params(axis="x", labelsize=6.8)
    ax.set_xlabel("Training lattice left out")
    if cax is not None:
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = plt.colorbar(sm, cax=cax, orientation="vertical")
        cb.set_ticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5])
        cb.set_ticklabels(["0.35×", "0.5×", "0.7×", "1×", "1.4×", "2×", "2.8×"])
        cb.set_label("MAE relative to the LOSO ensemble", fontsize=6.8)
        cb.ax.tick_params(labelsize=6.3)
        cb.outline.set_visible(False)


def fig4_mechanism(inp: Inputs, out: str):
    need = [("phase2", "report", "phase2_loso_all.csv"), ("phase2", "report", "phase2_l2o_all.csv")]
    need += [("phase2", "loso", f"{t}_loso{LAT2SG[l]}.npz") for t in TARGETS for l in LATTICE_ORDER]
    miss = _missing(inp, need)
    if miss:
        return None, f"missing {len(miss)} input file(s), e.g. {miss[0]}"
    A = inp.csv("phase2", "report", "phase2_loso_all.csv")
    B = inp.csv("phase2", "report", "phase2_l2o_all.csv")
    fig = new_figure(DOUBLE_COL_IN, 6.6, wspace=0.06, hspace=0.08)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.12], width_ratios=[1.12, 1])

    for c, (target, (tname, unit)) in enumerate(TARGETS.items()):
        ax = fig.add_subplot(gs[0, c])
        lo, hi = 0.0, 0.0
        stats = []
        for i, l in enumerate(LATTICE_ORDER):
            z = _loso_npz(inp, target, l)
            e = z["mu_te_l"].mean(0) - z["y_test"]
            q = np.percentile(e, [0.5, 99.5])
            lo, hi = min(lo, q[0]), max(hi, q[1])
            vp = ax.violinplot(e, positions=[i], vert=False, widths=0.82, showextrema=False, points=200)
            for b in vp["bodies"]:
                b.set_facecolor(LATTICE_COLORS[l])
                b.set_edgecolor(LATTICE_COLORS[l])
                b.set_alpha(0.42)
                b.set_linewidth(0.8)
            q1, med, q3 = np.percentile(e, [25, 50, 75])
            ax.plot([q1, q3], [i, i], color=PALETTE["ink"], lw=2.6, solid_capstyle="butt", zorder=3)
            ax.scatter(med, i, s=10, color="white", zorder=4)
            ax.scatter(e.mean(), i, s=40, marker="D", color=LATTICE_COLORS[l], edgecolor="white", lw=0.8, zorder=5)
            stats.append(abs(e.mean()) / np.abs(e).mean())
        ax.axvline(0, color=PALETTE["ink"], lw=1.0)
        span = hi - lo
        ax.set_xlim(lo - 0.05 * span, hi + 0.24 * span)
        for i, v in enumerate(stats):
            ax.text(0.995, i, f"{v:.2f}", transform=ax.get_yaxis_transform(), va="center", ha="right",
                    fontsize=6.8, color=PALETTE["grey_dark"])
        ax.text(0.995, -0.62, "|bias|/MAE", transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=6.3, color=PALETTE["grey_dark"])
        ax.set_yticks(range(len(LATTICE_ORDER)))
        ax.set_yticklabels(_lat_labels(LATTICE_ORDER))
        ax.set_ylim(len(LATTICE_ORDER) - 0.4, -0.85)
        ax.set_xlabel(f"Prediction − DFT ({unit})")
        panel_label(ax, "ab"[c], f"{tname}: error on the held-out lattice")

    sub = gs[1, 0].subgridspec(1, 2, width_ratios=[1, 0.05])
    ax = fig.add_subplot(sub[0, 0])
    cax = fig.add_subplot(sub[0, 1])
    M, cols, ens = _transfer_matrix(inp, "bandgap")
    _transfer_panel(ax, M, cols, ens, cax=cax)
    ax.set_ylabel("Held-out lattice")
    panel_label(ax, "c", "Band gap: MAE (eV) of each fold model")

    ax = fig.add_subplot(gs[1, 1])
    xs = [0, 1, 2]
    for target, ls, mk in (("bandgap", "-", "o"), ("formation", (0, (3, 2)), "s")):
        for l in [x for x in LATTICE_ORDER if x != "Fd-3m"]:
            def cov(D, m):
                r = D[(D.target == target) & (D.heldout == l) & (D.method == m) & (D.group == "heldout")
                      & np.isclose(D.alpha, ALPHA)]
                return float(r.coverage.iloc[0]) if len(r) else np.nan
            ys = [cov(A, "cvplus_lattice"), cov(A, "cvplus_lattice_noFd3m_cal"), cov(B, "cvplus_lattice")]
            emph = l == "P63/mmc"
            ax.plot(xs, ys, ls=ls, marker=mk, ms=5 if emph else 4, lw=2.2 if emph else 1.3,
                    color=LATTICE_COLORS[l], alpha=0.95 if emph else 0.75,
                    mfc=LATTICE_COLORS[l] if target == "bandgap" else "white", mew=1.2, zorder=4 if emph else 3)
    ax.axhline(NOMINAL, color=PALETTE["ink"], lw=1.0, ls=(0, (4, 3)))
    ax.set_xticks(xs)
    ax.set_xticklabels(["LOSO", "calibration\nwithout " + LATTICE_TEX["Fd-3m"],
                        "refit\nwithout " + LATTICE_TEX["Fd-3m"]])
    ax.set_xlim(-0.25, 2.25)
    ax.set_ylim(0.2, 1.03)
    ax.set_ylabel("Coverage of the held-out lattice")
    ax.annotate(LATTICE_TEX["P63/mmc"] + " loses its\nhigh-error relative", xy=(1.95, 0.335), xytext=(0.05, 0.62),
                fontsize=6.6, color=LATTICE_COLORS["P63/mmc"], ha="left",
                arrowprops=dict(arrowstyle="-|>", color=LATTICE_COLORS["P63/mmc"], lw=1.0, mutation_scale=7))
    hl = [Line2D([], [], color=LATTICE_COLORS[l], lw=1.8, label=LATTICE_TEX[l])
          for l in LATTICE_ORDER if l != "Fd-3m"]
    hl += [Line2D([], [], color=PALETTE["grey_dark"], marker="o", lw=1.3, ms=4, label="band gap"),
           Line2D([], [], color=PALETTE["grey_dark"], marker="s", mfc="white", lw=1.3, ls=(0, (3, 2)), ms=4,
                  label="formation")]
    ax.legend(handles=hl, loc="lower left", ncol=2, fontsize=6.4, handlelength=1.6, columnspacing=0.8)
    panel_label(ax, "d", "Who supplies the width?")
    return finalize_figure(fig, os.path.join(out, "fig4_mechanism_offset_transfer")), None


# --------------------------------------------------------------------------------------- Fig. 5
def fig5_overlap(inp: Inputs, out: str):
    miss = _missing(inp, [("phase2", "overlap", "overlap_weights.npz"), ("phase2", "overlap", "overlap_diagnostics.csv"),
                          ("phase2", "inputs", "train.csv")])
    if miss:
        return None, f"missing {miss}"
    W = inp.npz("phase2", "overlap", "overlap_weights.npz")
    O = inp.csv("phase2", "overlap", "overlap_diagnostics.csv")
    tr = inp.csv("phase2", "inputs", "train.csv").set_index("id")
    sg = tr.loc[W["id"], "spacegroup"].to_numpy()
    from scipy.stats import gaussian_kde

    fig = new_figure(DOUBLE_COL_IN, 3.4, wspace=0.06)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1])
    ax = fig.add_subplot(gs[0, 0])
    grid = np.linspace(-14, 14, 400)
    nL = len(LATTICE_ORDER)
    for i, l in enumerate(LATTICE_ORDER):
        p = np.clip(W[f"p_logistic_{LAT2SG[l]}"], 1e-6, 1 - 1e-6)
        lg = np.log(p / (1 - p))
        base = nL - 1 - i
        for mask, col, a in ((sg != LAT2SG[l], PALETTE["neutral"], 0.85), (sg == LAT2SG[l], LATTICE_COLORS[l], 0.7)):
            dens = gaussian_kde(lg[mask], bw_method=0.25)(grid)
            dens = 0.85 * dens / dens.max()
            ax.fill_between(grid, base, base + dens, color=col, alpha=a, lw=0, zorder=2 + i)
            ax.plot(grid, base + dens, color="white", lw=0.8, zorder=2 + i)
        row = O[O.heldout == l]
        if len(row):
            ax.text(13.8, base + 0.55, f"AUROC {row.auroc_logistic.iloc[0]:.3f}", fontsize=6.3, va="bottom",
                    ha="right", color=PALETTE["grey_dark"], zorder=30)
        ax.axhline(base, color=PALETTE["grey_mid"], lw=0.5, zorder=1)
    ax.axvline(0, color=PALETTE["ink"], lw=0.9, ls=(0, (4, 3)), zorder=20)
    ax.set_yticks([nL - 1 - i + 0.3 for i in range(nL)])
    ax.set_yticklabels(_lat_labels(LATTICE_ORDER))
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-14, 14)
    ax.set_ylim(-0.1, nL)
    ax.set_xlabel(r"Out-of-fold log-odds $\log[p/(1-p)]$ of 'held-out lattice'")
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=PALETTE["neutral"], label="other five lattices"),
                       Rectangle((0, 0), 1, 1, color=PALETTE["grey_dark"], alpha=0.7,
                                 label="held-out lattice (coloured)")],
              loc="upper left", ncol=1, fontsize=6.5, bbox_to_anchor=(0.0, 1.02))
    panel_label(ax, "a", "Domain classifier, SOAP PCA-50")

    ax = fig.add_subplot(gs[0, 1])
    for i, l in enumerate(LATTICE_ORDER):
        row = O[O.heldout == l]
        if row.empty:
            continue
        share, ess = row.median_test_weight_share.iloc[0], row.ess_frac.iloc[0]
        ax.plot([0, share], [i, i], color=PALETTE["neutral"], lw=2.4, zorder=1, solid_capstyle="round")
        ax.scatter(share, i, s=60, color=LATTICE_COLORS[l], edgecolor="white", lw=0.8, zorder=3)
        ax.scatter(ess, i, s=36, marker="D", facecolor="white", edgecolor=PALETTE["grey_dark"], lw=1.2, zorder=3)
    ax.axvline(ALPHA, color=PALETTE["red_strong"], lw=1.2, ls=(0, (4, 3)))
    ax.text(ALPHA + 0.015, nL - 0.5, r"$\alpha$ = 0.1", color=PALETTE["red_strong"], fontsize=6.8, va="bottom")
    ax.set_yticks(range(nL))
    ax.set_yticklabels(_lat_labels(LATTICE_ORDER))
    ax.set_ylim(nL - 0.4, -0.6)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Fraction")
    ax.legend(handles=[Line2D([], [], marker="o", color=PALETTE["grey_dark"], lw=0, ms=6,
                              label="median weight on the test point"),
                       Line2D([], [], marker="D", mfc="white", mec=PALETTE["grey_dark"], lw=0, ms=5,
                              label=r"Kish ESS / $n$ of the calibration weights")],
              loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=1, fontsize=6.6)
    panel_label(ax, "b", r"Weighted CP: unbounded when test weight > $\alpha$")
    return finalize_figure(fig, os.path.join(out, "fig5_no_covariate_overlap")), None


# --------------------------------------------------------------------------------------- Fig. 6
def fig6_detection(inp: Inputs, out: str):
    rep = ("phase4", "report")
    names = ["phase4_detection_metrics.csv", "phase4_abstention.csv", "phase4_lattice_level.csv",
             "phase4_lattice_correlations.csv", "phase4_flags.csv"]
    miss = _missing(inp, [rep + (n,) for n in names])
    if miss:
        return None, f"missing {miss} (run phase4 and report4 first)"
    Mt = inp.csv(*rep, names[0])
    Ab = inp.csv(*rep, names[1])
    Lv = inp.csv(*rep, names[2])
    Co = inp.csv(*rep, names[3])
    Fl = inp.csv(*rep, names[4])
    fig = new_figure(DOUBLE_COL_IN, 6.0, wspace=0.05, hspace=0.07)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.2, 1, 1])

    for r, (target, (tname, _)) in enumerate(TARGETS.items()):
        ax = fig.add_subplot(gs[r, 0])
        M = Mt[(Mt.target == target) & (Mt.outcome == "miscovered")]
        for i, (s, lab, col) in enumerate(DETECT_SCORES):
            if s in ("knn", "kde"):
                ax.axhspan(i - 0.48, i + 0.48, color=PALETTE["highlight"], alpha=0.13, lw=0, zorder=0)
            for k, (lev, mk, fc) in enumerate((("pooled", "o", col), ("within (n-weighted)", "s", "white"))):
                row = M[(M.score == s) & (M.level == lev)]
                if row.empty or not np.isfinite(row.auroc.iloc[0]):
                    continue
                y = i + (k - 0.5) * 0.34
                a, lo, hi = row.auroc.iloc[0], row.auroc_lo.iloc[0], row.auroc_hi.iloc[0]
                if np.isfinite(lo) and np.isfinite(hi):
                    ax.plot([lo, hi], [y, y], color=col, lw=1.6, solid_capstyle="round", zorder=2)
                ax.scatter(a, y, s=30, marker=mk, facecolor=fc, edgecolor=col, lw=1.4, zorder=3)
        ax.axvline(0.5, color=PALETTE["ink"], lw=1.0)
        ax.axvline(0.7, color=PALETTE["grey_mid"], lw=0.9, ls=(0, (2, 2)))
        ax.set_yticks(range(len(DETECT_SCORES)))
        ax.set_yticklabels([s[1] for s in DETECT_SCORES], fontsize=6.8)
        ax.set_ylim(len(DETECT_SCORES) - 0.5, -0.5)
        ax.set_xlim(0, 1)
        ax.set_xlabel("AUROC for a missed interval")
        panel_label(ax, "ad"[r], f"{tname}: detection")

    for r, (target, (tname, _)) in enumerate(TARGETS.items()):
        ax = fig.add_subplot(gs[r, 1])
        A = Ab[Ab.target == target]
        for s, lab, col in DETECT_SCORES:
            S = A[A.score == s].sort_values("abstained")
            if S.empty:
                continue
            ax.plot(S.abstained, S.coverage_retained, color=col, lw=1.8, marker="o", ms=2.6, alpha=0.9,
                    zorder=3 if s in ("knn", "kde") else 2)
        S = A[A.score == "oracle_abs_err"].sort_values("abstained")
        if len(S):
            ax.plot(S.abstained, S.coverage_retained, color=PALETTE["ink"], lw=1.2, ls=(0, (3, 2)), zorder=1)
        ax.axhline(NOMINAL, color=PALETTE["ink"], lw=0.9, ls=(0, (4, 3)))
        ax.set_xlabel("Fraction abstained")
        ax.set_ylabel("Coverage of the retained")
        ax.set_xlim(0, 0.9)
        ax.set_ylim(max(0.0, np.floor(A.coverage_retained.min() * 20) / 20 - 0.02), 1.01)
        panel_label(ax, "be"[r], "Selective prediction")

    ax = fig.add_subplot(gs[0, 2])
    for target, fill in (("bandgap", True), ("formation", False)):
        for _, row in Lv[Lv.target == target].iterrows():
            col = LATTICE_COLORS[row.heldout]
            ax.scatter(row.kde_mean, row.miscoverage, s=34, marker="o" if row.design == "loso" else "^",
                       facecolor=col if fill else "white", edgecolor=col, lw=1.3, zorder=3)
    lines = []
    for target, (tname, _) in TARGETS.items():
        c = Co[(Co.target == target) & (Co.configs == "all 11") & (Co.score == "kde_mean")
               & (Co.outcome == "miscoverage")]
        if len(c):
            lines.append(f"{'BG' if target == 'bandgap' else 'FE'} ρ = {c.spearman.iloc[0]:.2f}, "
                         f"p = {c.p_perm.iloc[0]:.2g}")
    ax.axhline(ALPHA, color=PALETTE["ink"], lw=0.9, ls=(0, (4, 3)))
    ax.set_xlabel("Mean KDE dissimilarity")
    ax.set_ylabel("Miscoverage of the lattice")
    panel_label(ax, "c", "Lattice level (n = 11)\n" + "; ".join(lines))

    ax = fig.add_subplot(gs[1, 2])
    F = Fl.copy()
    F["order"] = F.heldout.map({l: i for i, l in enumerate(LATTICE_ORDER)}) + (F.design != "loso") * 0.5
    F = F.sort_values("order").reset_index(drop=True)
    for i, row in F.iterrows():
        col = LATTICE_COLORS[row.heldout]
        ax.plot([row.ffr_ci_lo, row.ffr_ci_hi], [i, i], color=PALETTE["grey_mid"], lw=1.4, zorder=1)
        ax.scatter(row.false_flag_rate, i, s=16, color=PALETTE["grey_mid"], zorder=2, clip_on=False)
        ax.scatter(row.flag_rate_heldout, i, s=40, color=col, marker="o" if row.design == "loso" else "^",
                   edgecolor="white", lw=0.6, zorder=3, clip_on=False)
    ax.axvline(0.05, color=PALETTE["red_strong"], lw=1.0, ls=(0, (4, 3)))
    ax.set_yticks(range(len(F)))
    ax.set_yticklabels([LATTICE_TEX[h] for h in F.heldout], fontsize=6.5)
    ax.set_ylim(len(F) - 0.5, -0.5)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel(r"Flag rate (conformal $p \leq 0.05$)")
    panel_label(ax, "f", "Outlier flag (kNN score)")

    h = [Line2D([], [], color=c, lw=1.8, label=l) for _, l, c in DETECT_SCORES]
    h += [Line2D([], [], color=PALETTE["ink"], lw=1.2, ls=(0, (3, 2)), label="oracle (true |error|)"),
          Line2D([], [], marker="o", color=PALETTE["grey_dark"], lw=1.4, ms=4.5, label="AUROC pooled, 95 % CI"),
          Line2D([], [], marker="s", mfc="white", color=PALETTE["grey_dark"], lw=1.4, ms=4.5,
                 label="AUROC within lattice"),
          Line2D([], [], marker="o", color=PALETTE["grey_dark"], lw=0, ms=5, label="band gap / LOSO"),
          Line2D([], [], marker="o", mfc="white", color=PALETTE["grey_dark"], lw=0, ms=5, label="formation"),
          Line2D([], [], marker="^", color=PALETTE["grey_dark"], lw=0, ms=5,
                 label="refit without " + LATTICE_TEX["Fd-3m"]),
          Line2D([], [], marker="o", color=PALETTE["grey_mid"], lw=1.4, ms=3.5,
                 label="in-distribution false-flag rate, binomial 95 % CI")]
    fig.legend(handles=h, loc="outside lower center", ncol=4, fontsize=6.5, columnspacing=1.0, handletextpad=0.4)
    return finalize_figure(fig, os.path.join(out, "fig6_failure_detection")), None


# --------------------------------------------------------------------------------------- supplementary
def figS1_calibration_tail(inp: Inputs, out: str):
    miss = _missing(inp, [("phase2", "report", "phase2_calibration_composition.csv")])
    if miss:
        return None, f"missing {miss}"
    C = inp.csv("phase2", "report", "phase2_calibration_composition.csv")
    C = C[np.isclose(C.alpha, ALPHA) & (C.folds == "lattice")]
    fig = new_figure(DOUBLE_COL_IN, 3.1, wspace=0.04)
    axes = fig.subplots(1, 2)
    norm = LogNorm(vmin=0.05, vmax=5)
    cmap = sequential_red()
    for c, (target, (tname, _)) in enumerate(TARGETS.items()):
        V = np.full((len(LATTICE_ORDER), len(LATTICE_ORDER)), np.nan)
        txt = [[""] * len(LATTICE_ORDER) for _ in LATTICE_ORDER]
        for i, h in enumerate(LATTICE_ORDER):
            for j, cl in enumerate(LATTICE_ORDER):
                r = C[(C.target == target) & (C.heldout == h) & (C.cal_lattice == cl)]
                if len(r):
                    V[i, j] = max(float(r.ratio.iloc[0]), 0.05)
                    txt[i][j] = f"{float(r.ratio.iloc[0]):.2f}"
        _heatmap(axes[c], V, _lat_labels(LATTICE_ORDER), _lat_labels(LATTICE_ORDER), cmap, norm, txt,
                 na_text="", fontsize=6.6)
        axes[c].set_xlabel("Lattice of the calibration residual")
        if c == 0:
            axes[c].set_ylabel("Held-out lattice")
        else:
            axes[c].set_yticklabels([])
        panel_label(axes[c], "ab"[c], tname)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=axes, fraction=0.04, aspect=30, pad=0.02)
    cb.set_label("Share of the top 10 % / share of calibration points", fontsize=6.8)
    cb.set_ticks([0.1, 0.5, 1, 2, 5])
    cb.set_ticklabels(["0.1", "0.5", "1", "2", "5"])
    cb.outline.set_visible(False)
    return finalize_figure(fig, os.path.join(out, "figS1_calibration_tail_composition")), None


def figS2_transfer_formation(inp: Inputs, out: str):
    need = [("phase2", "loso", f"formation_loso{LAT2SG[l]}.npz") for l in LATTICE_ORDER]
    miss = _missing(inp, need)
    if miss:
        return None, f"missing {len(miss)} file(s)"
    fig = new_figure(4.6, 3.1)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.05])
    ax, cax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    M, cols, ens = _transfer_matrix(inp, "formation")
    _transfer_panel(ax, M, cols, ens, cax=cax, fmt="{:.3f}")
    ax.set_ylabel("Held-out lattice")
    panel_label(ax, "", "Formation energy: MAE (eV/atom) of each fold model")
    return finalize_figure(fig, os.path.join(out, "figS2_transfer_formation")), None


def figS3_detection_radar(inp: Inputs, out: str):
    miss = _missing(inp, [("phase4", "report", "phase4_detection_metrics.csv")])
    if miss:
        return None, f"missing {miss}"
    Mt = inp.csv("phase4", "report", "phase4_detection_metrics.csv")
    settings = [(t, o, lv) for t in TARGETS for o in ("miscovered", "large_error")
                for lv in ("pooled", "within (n-weighted)")]
    slabels = [f"{'Band gap' if t == 'bandgap' else 'Formation'}\n{'missed' if o == 'miscovered' else 'large error'}"
               f", {'pooled' if lv == 'pooled' else 'within'}" for t, o, lv in settings]
    ang = np.linspace(0, 2 * np.pi, len(settings), endpoint=False)
    fig = new_figure(4.8, 4.9)
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    for s, lab, col in DETECT_SCORES:
        vals = []
        for t, o, lv in settings:
            r = Mt[(Mt.target == t) & (Mt.outcome == o) & (Mt.level == lv) & (Mt.score == s)]
            vals.append(float(r.auroc.iloc[0]) if len(r) and np.isfinite(r.auroc.iloc[0]) else np.nan)
        v = np.array(vals + vals[:1])
        a = np.r_[ang, ang[:1]]
        key = s in ("knn", "kde")
        ax.plot(a, v, color=col, lw=2.0 if key else 1.4, alpha=0.9, label=lab)
        ax.fill(a, np.nan_to_num(v, nan=0.5), color=col, alpha=0.12 if key else 0.05)
        ax.scatter(ang, vals, s=12, color=col, zorder=3)
    ax.plot(np.linspace(0, 2 * np.pi, 200), np.full(200, 0.5), color=PALETTE["ink"], lw=1.0, ls=(0, (3, 2)))
    ax.set_xticks(ang)
    ax.set_xticklabels(slabels, fontsize=6.3)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0.25", "0.5", "0.75", "1"], fontsize=6)
    ax.set_rlabel_position(22)
    ax.spines["polar"].set_color(PALETTE["neutral"])
    ax.grid(color=PALETTE["neutral"], lw=0.6)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncol=2, fontsize=6.6)
    return finalize_figure(fig, os.path.join(out, "figS3_detection_radar")), None


FIGURES = [fig1_indistribution, fig2_group_mondrian, fig3_unseen_lattice, fig4_mechanism, fig5_overlap,
           fig6_detection, figS1_calibration_tail, figS2_transfer_formation, figS3_detection_radar]


def make_figures(phase1: str, phase2: str, phase4: str, out: str, only: list[str] | None = None) -> dict:
    apply_publication_style()
    os.makedirs(out, exist_ok=True)
    inp = Inputs(phase1, phase2, phase4)
    done, skipped = {}, {}
    for f in FIGURES:
        if only and not any(f.__name__.startswith(o) for o in only):
            continue
        try:
            saved, why = f(inp, out)
        except Exception as e:  # report and carry on: one bad input must not block the rest
            plt.close("all")
            saved, why = None, f"error: {type(e).__name__}: {e}"
        if saved:
            done[f.__name__] = [os.path.basename(str(p)) for p in saved]
            print(f"[figures] {f.__name__}: {', '.join(done[f.__name__])}", flush=True)
        else:
            skipped[f.__name__] = why
            print(f"[figures] {f.__name__}: skipped ({why})", flush=True)
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    style="figures4papers house style (Chen Liu, CC BY-NC 4.0), adapted in cpaudit/pubstyle.py",
                    figures=done, skipped=skipped, inputs_sha256=inp.used,
                    phase4_note="Fig. 6 and S3 are drawn from phase04 outputs; check phase4_decisions.json first.")
    with open(os.path.join(out, "FIGURES_MANIFEST.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    return manifest


def run_figures(cfg: dict, only: list[str] | None = None) -> dict:
    p = cfg["paths"]
    root = p["drive_root"]
    if not os.path.isdir(root):
        raise FileNotFoundError(f"Drive folder not found: {root}. Mount Drive first.")
    out = os.path.join(root, p.get("figures_dir", "paper_figures_v1"))
    return make_figures(os.path.join(root, p.get("phase1_dir", "phase01_v1")), os.path.join(root, p["phase2_dir"]),
                        os.path.join(root, p["out_dir"]), out, only)
