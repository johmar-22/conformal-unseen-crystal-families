"""Paper figures for the JARVIS replication, v8.1 (presentation fixes found when the v8 figures were
inspected; the numbers come from the same files). Written to paper_figures_v4/ (v3 is left as it is).

Fig. 8   as v8, but every worst-group label is joined to its point by a leader line (in v8 the P4mm
         label sat next to the Pm-3m point).
Fig. 9   as v8, but 0.70 in panel c is labelled as the G4-J threshold (not "nominal"), 0.50 as chance,
         and 0.80 in panels a, b as the bound of the hierarchical jackknife+ (Theorem 5).
Fig. 10  redesigned. On JARVIS the mean few-shot coverage equals k/(m+1) by construction, so the
         v8 figure showed flat lines only. Now: coverage over draws (mean and 5-95 % band) for the
         three worst groups, and the width of the model intervals relative to a no-model interval
         built from the same labels (posthoc8; POST HOC, not pre-registered).
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .an8 import sg_tex
from .figs8 import TNAME, fig8_jarvis_groups, fig9_jarvis_guarantees_detection
from .pubstyle import DOUBLE_COL_IN, OVERLAP_LOG, PALETTE, apply_publication_style, finalize_figure, new_figure, panel_label

GROUP_COLORS = (PALETTE["blue_main"], PALETTE["teal"], PALETTE["violet"])
GROUP_MARKERS = ("o", "s", "D")


def nominal_split(m: int, alpha: float) -> float:
    k = math.ceil((m + 1) * (1 - alpha) - 1e-9)
    return min(k, m + 1) / (m + 1)


def fig10_jarvis_fewshot_v81(d8, posthoc, out, alpha=0.10):
    f = os.path.join(posthoc, "fewshot_by_m.csv")
    if not os.path.exists(f):
        return None, "run posthoc8 first"
    M = pd.read_csv(f)
    M = M[M.scheme == "one_per_group"]
    W = pd.read_csv(os.path.join(d8, "fewshot", "worst_groups.csv"))
    targets = [t for t in ("formation", "bandgap") if (W.target == t).any()]
    ncol = int(W.groupby("target").size().max())
    fig = new_figure(DOUBLE_COL_IN, 2.25 * len(targets) + 0.75, wspace=0.05, hspace=0.12)
    gs = fig.add_gridspec(len(targets), ncol + 1, width_ratios=[1] * ncol + [1.2])
    letters = iter("abcdefghij")
    mg = sorted(M.m.unique())
    for i, t in enumerate(targets):
        Wt = W[W.target == t].reset_index(drop=True)
        axr = None
        for j, w in Wt.iterrows():
            ax = fig.add_subplot(gs[i, j])
            d = M[(M.target == t) & (M.group == w.group)].sort_values("m")
            ax.fill_between(d.m, d.cov_q05_split, d.cov_q95_split, color=PALETTE["blue_main"], alpha=0.16, lw=0)
            ax.plot(d.m, d.cov_mean_split, color=PALETTE["blue_main"], marker="o", ms=3, lw=1.3)
            ax.plot(d.m, d.cov_mean_offset_full, color=PALETTE["green_3"], marker="o", ms=2.6, lw=1.1, ls=(0, (3, 1.5)))
            ax.plot(mg, [nominal_split(m, alpha) for m in mg], color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)))
            ax.axhline(0.80, color=PALETTE["red_strong"], lw=0.6, ls=(0, (1, 2)))
            ax.axhline(float(w.primary_coverage), color=PALETTE["grey_dark"], lw=0.8, ls=(0, (1, 1)))
            ax.set_ylim(min(0.38, float(w.primary_coverage) - 0.06), 1.04)
            ax.set_yticks([v for v in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) if v >= ax.get_ylim()[0]])
            ax.set_xticks(mg)
            ax.set_xlim(min(mg) - 3, max(mg) + 3)
            ax.set_xlabel("Labelled structures m")
            if j == 0:
                ax.set_ylabel(f"Coverage, {TNAME[t].split(' (')[0].lower()}")
            panel_label(ax, next(letters), f"{sg_tex(w.symbol)} (unseen: {w.primary_coverage:.2f})")
        axr = fig.add_subplot(gs[i, ncol])
        top = 1.3
        for j, w in Wt.iterrows():
            d = M[(M.target == t) & (M.group == w.group)].sort_values("m")
            c, mk = GROUP_COLORS[j % 3], GROUP_MARKERS[j % 3]
            axr.plot(d.m, d.ratio_offset_to_nomodel, color=c, marker=mk, ms=3, lw=1.3, label=sg_tex(w.symbol))
            axr.plot(d.m, d.ratio_split_to_nomodel, color=c, marker=mk, ms=2.6, mfc="white", lw=0.9, ls=(0, (3, 1.5)))
            top = max(top, float(np.nanmax(d[["ratio_offset_to_nomodel", "ratio_split_to_nomodel"]].to_numpy())) + 0.12)
        axr.axhline(1.0, color=PALETTE["grey_dark"], lw=0.8)
        axr.set_ylim(0, top)
        axr.set_xticks(mg)
        axr.set_xlim(min(mg) - 3, max(mg) + 3)
        axr.set_xlabel("Labelled structures m")
        axr.set_ylabel("Width / no-model width")
        axr.legend(loc="lower left", fontsize=6.2, handlelength=1.6, borderaxespad=0.4, frameon=False)
        panel_label(axr, next(letters), "Interval width vs no model")
    h = [(Patch(facecolor=PALETTE["blue_main"], alpha=0.16, edgecolor="none"), "split: 5-95 % of draws"),
         (Line2D([], [], color=PALETTE["blue_main"], marker="o", ms=3), "split: mean coverage"),
         (Line2D([], [], color=PALETTE["green_3"], marker="o", ms=2.6, ls=(0, (3, 1.5))), "offset-corrected: mean coverage"),
         (Line2D([], [], color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3))), "exact mean k/(m+1)"),
         (Line2D([], [], color=PALETTE["red_strong"], lw=0.6, ls=(0, (1, 2))), "0.80"),
         (Line2D([], [], color=PALETTE["grey_dark"], lw=0.8, ls=(0, (1, 1))), "coverage without labels (unseen)"),
         (Line2D([], [], color=PALETTE["grey_dark"], lw=1.3), "ratio: offset-corrected / no model"),
         (Line2D([], [], color=PALETTE["grey_dark"], lw=0.9, ls=(0, (3, 1.5)), marker="o", ms=2.6, mfc="white"),
          "ratio: split / no model"),
         (Line2D([], [], color=PALETTE["grey_dark"], lw=0.8), "1 = as wide as a no-model interval")]
    fig.legend(handles=[a for a, _ in h], labels=[b for _, b in h], loc="outside lower center", ncol=3, fontsize=6.6)
    return finalize_figure(fig, os.path.join(out, "fig10_jarvis_fewshot")), None


def run_figures81(cfg):
    import matplotlib.pyplot as plt
    from .fit7 import dirs7
    from .posthoc8 import ph8
    from .run8 import dirs8
    from .fit7 import MOCK
    d7, d8 = dirs7(cfg)["root"], dirs8(cfg)
    P = ph8(cfg)
    root = os.path.dirname(d8["figures"])
    out = os.path.join(root, "paper_figures_mock81" if MOCK else P.get("figures_dir", "paper_figures_v4"))
    man = os.path.join(out, "FIGURES81_MANIFEST.json")
    if os.path.exists(man):
        raise SystemExit(f"[figures81] {out} already holds a figure set (FIGURES81_MANIFEST.json); nothing was "
                         "overwritten. Rename that folder first if the figures should be redrawn.")
    os.makedirs(out, exist_ok=True)
    posthoc = os.path.join(d8["root"], P["out_subdir"])
    apply_publication_style()
    OVERLAP_LOG.clear()
    a = float(cfg["models7"]["alpha_primary"])
    jobs = [("fig8_jarvis_groups", lambda: fig8_jarvis_groups(d7, d8["root"], out, a, leader_lines=True)),
            ("fig9_jarvis_guarantees_detection", lambda: fig9_jarvis_guarantees_detection(d7, d8["root"], out, a, v81=True)),
            ("fig10_jarvis_fewshot_v81", lambda: fig10_jarvis_fewshot_v81(d8["root"], posthoc, out, a))]
    done, skipped = {}, {}
    for name, fn in jobs:
        try:
            saved, why = fn()
        except Exception as e:
            plt.close("all")
            saved, why = None, f"error: {type(e).__name__}: {e}"
        if saved:
            stem = os.path.splitext(os.path.basename(str(saved[0])))[0]
            issues = OVERLAP_LOG.get(stem, [])
            done[name] = [os.path.basename(str(p)) for p in saved]
            print(f"[figures81] {name}: {', '.join(done[name])} | overlap check: "
                  f"{'clean' if not issues else f'{len(issues)} ISSUE(S)'}", flush=True)
            for d in issues:
                print(f"           {d['kind']}: {d['a']}  <->  {d['b']}", flush=True)
        else:
            skipped[name] = why
            print(f"[figures81] {name}: skipped ({why})", flush=True)
    if not skipped:        # the manifest marks a complete set; an incomplete run can simply be repeated
        with open(man, "w") as fh:
            json.dump(dict(figures=done, skipped=skipped, overlap_issues=OVERLAP_LOG), fh, indent=2, default=str)
    n = sum(len(v) for v in OVERLAP_LOG.values())
    print(f"[figures81] {len(done)} figures in {out}; overlap issues: {n}")
    return done
