"""Paper figures for the JARVIS replication (Phase 6, v8), drawn from phase07_v1 and phase08_v1 files only.

Fig. 8   unseen space groups: group coverage vs distance to training (a, b) and per-group coverage
         by method (c, d)
Fig. 9   guarantees and detection: H3 (HCP, single subsampling) and hierarchical jackknife+ (a);
         failure-detection AUROC, pooled and within group (b)
Fig. 10  few-shot recalibration of the three worst-covered groups per target (H5-J)
Every figure goes through the automatic overlap check (figcheck) when it is saved.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.text import Text

from .an8 import sg_tex
from .pubstyle import (DOUBLE_COL_IN, OVERLAP_LOG, PALETTE, apply_publication_style, coverage_cmap,
                       coverage_norm, finalize_figure, new_figure, panel_label)

TNAME = {"formation": "Formation energy", "bandgap": "Band gap (non-metals)", "tbmbj": "TBmBJ gap"}
METHODS = [("cvplus", "CV+", PALETTE["blue_main"]), ("cvplus_gw", "group-weighted CV+", PALETTE["teal"]),
           ("ncvplus", "normalized CV+", PALETTE["violet"]), ("cqr_cvplus", "CQR-CV+", PALETTE["ochre"]),
           ("randinner", "CV+, random inner folds", PALETTE["red_2"]), ("enssd", "ensemble SD", PALETTE["red_strong"])]
SCORE_LABEL = {"width": "CV+ width", "knn": "kNN distance", "kde": "KDE dissim.", "committee": "inner-model SD",
               "classifier": "domain classifier"}


def _place_labels(fig, ax, X, Y, S, items, color, hlines=()):
    """Annotate points without covering markers, horizontal reference lines, other labels or the axes
    edge: candidate offsets are tried in order on the drawn layout. As in the overlap check
    (figcheck), every marker counts with the radius of the largest one (sizes S in pt^2). A thin
    leader line joins labels placed away from their point."""
    prev = None
    for _ in range(6):                      # constrained layout can need a few passes to settle
        fig.draw_without_rendering()
        now = tuple(tuple(np.round(a_.get_position().bounds, 6)) for a_ in fig.axes)
        if now == prev:
            break
        prev = now
    r = fig.canvas.get_renderer()
    P = ax.transData.transform(np.c_[X, Y])
    rad = np.sqrt(np.max(np.asarray(S, float))) / 2 * r.points_to_pixels(1) + 2.0
    HL = [ax.transData.transform((X[0], h))[1] for h in hlines]
    axbb = ax.get_window_extent(r)
    cands = [(7, 0), (-7, 0), (8, -11), (8, 11), (-8, -11), (-8, 11), (0, -15), (0, 15)]
    for rr in (24, 36, 50, 66):             # farther spots (with a leader line), 16 directions each
        for k in range(16):
            ang = -np.pi / 2 + 2 * np.pi * k / 16
            cands.append((int(round(rr * np.cos(ang))), int(round(rr * np.sin(ang)))))
    placed = []

    def clear(bb):
        if not (axbb.x0 <= bb.x0 and bb.x1 <= axbb.x1 and axbb.y0 <= bb.y0 and bb.y1 <= axbb.y1):
            return False
        cx = np.clip(P[:, 0], bb.x0, bb.x1)
        cy = np.clip(P[:, 1], bb.y0, bb.y1)
        if np.any((P[:, 0] - cx) ** 2 + (P[:, 1] - cy) ** 2 < rad ** 2):
            return False
        if any(bb.y0 - 2 < h < bb.y1 + 2 for h in HL):
            return False
        return not any(bb.x0 < q.x1 and q.x0 < bb.x1 and bb.y0 < q.y1 and q.y0 < bb.y1 for q in placed)

    def make(x, y, t, dx, dy, arrow):
        ha = "left" if dx > 0 else ("right" if dx < 0 else "center")
        va = "center" if dy == 0 else ("bottom" if dy > 0 else "top")
        a = ax.annotate(t, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va=va, fontsize=7, color=color,
                        arrowprops=(dict(arrowstyle="-", lw=0.5, color=color, shrinkA=0, shrinkB=3) if arrow else None))
        a.set_in_layout(False)              # a label must not make the constrained layout move the axes
        return a

    for x, y, t in items:
        chosen, spot = None, cands[0]
        for dx, dy in cands:
            a = make(x, y, t, dx, dy, False)   # measured without a leader line: the box is the text only
            bb = a.get_window_extent(r).expanded(1.04, 1.08)
            a.remove()
            if clear(bb):
                chosen, spot = bb, (dx, dy)
                break
        a = make(x, y, t, spot[0], spot[1], abs(spot[0]) > 10 or abs(spot[1]) > 12)
        if chosen is None:                  # no clear spot: first candidate (the overlap check will report it)
            chosen = Text.get_window_extent(a, r)
        placed.append(chosen)


def _method_coverage(d7, target, alpha):
    rows = []
    for design, rep in (("primary", 1), ("rep2", 2), ("randinner", 1)):
        for o in range(10):
            f = os.path.join(d7, "unitsB", f"{target}_{design}_r{rep}_o{o:02d}.csv")
            if not os.path.exists(f):
                continue
            d = pd.read_csv(f)
            d = d[(d.group != "ALL") & np.isclose(d.alpha, alpha)]
            if design == "randinner":
                d = d[d.method == "cvplus"].assign(method="randinner")
            rows.append(d[["group", "method", "coverage"]].assign(rep=rep))
    if not rows:
        return None
    D = pd.concat(rows, ignore_index=True)
    return D.groupby(["method", "group"]).coverage.mean().reset_index()


def fig8_jarvis_groups(d7, d8, out, alpha=0.10):
    need = [os.path.join(d8, "report", f"{t}_groups.csv") for t in ("formation", "bandgap")]
    if not all(os.path.exists(f) for f in need):
        return None, "run report8 first"
    fig = new_figure(DOUBLE_COL_IN, 5.4, wspace=0.08, hspace=0.10)
    gs = fig.add_gridspec(2, 2)
    worst = json.load(open(os.path.join(d8, "phase0", "PREREGISTRATION_phase6_v8.json")))["worst_groups"]
    dec = json.load(open(os.path.join(d8, "report", "phase6_decisions_v8.json")))
    cmap, norm = coverage_cmap(), coverage_norm(1 - alpha)
    label_jobs = []
    for j, t in enumerate(("formation", "bandgap")):
        G = pd.read_csv(os.path.join(d8, "report", f"{t}_groups.csv"))
        ax = fig.add_subplot(gs[0, j])
        s = 8 + 60 * np.sqrt(G.n / G.n.max())
        ax.scatter(G.knn_median, G.coverage, s=s, c=G.coverage, cmap=cmap, norm=norm, edgecolor=PALETTE["ink"],
                   lw=0.4, zorder=3)
        ax.axhline(1 - alpha, color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)))
        ax.axhline(0.80, color=PALETTE["red_strong"], lw=0.7, ls=(0, (1, 2)))
        ax.set_ylim(min(0.45, G.coverage.min() - 0.08), 1.03)
        x0, x1 = G.knn_median.min(), G.knn_median.max()
        ax.set_xlim(x0 - 0.06 * (x1 - x0), x1 + 0.25 * (x1 - x0))
        ax.set_xlabel("Median distance to training (kNN)")
        ax.set_ylabel("Group coverage (CV+)")
        sp = dec[t]["across_groups"].get("spearman_coverage_knn_median")
        sub = (rf"Spearman $\rho$ = {sp['rho']:.2f} [{sp['ci95'][0]:.2f}, {sp['ci95'][1]:.2f}]" if sp else "")
        panel_label(ax, "ab"[j], f"{TNAME[t]}: {len(G)} held-out groups\n{sub}")
        wl = [w["group"] for w in worst[t]]
        W_ = G[G.group.isin(wl)].sort_values("coverage")
        label_jobs.append((ax, G.knn_median.to_numpy(), G.coverage.to_numpy(), s.to_numpy(),
                           [(r.knn_median, r.coverage, sg_tex(r.symbol)) for _, r in W_.iterrows()]))
        ax = fig.add_subplot(gs[1, j])
        M = _method_coverage(d7, t, alpha)
        meths = [m for m in METHODS if M is not None and m[0] in set(M.method)]
        rng = np.random.default_rng(0)
        for i, (mk, lab, col) in enumerate(meths):
            v = M[M.method == mk].coverage.to_numpy()
            ax.scatter(i + rng.uniform(-0.18, 0.18, v.size), v, s=7, color=col, alpha=0.75, lw=0, zorder=3)
            ax.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color=PALETTE["ink"], lw=1.4, zorder=4)
        ax.axhline(1 - alpha, color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)))
        ax.axhline(0.80, color=PALETTE["red_strong"], lw=0.7, ls=(0, (1, 2)))
        ax.set_xticks(range(len(meths)), [m[1] for m in meths], rotation=30, ha="right")
        ax.set_ylim(-0.02, 1.04)
        ax.set_ylabel("Group coverage")
        panel_label(ax, "cd"[j], f"{TNAME[t]}: per-group coverage by method")
    h = [Line2D([], [], color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)), label=f"nominal {1 - alpha:.2f}"),
         Line2D([], [], color=PALETTE["red_strong"], lw=0.7, ls=(0, (1, 2)), label="0.80"),
         Line2D([], [], color=PALETTE["ink"], lw=1.4, label="median over groups (c, d)"),
         Line2D([], [], marker="o", lw=0, color=PALETTE["grey_mid"], ms=4, label="marker area ~ group size (a, b)")]
    fig.legend(handles=h, loc="outside lower center", ncol=4, fontsize=7)
    for ax_, X_, Y_, S_, items in label_jobs:      # last, on the final layout
        _place_labels(fig, ax_, X_, Y_, S_, items, PALETTE["red_strong"], hlines=(1 - alpha, 0.80))
    return finalize_figure(fig, os.path.join(out, "fig8_jarvis_unseen_groups")), None


def fig9_jarvis_guarantees_detection(d7, d8, out, alpha=0.10):
    need = [os.path.join(d8, "report", f"{t}_h4_pooled.csv") for t in ("formation", "bandgap")]
    if not all(os.path.exists(f) for f in need):
        return None, "run report8 first"
    fig = new_figure(DOUBLE_COL_IN, 3.2, wspace=0.06)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.62, 0.62, 1.45])
    lows = []
    for j, t in enumerate(("formation", "bandgap")):
        ax = fig.add_subplot(gs[0, j])
        files = sorted(f for f in os.listdir(os.path.join(d7, "unitsH3")) if f.startswith(f"{t}_H3_") and f.endswith(".csv"))
        H = pd.concat([pd.read_csv(os.path.join(d7, "unitsH3", f)) for f in files], ignore_index=True) if files else None
        xt, xl, x = [], [], 0
        for mk, lab, col in (("hcp", "HCP", PALETTE["blue_main"]), ("single_subsample", "single subs.", PALETTE["teal"])):
            if H is None:
                continue
            gw = H[(H.method == mk) & np.isclose(H.alpha, alpha)].groupby(["partition", "reshuffle"]).coverage.mean()
            q = np.quantile(gw, [0.05, 0.95])
            lows.append(q[0])
            ax.plot([x, x], [q[0], q[1]], color=col, lw=2.2, solid_capstyle="butt")
            ax.plot(x, gw.mean(), "o", color=col, ms=4.5, mec=PALETTE["ink"], mew=0.4, zorder=4)
            xt.append(x)
            xl.append(lab)
            x += 1
        f = os.path.join(d8, "report", f"{t}_hjplus_summary.csv")
        if os.path.exists(f):
            S = pd.read_csv(f)
            s_ = S[np.isclose(S.alpha, alpha)]
            if len(s_):
                ax.plot(x, float(s_.group_weighted.iloc[0]), "D", color=PALETTE["violet"], ms=4.5, mec=PALETTE["ink"],
                        mew=0.4, zorder=4)
                xt.append(x)
                xl.append("HJ+" + ("" if int(s_.n_units.iloc[0]) == 10 else f" ({int(s_.n_units.iloc[0])}/10)"))
                x += 1
        ax.axhline(1 - alpha, color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)))
        ax.axhline(1 - 2 * alpha, color=PALETTE["violet"], lw=0.7, ls=(0, (1, 2)))
        ax.set_xticks(xt, xl, rotation=30, ha="right")
        ax.set_xlim(-0.6, max(x - 0.4, 0.6))
        if j == 0:
            ax.set_ylabel("Group-weighted coverage")
        panel_label(ax, "ab"[j], TNAME[t].split(" (")[0])
    lo_ = min([0.6] + [v - 0.03 for v in lows])
    for a_ in fig.axes:
        a_.set_ylim(lo_, 1.03)
    ax = fig.add_subplot(gs[0, 2])
    scores = list(SCORE_LABEL)
    xs = np.arange(len(scores))
    for k, (t, mk) in enumerate((("formation", "o"), ("bandgap", "s"))):
        P = pd.read_csv(os.path.join(d8, "report", f"{t}_h4_pooled.csv")).set_index("score").loc[scores]
        W = pd.read_csv(os.path.join(d8, "report", f"{t}_h4_within_summary.csv")).set_index("score").loc[scores]
        off = -0.21 + 0.14 * k
        col = PALETTE["blue_main"] if t == "formation" else PALETTE["ochre"]
        ax.errorbar(xs + off, P.auroc, yerr=[P.auroc - P.lo, P.hi - P.auroc], fmt=mk, ms=3.5, color=col, capsize=1.5, lw=0.9)
        ax.errorbar(xs + off + 0.28, W.weighted_auroc, yerr=[W.weighted_auroc - W.lo, W.hi - W.weighted_auroc], fmt=mk,
                    ms=3.5, mfc="white", color=col, capsize=1.5, lw=0.9)
    ax.axhline(0.5, color=PALETTE["grey_mid"], lw=0.7, ls=":")
    ax.axhline(0.7, color=PALETTE["grey_mid"], lw=0.7, ls=(0, (4, 3)))
    ax.set_xticks(xs, [SCORE_LABEL[s] for s in scores], rotation=30, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("AUROC for miscoverage")
    panel_label(ax, "c", "Failure detection (held-out structures)")
    h = [Line2D([], [], color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)), label=f"nominal {1 - alpha:.2f} (a, b); 0.70 (c)"),
         Line2D([], [], color=PALETTE["violet"], lw=0.7, ls=(0, (1, 2)), label=f"1 - 2 alpha = {1 - 2 * alpha:.2f} (a, b)"),
         Line2D([], [], color=PALETTE["grey_dark"], lw=2.2, marker="o", ms=4, mec=PALETTE["ink"],
                label="mean, 5-95 % over partitions x reshuffles (a, b)"),
         Line2D([], [], color=PALETTE["violet"], marker="D", lw=0, ms=4, mec=PALETTE["ink"],
                label="hierarchical jackknife+ (a, b)"),
         Line2D([], [], color=PALETTE["blue_main"], marker="o", lw=0, ms=4, label="formation, pooled (c)"),
         Line2D([], [], color=PALETTE["ochre"], marker="s", lw=0, ms=4, label="band gap, pooled (c)"),
         Line2D([], [], color=PALETTE["grey_dark"], marker="o", mfc="white", lw=0, ms=4, label="open: within group (c)")]
    fig.legend(handles=h, loc="outside lower center", ncol=3, fontsize=6.8)
    return finalize_figure(fig, os.path.join(out, "fig9_jarvis_guarantees_detection")), None


def fig10_jarvis_fewshot(d7, d8, out, alpha=0.10):
    f = os.path.join(d8, "report", "h5_summary.csv")
    if not os.path.exists(f):
        return None, "run fewshot8 and report8 first"
    S = pd.read_csv(f)
    W = pd.read_csv(os.path.join(d8, "fewshot", "worst_groups.csv"))
    targets = [t for t in ("formation", "bandgap") if (W.target == t).any()]
    ncol = int(W.groupby("target").size().max())
    fig = new_figure(DOUBLE_COL_IN, 2.15 * len(targets) + 0.5, wspace=0.05, hspace=0.10)
    gs = fig.add_gridspec(len(targets), ncol)
    letters = iter("abcdefghij")
    for i, t in enumerate(targets):
        for j, (_, w) in enumerate(W[W.target == t].iterrows()):
            ax = fig.add_subplot(gs[i, j])
            for sch, meth, col, ls in (("one_per_group", "split", PALETTE["blue_main"], "-"),
                                       ("one_per_group", "offset_full", PALETTE["green_3"], "-"),
                                       ("clustered", "split", PALETTE["red_strong"], (0, (4, 2)))):
                d = S[(S.target == t) & (S.group == w.group) & (S.scheme == sch) & (S.method == meth) & (~S.null)].sort_values("m")
                if d.empty:
                    continue
                ax.plot(d.m, d.mean_cov, color=col, ls=ls, marker="o", ms=3, lw=1.3)
            d = S[(S.target == t) & (S.group == w.group) & (S.scheme == "one_per_group") & (S.method == "split") & (~S.null)].sort_values("m")
            if len(d):
                ax.plot(d.m, d.nominal, color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)))
            ax.axhline(0.80, color=PALETTE["red_strong"], lw=0.6, ls=(0, (1, 2)))
            ax.axhline(w.primary_coverage, color=PALETTE["grey_dark"], lw=0.8, ls=(0, (1, 1)))
            ax.set_ylim(min(0.38, float(w.primary_coverage) - 0.06), 1.04)
            ax.set_yticks([v for v in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) if v >= ax.get_ylim()[0]])
            ax.set_xticks(sorted(S.m.unique()))
            ax.set_xlabel("Labelled structures m")
            if j == 0:
                ax.set_ylabel(f"Coverage, {TNAME[t].split(' (')[0].lower()}")
            panel_label(ax, next(letters), f"{sg_tex(w.symbol)} (unseen: {w.primary_coverage:.2f})")
    h = [Line2D([], [], color=PALETTE["blue_main"], marker="o", ms=3, label="split, one per composition"),
         Line2D([], [], color=PALETTE["green_3"], marker="o", ms=3, label="offset-corrected, one per composition"),
         Line2D([], [], color=PALETTE["red_strong"], ls=(0, (4, 2)), marker="o", ms=3, label="split, clustered"),
         Line2D([], [], color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)), label="exact mean k/(m+1)"),
         Line2D([], [], color=PALETTE["grey_dark"], lw=0.8, ls=(0, (1, 1)), label="coverage without labels")]
    fig.legend(handles=h, loc="outside lower center", ncol=3, fontsize=6.8)
    return finalize_figure(fig, os.path.join(out, "fig10_jarvis_fewshot")), None


FIGURES8 = [fig8_jarvis_groups, fig9_jarvis_guarantees_detection, fig10_jarvis_fewshot]


def run_figures8(cfg):
    import matplotlib.pyplot as plt
    from .fit7 import dirs7
    from .run8 import dirs8
    d7, d8 = dirs7(cfg)["root"], dirs8(cfg)
    out = d8["figures"]
    os.makedirs(out, exist_ok=True)
    apply_publication_style()
    OVERLAP_LOG.clear()
    done, skipped = {}, {}
    for f in FIGURES8:
        try:
            saved, why = f(d7, d8["root"], out, float(cfg["models7"]["alpha_primary"]))
        except Exception as e:
            plt.close("all")
            saved, why = None, f"error: {type(e).__name__}: {e}"
        if saved:
            name = os.path.splitext(os.path.basename(str(saved[0])))[0]
            issues = OVERLAP_LOG.get(name, [])
            done[f.__name__] = [os.path.basename(str(p)) for p in saved]
            print(f"[figures8] {f.__name__}: {', '.join(done[f.__name__])} | overlap check: "
                  f"{'clean' if not issues else f'{len(issues)} ISSUE(S)'}", flush=True)
            for d in issues:
                print(f"           {d['kind']}: {d['a']}  <->  {d['b']}", flush=True)
        else:
            skipped[f.__name__] = why
            print(f"[figures8] {f.__name__}: skipped ({why})", flush=True)
    with open(os.path.join(out, "FIGURES8_MANIFEST.json"), "w") as fh:
        json.dump(dict(figures=done, skipped=skipped, overlap_issues=OVERLAP_LOG), fh, indent=2, default=str)
    n = sum(len(v) for v in OVERLAP_LOG.values())
    print(f"[figures8] {len(done)} figures in {out}; overlap issues: {n}")
    return done
