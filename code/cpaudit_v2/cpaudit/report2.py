"""Phase 2/3 report: tables, pre-registered decisions, conditional weighted CV+, figure."""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

from . import conformal as cp
from .config import atomic_write_csv, atomic_write_json, atomic_write_text, get_paths
from .data import LATTICE_NAME, LATTICE_ORDER
from .metrics import summarize

REF_NAME = "Fd-3m"
METHODS = ["cvplus_lattice", "cvplus_random_in", "ncvplus_knn", "cqr_cvplus", "cvplus_lattice_noFd3m_cal", "ncvplus_knn_noFd3m_cal", "cqr_cvplus_noFd3m_cal", "weighted_cvplus"]
LABEL = {"cvplus_lattice": "CV+ (lattice folds)", "cvplus_random_in": "CV+ (random folds, in-distribution)",
         "ncvplus_knn": "Normalized CV+ (kNN distance)", "cqr_cvplus": "CQR-CV+",
         "cvplus_lattice_noFd3m_cal": "CV+ without Fd-3m calibration terms"}


def _load(pattern):
    files = sorted(glob.glob(pattern))
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame()


def weighted_rows(paths, cfg, overlap: pd.DataFrame) -> pd.DataFrame:
    """Weighted (JAW-style) lattice-fold CV+, only for holdouts that passed the overlap gate."""
    ok = overlap[overlap.get("overlap_ok", False) == True]  # noqa: E712
    if ok.empty:
        return pd.DataFrame()
    wz = np.load(os.path.join(paths.overlap, "overlap_weights.npz"))
    wid = pd.Series(np.arange(len(wz["id"])), index=wz["id"])
    rows = []
    for _, o in ok.iterrows():
        g = int(o.sg)
        for f in sorted(glob.glob(os.path.join(paths.loso, f"*_loso{g}.npz"))):
            target = os.path.basename(f).split("_loso")[0]
            z = np.load(f)
            w = wz[f"w_{g}"]
            w_tr, w_te = w[wid.loc[z["train_id"]].to_numpy()], w[wid.loc[z["test_id"]].to_numpy()]
            mubar = z["mu_te_l"].mean(axis=0)
            for a in cfg["alpha"]:
                lo, hi = cp.weighted_cvplus_interval(z["mu_te_l"], z["fold_l"], z["R_l"], w_tr, w_te, float(a))
                r = dict(target=target, design="loso", heldout=LATTICE_NAME[g], sg=g, method="weighted_cvplus",
                         alpha=float(a), group="heldout")
                r.update(summarize(z["y_test"], mubar, lo, hi, float(a)))
                rows.append(r)
    return pd.DataFrame(rows)


def decisions(loso, l2o, comp, a_p) -> dict:
    out = {}
    sel = lambda d, m, g="heldout": d[(d.method == m) & (d.alpha.round(4) == round(a_p, 4)) & (d.group == g)]
    for tgt in ("bandgap", "formation"):
        res = {}
        # P2.1 calibration composition
        c = comp[(comp.target == tgt) & (comp.folds == "lattice") & (comp.alpha.round(4) == round(a_p, 4))
                 & (comp.cal_lattice == REF_NAME) & (comp.heldout != REF_NAME)]
        if not c.empty:
            ratios = dict(zip(c.heldout, c.ratio.round(3)))
            n_ok = int((c.ratio >= 2).sum())
            res["P2.1_fd3m_ratio_by_holdout"] = ratios
            res["P2.1_supported"] = bool(n_ok >= 4)
        # H2 ablation and leave-two-out
        for label, d_alt, m_alt in (("ablation", loso, "cvplus_lattice_noFd3m_cal"),
                                    ("leave_two_out", l2o, "cvplus_lattice")):
            if d_alt.empty:
                continue
            base = sel(loso[loso.target == tgt], "cvplus_lattice").set_index("heldout")
            alt = sel(d_alt[d_alt.target == tgt], m_alt).set_index("heldout")
            common = [x for x in LATTICE_ORDER if x != REF_NAME and x in base.index and x in alt.index]
            if not common:
                continue
            dcov = (alt.loc[common, "coverage"] - base.loc[common, "coverage"])
            r = (alt.loc[common, "width_median"] / base.loc[common, "width_median"])
            sup = int(((r <= 0.80) & (dcov <= -0.05)).sum())
            rej = int(((r >= 0.90) & (dcov.abs() <= 0.03)).sum())
            res[f"H2_{label}"] = dict(delta_coverage=dcov.round(3).to_dict(), width_ratio=r.round(3).to_dict(),
                                      n_support=sup, n_reject=rej,
                                      decision=("H2 supported" if sup >= 4 else
                                                "H2 rejected (gate G2)" if rej >= 4 else "partial / mixed"))
        # G3
        d = loso[(loso.target == tgt) & (loso.heldout == REF_NAME)]
        if not d.empty:
            cov = {m: float(sel(d, m).coverage.iloc[0]) for m in ("cvplus_lattice", "ncvplus_knn", "cqr_cvplus")
                   if not sel(d, m).empty}
            wid = {m: float(sel(d, m).width_median.iloc[0]) for m in cov}
            res["G3_fd3m_coverage"] = cov
            res["G3_fd3m_width_median"] = wid
            res["G3_positive"] = bool(max(cov.get("ncvplus_knn", 0), cov.get("cqr_cvplus", 0)) >= 0.85)
            sub = {s: float(sel(d, "cvplus_lattice", s).coverage.iloc[0])
                   for s in ("Fd-3m (gamma~30)", "Fd-3m (gamma~120)") if not sel(d, "cvplus_lattice", s).empty}
            res["T2.3_subfamily_coverage"] = sub
            if len(sub) == 2:
                lo_s, hi_s = min(sub, key=sub.get), max(sub, key=sub.get)
                res["T2.3_one_subfamily_carries_failure"] = bool(sub[lo_s] < 0.75 and sub[hi_s] >= 0.85)
        out[tgt] = res
    return out


def figure(loso, paths, a_p):
    sys.path.insert(0, paths.drive_root)
    try:
        import figstyle as fs
        fig, axes = fs.fig_double(0.75, 2, 1)
        style = lambda i: fs.series(i, line=False)
        panel = fs.panel
        save = lambda f, name: fs.save_fig(f, name, outdir=paths.report)
    except Exception as e:
        print(f"  figstyle not usable ({e}); using matplotlib defaults")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(5.9, 5.3), constrained_layout=True)
        style = lambda i: dict(marker="osD^v"[i % 5], linestyle="none")
        panel = lambda ax, l: ax.text(0, 1.02, f"({l})", transform=ax.transAxes)
        save = lambda f, name: f.savefig(os.path.join(paths.report, name + ".png"), dpi=600)
    methods = ["cvplus_lattice", "cvplus_lattice_noFd3m_cal", "cvplus_random_in", "ncvplus_knn", "cqr_cvplus"]
    d = loso[(loso.alpha.round(4) == round(a_p, 4)) & (loso.group == "heldout")]
    for ax, tgt, letter in zip(axes, ("bandgap", "formation"), "ab"):
        x0 = np.arange(len(LATTICE_ORDER))
        for i, m in enumerate(methods):
            g = d[(d.target == tgt) & (d.method == m)].set_index("heldout")
            xs = [j for j, lat in enumerate(LATTICE_ORDER) if lat in g.index]
            ys = [float(g.loc[LATTICE_ORDER[j], "coverage"]) for j in xs]
            if xs:
                ax.plot(np.asarray(xs) + (i - 2) * 0.13, ys, label=LABEL.get(m, m), **style(i))
        ax.axhline(1 - a_p, color="0.4", lw=0.7, ls="--")
        ax.set_xticks(x0, LATTICE_ORDER)
        ax.set_ylabel(f"{'Band gap' if tgt == 'bandgap' else 'Formation energy'} coverage")
        ax.set_ylim(0, 1.03)
        panel(ax, letter)
    placer = getattr(sys.modules.get("figstyle"), "place_legend", None)
    if placer is not None:
        placer(axes[1], prefer="lower left", ncol=2)
    else:
        axes[1].legend(ncol=2, loc="lower left", fontsize=6)
    fig.text(0.5, -0.01, "Held-out lattice (LOSO); dashed: nominal 0.90", ha="center", va="top", fontsize=7)
    save(fig, "fig_phase2_loso_coverage")


def run_report2(cfg):
    paths = get_paths(cfg)
    a_p = float(cfg["alpha_primary"])
    loso = _load(os.path.join(paths.loso, "*_loso*[0-9].csv"))
    if loso.empty:
        raise RuntimeError(f"No Phase 2 results in {paths.loso}")
    l2o = _load(os.path.join(paths.loso, "*_l2o*[0-9].csv"))
    comp = _load(os.path.join(paths.loso, "*_calcomp.csv"))
    grp = _load(os.path.join(paths.loso, "*_grouplevel.csv"))
    ov_f = os.path.join(paths.overlap, "overlap_diagnostics.csv")
    overlap = pd.read_csv(ov_f) if os.path.exists(ov_f) else pd.DataFrame()
    wrows = weighted_rows(paths, cfg, overlap) if not overlap.empty else pd.DataFrame()
    if not wrows.empty:
        loso = pd.concat([loso, wrows], ignore_index=True)
    atomic_write_csv(os.path.join(paths.report, "phase2_loso_all.csv"), loso)
    if not l2o.empty:
        atomic_write_csv(os.path.join(paths.report, "phase2_l2o_all.csv"), l2o)
    atomic_write_csv(os.path.join(paths.report, "phase2_calibration_composition.csv"), comp)
    atomic_write_csv(os.path.join(paths.report, "phase2_group_level.csv"), grp)
    dec = decisions(loso, l2o, comp, a_p)
    atomic_write_json(os.path.join(paths.report, "phase2_decisions.json"), dec)
    try:
        figure(loso, paths, a_p)
    except Exception as e:
        print(f"  figure skipped: {e}")

    lines = [f"# Phase 2 / Phase 3 summary (alpha = {a_p})", ""]
    if not overlap.empty:
        lines += ["## Phase 3 overlap diagnostics (computed before LOSO)", "",
                  overlap.round(3).to_markdown(index=False), ""]
    order_m = [m for m in METHODS if m in set(loso.method)]
    for tgt in ("bandgap", "formation"):
        d = loso[(loso.target == tgt) & (loso.alpha.round(4) == round(a_p, 4)) & (loso.group == "heldout")]
        if d.empty:
            continue
        idx = [x for x in LATTICE_ORDER if x in set(d.heldout)]
        cov = d.pivot_table(index="heldout", columns="method", values="coverage").reindex(idx)
        wid = d.pivot_table(index="heldout", columns="method", values="width_median").reindex(idx)
        mae = d.groupby("heldout")[["mae_lattice_ensemble", "mae_random_ensemble"]].first().reindex(idx)
        lines += [f"## {tgt}: LOSO coverage of the held-out lattice", "",
                  cov[[m for m in order_m if m in cov.columns]].round(3).to_markdown(), "",
                  f"## {tgt}: LOSO median width", "",
                  wid[[m for m in order_m if m in wid.columns]].round(3).to_markdown(), "",
                  f"## {tgt}: LOSO MAE (fold-model ensembles)", "", mae.round(4).to_markdown(), ""]
        dn = loso[(loso.target == tgt) & (loso.alpha.round(4) == round(a_p, 4)) & (loso.group == "heldout_noconflict")]
        if not dn.empty:
            lines += [f"## {tgt}: coverage with conflicting records excluded from evaluation", "",
                      dn.pivot_table(index="heldout", columns="method", values="coverage").reindex(idx)
                      [[m for m in order_m if m in set(dn.method)]].round(3).to_markdown(), ""]
        if not l2o.empty:
            e = l2o[(l2o.target == tgt) & (l2o.alpha.round(4) == round(a_p, 4)) & (l2o.group == "heldout")]
            if not e.empty:
                piv = e.pivot_table(index="heldout", columns="method", values=["coverage", "width_median"])
                piv = piv.reindex([x for x in LATTICE_ORDER if x in piv.index])
                lines += [f"## {tgt}: leave-two-out (Fd-3m removed from training and calibration)", "",
                          piv.round(3).to_markdown(), ""]
        c = comp[(comp.target == tgt) & (comp.folds == "lattice") & (comp.alpha.round(4) == round(a_p, 4))]
        if not c.empty:
            lines += [f"## {tgt}: calibration composition, share of residuals >= pooled quantile / share of points",
                      "", c.pivot_table(index="heldout", columns="cal_lattice", values="ratio")
                      .reindex(idx)[[x for x in LATTICE_ORDER if x in set(c.cal_lattice)]].round(2).to_markdown(), ""]
        gl = grp[grp.target == tgt] if not grp.empty else grp
        if not gl.empty:
            lines += [f"## {tgt}: group-level conformal (K = 5 families, illustrative; not in the tables above)", "",
                      gl[["heldout", "alpha", "k_needed", "threshold", "coverage_of_heldout", "family_event_covered"]]
                      .round(3).to_markdown(index=False), ""]
    lines += ["## Pre-registered decisions", "", "```", json.dumps(dec, indent=2), "```", ""]
    atomic_write_text(os.path.join(paths.report, "PHASE2_SUMMARY.md"), "\n".join(lines))
    print("\n".join(lines))
    print(f"\nReport written to {paths.report}")
