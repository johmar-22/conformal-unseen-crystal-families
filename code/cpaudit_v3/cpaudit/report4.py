"""Phase 4 report: detection metrics with cluster-bootstrap CIs, selective prediction, lattice-level
analysis, flag rates, pre-registered decisions, the Phase 2 bias supplement and a figure."""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from .config import atomic_write_csv, atomic_write_json, atomic_write_text, get_paths
from .data import LATTICE_NAME
from .detect import (SCORES, ClusterBootstrap, abstention_curve, auroc, avg_precision,
                     percentile_ci)
from .metrics import clopper_pearson

LABEL = {"width": "CV+ width", "knn": "kNN distance", "kde": "KDE dissimilarity", "committee": "Lattice-fold SD",
         "subsample": "Subsample-ensemble SD", "classifier": "Domain classifier"}


def _scores_for(df):
    return [s for s in SCORES if df[s].notna().all()]


# ----------------------------------------------------------------------------- structure level
def pooled_metrics(L, outcome, B, seed):
    y = L[outcome].astype(bool).to_numpy()
    scores = _scores_for(L)
    S = {s: L[s].to_numpy(float) for s in scores}
    point = {s: auroc(y, S[s]) for s in scores}
    boot = {s: [] for s in scores}
    diff = []
    cb = ClusterBootstrap(L.group.to_numpy(), seed)
    for _ in range(B):
        i = cb.draw()
        a = {s: auroc(y[i], S[s][i]) for s in scores}
        for s in scores:
            boot[s].append(a[s])
        diff.append(max(a["knn"], a["kde"]) - a["width"])
    rows = []
    for s in scores:
        lo, hi = percentile_ci(boot[s])
        rows.append(dict(score=s, auroc=point[s], auroc_lo=lo, auroc_hi=hi, ap=avg_precision(y, S[s]),
                         prevalence=float(y.mean()), n=int(y.size), n_pos=int(y.sum())))
    d_point = max(point["knn"], point["kde"]) - point["width"]
    return pd.DataFrame(rows), dict(diff=d_point, diff_lo=percentile_ci(diff)[0], diff_hi=percentile_ci(diff)[1],
                                     best_distance=max(point["knn"], point["kde"]))


def within_metrics(L, outcome, B, seed, min_class):
    scores = _scores_for(L)
    elig, rows = [], []
    for lat, g in L.groupby("heldout", sort=False):
        y = g[outcome].astype(bool).to_numpy()
        ok = y.sum() >= min_class and (~y).sum() >= min_class
        r = dict(heldout=lat, n=len(g), n_pos=int(y.sum()), eligible=bool(ok))
        for s in scores:
            r[s] = auroc(y, g[s].to_numpy(float)) if ok else np.nan
        rows.append(r)
        if ok:
            elig.append(g)
    per = pd.DataFrame(rows)
    if not elig:
        return per, pd.DataFrame([dict(score=s, weighted_auroc=np.nan, lo=np.nan, hi=np.nan) for s in scores])
    w = np.array([len(g) for g in elig], float)
    cbs = [ClusterBootstrap(g.group.to_numpy(), seed + k) for k, g in enumerate(elig)]
    ys = [g[outcome].astype(bool).to_numpy() for g in elig]
    point = {s: float(np.average([auroc(y, g[s].to_numpy(float)) for y, g in zip(ys, elig)], weights=w))
             for s in scores}
    boot = {s: [] for s in scores}
    for _ in range(B):
        idx = [cb.draw() for cb in cbs]
        for s in scores:
            boot[s].append(np.average([auroc(y[i], g[s].to_numpy(float)[i]) for y, g, i in zip(ys, elig, idx)],
                                      weights=w))
    summ = pd.DataFrame([dict(score=s, weighted_auroc=point[s], lo=percentile_ci(boot[s])[0],
                              hi=percentile_ci(boot[s])[1], n_lattices=len(elig)) for s in scores])
    return per, summ


def abstention(L, grid, seed):
    rows = []
    cov = ~L.miscovered.astype(bool).to_numpy()
    err = L.abs_err.to_numpy()
    for s in _scores_for(L) + ["oracle_abs_err"]:
        sc = err if s == "oracle_abs_err" else L[s].to_numpy(float)
        for f, c, e in abstention_curve(sc, cov, err, grid, seed):
            rows.append(dict(score=s, abstained=f, coverage_retained=c, mae_retained=e))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- lattice level
def perm_spearman(x, y, n_perm: int, rng):
    """Spearman rho with a two-sided permutation p-value: exact over all n! pairings for n <= 8,
    otherwise Monte Carlo with the (1 + count) / (n_perm + 1) correction (Phipson & Smyth 2010)."""
    import itertools
    rx, ry = rankdata(x), rankdata(y)
    rx, ry = (rx - rx.mean()) / rx.std(), (ry - ry.mean()) / ry.std()
    rho = float(np.mean(rx * ry))
    n = len(rx)
    if n <= 8:
        perms = np.array(list(itertools.permutations(range(n))))
        r = (ry[perms] * rx).mean(axis=1)
        return rho, float(np.mean(np.abs(r) >= abs(rho) - 1e-12))
    P = np.argsort(rng.random((n_perm, n)), axis=1)
    r = (ry[P] * rx).mean(axis=1)
    return rho, float((1 + np.sum(np.abs(r) >= abs(rho) - 1e-12)) / (n_perm + 1))


def lattice_level(T, n_perm, seed):
    rows = []
    for (design, target, lat), g in T.groupby(["design", "target", "heldout"], sort=False):
        rows.append(dict(design=design, target=target, heldout=lat, n=len(g),
                         miscoverage=float(g.miscovered.mean()), abs_bias=float(abs(g.err.mean())),
                         width_median=float(g.width.median()), knn_mean=float(g.knn.mean()),
                         kde_mean=float(g.kde.mean()), classifier_mean=float(g.classifier.mean()),
                         committee_mean=float(g.committee.mean()), flagged=float(g.flagged.mean())))
    C = pd.DataFrame(rows)
    cors = []
    rng = np.random.default_rng(seed)
    for target in ("bandgap", "formation"):
        for subset, D in (("all 11", C[C.target == target]), ("LOSO 6", C[(C.target == target) & (C.design == "loso")])):
            for sc in ("width_median", "knn_mean", "kde_mean", "classifier_mean", "committee_mean"):
                for out in ("miscoverage", "abs_bias"):
                    rho, p = perm_spearman(D[sc].to_numpy(float), D[out].to_numpy(float), n_perm, rng)
                    cors.append(dict(target=target, configs=subset, score=sc, outcome=out, n=len(D),
                                     spearman=rho, p_perm=p))
    return C, pd.DataFrame(cors)


# ----------------------------------------------------------------------------- Phase 2 supplement
def phase2_bias_supplement(paths, cfg):
    a = float(cfg["alpha_primary"])
    rows = []
    for design, gs in (("loso", (227, 12, 33, 167, 194, 206)), ("l2o", (12, 33, 167, 194, 206))):
        for g in gs:
            for t in ("bandgap", "formation"):
                z = np.load(os.path.join(paths.phase2_loso, f"{t}_{design}{g}.npz"))
                y = z["y_test"]
                base = z["mu_te_l"].mean(0) if design == "loso" else z["mu_te"].mean(0)
                meths = (["cvplus_lattice", "cvplus_random_in", "ncvplus_knn", "cqr_cvplus", "cvplus_lattice_noFd3m_cal"]
                         if design == "loso" else ["cvplus_lattice", "ncvplus_knn"])
                for m in meths:
                    k = f"{m}__a{a:.2f}__lo"
                    if k not in z:
                        continue
                    lo, hi = z[k], z[f"{m}__a{a:.2f}__hi"]
                    pred = z["mu_te_r"].mean(0) if m == "cvplus_random_in" else base
                    e = pred - y
                    rows.append(dict(design=design, target=t, heldout=LATTICE_NAME[g], method=m,
                                     bias=float(e.mean()), mae=float(np.abs(e).mean()),
                                     mae_without_offset=float(np.abs(e - e.mean()).mean()),
                                     coverage=float(np.mean((y >= lo) & (y <= hi))),
                                     miss_below=float(np.mean(y < lo)), miss_above=float(np.mean(y > hi))))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- figure
def figure(results, abst, paths, a):
    sys.path.insert(0, paths.drive_root)
    try:
        import figstyle as fs
        fig, axes = fs.fig_double(0.75, 2, 2)
        save = lambda f, name: fs.save_fig(f, name, outdir=paths.report)
        panel = fs.panel
    except Exception as e:
        print(f"  figstyle not usable ({e}); using matplotlib defaults")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.3), constrained_layout=True)
        save = lambda f, name: f.savefig(os.path.join(paths.report, name + ".png"), dpi=600)
        panel = lambda ax, l: ax.text(0, 1.02, f"({l})", transform=ax.transAxes)
    axes = np.asarray(axes).reshape(2, 2)
    letters = iter("abcd")
    for r, target in enumerate(("bandgap", "formation")):
        ax = axes[r, 0]
        P = results[target]["pooled"]["miscovered"].set_index("score")
        W = results[target]["within_summary"]["miscovered"].set_index("score")
        sc = [s for s in SCORES if s in P.index]
        x = np.arange(len(sc))
        ax.errorbar(x - 0.12, P.loc[sc, "auroc"], yerr=[P.loc[sc, "auroc"] - P.loc[sc, "auroc_lo"],
                                                       P.loc[sc, "auroc_hi"] - P.loc[sc, "auroc"]],
                    fmt="o", ms=3, capsize=2, label="pooled over LOSO holdouts")
        if W["weighted_auroc"].notna().any():
            ax.errorbar(x + 0.12, W.loc[sc, "weighted_auroc"], yerr=[W.loc[sc, "weighted_auroc"] - W.loc[sc, "lo"],
                                                                    W.loc[sc, "hi"] - W.loc[sc, "weighted_auroc"]],
                        fmt="s", ms=3, capsize=2, label="within held-out lattice")
        ax.axhline(0.5, color="0.5", lw=0.6, ls=":")
        ax.axhline(0.7, color="0.5", lw=0.6, ls="--")
        ax.set_xticks(x, [LABEL[s] for s in sc], rotation=35, ha="right")
        ax.set_ylim(0, 1)
        ax.set_ylabel(f"AUROC, {'band gap' if target == 'bandgap' else 'formation'}")
        if r == 0:
            ax.legend(loc="best", fontsize=6, frameon=False)
        panel(ax, next(letters))
        ax = axes[r, 1]
        A = abst[target]
        for s in ["width", "knn", "kde", "committee", "oracle_abs_err"]:
            d = A[A.score == s]
            ax.plot(d.abstained, d.coverage_retained, lw=1.0, ls="--" if s == "oracle_abs_err" else "-",
                    label="oracle (|error|)" if s == "oracle_abs_err" else LABEL[s])
        ax.axhline(1 - a, color="0.5", lw=0.6, ls=":")
        ax.set_xlabel("Fraction of structures abstained")
        ax.set_ylabel("Retained coverage")
        ax.set_ylim(0, 1.02)
        if r == 0:
            ax.legend(loc="lower right", fontsize=6, frameon=False)
        panel(ax, next(letters))
    save(fig, "fig_phase4_detection")


# ----------------------------------------------------------------------------- driver
def run_report4(cfg):
    paths = get_paths(cfg)
    f = os.path.join(paths.scores, "phase4_structures.csv")
    if not os.path.exists(f):
        raise RuntimeError("Run phase4 first.")
    T = pd.read_csv(f)
    flags = pd.read_csv(os.path.join(paths.scores, "phase4_flags_by_config.csv"))
    d, r = cfg["detect"], cfg["decisions"]
    a = float(cfg["alpha_primary"])
    seed = int(cfg["seed_base"]) + 4000
    results, abst, decisions = {}, {}, {}
    out_rows = []
    for target in ("bandgap", "formation"):
        L = T[(T.design == "loso") & (T.target == target)].reset_index(drop=True)
        res = dict(pooled={}, within={}, within_summary={}, h4=None)
        for outcome in ("miscovered", "large_error"):
            P, dd = pooled_metrics(L, outcome, d["bootstrap"], seed)
            Wp, Ws = within_metrics(L, outcome, d["bootstrap"], seed + 1, d["min_class"])
            res["pooled"][outcome], res["within"][outcome], res["within_summary"][outcome] = P, Wp, Ws
            if outcome == "miscovered":
                res["h4"] = dd
            for _, row in P.iterrows():
                out_rows.append(dict(target=target, outcome=outcome, level="pooled", **row.to_dict()))
            for _, row in Ws.iterrows():
                out_rows.append(dict(target=target, outcome=outcome, level="within (n-weighted)",
                                     score=row.score, auroc=row.weighted_auroc, auroc_lo=row.lo, auroc_hi=row.hi))
            atomic_write_csv(os.path.join(paths.report, f"phase4_within_{target}_{outcome}.csv"), Wp)
        results[target] = res
        abst[target] = abstention(L, d["abstention_grid"], seed + 2)
        h = res["h4"]
        P = res["pooled"]["miscovered"].set_index("score")
        Ws = res["within_summary"]["miscovered"].set_index("score")
        h4 = bool(h["diff"] >= r["h4_min_diff"] and h["diff_lo"] > 0 and h["best_distance"] >= r["h4_min_auroc"])
        g4 = bool((P.auroc < r["g4_auroc"]).all() and (Ws.weighted_auroc.fillna(0) < r["g4_auroc"]).all())
        decisions[target] = dict(H4_supported=h4, H4_detail={k: round(float(v), 4) for k, v in h.items()},
                                 G4_triggered=g4, pooled_auroc=P.auroc.round(3).to_dict(),
                                 within_weighted_auroc=Ws.weighted_auroc.round(3).to_dict())
    metrics = pd.DataFrame(out_rows)
    atomic_write_csv(os.path.join(paths.report, "phase4_detection_metrics.csv"), metrics)
    atomic_write_csv(os.path.join(paths.report, "phase4_abstention.csv"),
                     pd.concat([v.assign(target=k) for k, v in abst.items()], ignore_index=True))
    C, cors = lattice_level(T, d["lattice_permutations"], seed + 3)
    atomic_write_csv(os.path.join(paths.report, "phase4_lattice_level.csv"), C)
    atomic_write_csv(os.path.join(paths.report, "phase4_lattice_correlations.csv"), cors)
    fl = flags.copy()
    ci = [clopper_pearson(int(round(x * n)), int(n)) for x, n in zip(fl.false_flag_rate, fl.n_idtest)]
    fl["ffr_ci_lo"], fl["ffr_ci_hi"] = [c[0] for c in ci], [c[1] for c in ci]
    prec = (T.groupby(["design", "target", "heldout", "flagged"]).miscovered.mean().unstack("flagged")
            .rename(columns={True: "miscov_if_flagged", False: "miscov_if_not_flagged"}).reset_index())
    atomic_write_csv(os.path.join(paths.report, "phase4_flags.csv"), fl)
    atomic_write_csv(os.path.join(paths.report, "phase4_flag_precision.csv"), prec)
    sup = phase2_bias_supplement(paths, cfg)
    atomic_write_csv(os.path.join(paths.report, "phase2_bias_supplement.csv"), sup)
    atomic_write_json(os.path.join(paths.report, "phase4_decisions.json"), decisions)
    try:
        figure(results, abst, paths, a)
    except Exception as e:
        print(f"  figure skipped: {e}")

    lines = [f"# Phase 4 summary (primary outcome: miscoverage by {cfg['reference_interval']}, alpha = {a})", ""]
    for target in ("bandgap", "formation"):
        res = results[target]
        for outcome in ("miscovered", "large_error"):
            P = res["pooled"][outcome].copy()
            Ws = res["within_summary"][outcome].set_index("score")
            P["within_weighted"] = P.score.map(Ws.weighted_auroc)
            P["within_lo"] = P.score.map(Ws.lo)
            P["within_hi"] = P.score.map(Ws.hi)
            lines += [f"## {target}: {outcome} (pooled over the 6 LOSO holdouts; 95 % cluster-bootstrap CI)", "",
                      P.round(3).to_markdown(index=False), "",
                      "Within-lattice AUROC per holdout (eligible if >= "
                      f"{d['min_class']} structures in each class):", "",
                      res["within"][outcome].round(3).to_markdown(index=False), ""]
        A = abst[target]
        piv = A[A.abstained.isin([0.0, 0.1, 0.25, 0.5])].pivot_table(index="score", columns="abstained",
                                                                       values="coverage_retained")
        lines += [f"## {target}: coverage of retained structures by abstention fraction", "",
                  piv.round(3).to_markdown(), ""]
    lines += ["## Flag (kNN-distance conformal p-values, p <= "
              f"{d['flag']['alpha']})", "", fl.round(3).to_markdown(index=False), "",
              "Miscoverage among flagged / unflagged held-out structures:", "", prec.round(3).to_markdown(index=False), "",
              "## Lattice level (descriptive; 11 configurations reuse 6 test lattices)", "",
              C.round(3).to_markdown(index=False), "", cors.round(3).to_markdown(index=False), "",
              "## Phase 2 supplement: bias and one-sided misses", "",
              sup[sup.method == "cvplus_lattice"].round(3).to_markdown(index=False), "",
              "## Pre-registered decisions", "", "```", json.dumps(decisions, indent=2), "```", ""]
    atomic_write_text(os.path.join(paths.report, "PHASE4_SUMMARY.md"), "\n".join(lines))
    print("\n".join(lines))
    print(f"\nReport written to {paths.report}")
