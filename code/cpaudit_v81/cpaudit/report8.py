"""Phase 6 Task 6.5 report (v8): per-group table, across-group statistics, H3 (valid design) and
the hierarchical jackknife+, H4-J (failure detection), H5-J (few-shot), H2 (descriptive), oxide
subset, and the side-by-side H6 table. Reads phase07_v1 and phase08_v1; writes phase08_v1/report.
"""
from __future__ import annotations

import glob
import json
import math
import os

import numpy as np
import pandas as pd

from . import cvfast as cf
from .an8 import boot_spearman, perm_spearman_mc, spearman, tail_enrichment
from .config import atomic_write_csv, atomic_write_json, atomic_write_text
from .detect import ClusterBootstrap, abstention_curve, auroc, percentile_ci
from .fewshot import nominal_mean
from .metrics import clopper_pearson
from .run8 import DETECT_TARGETS, Ctx8, primary_group_coverage

SCORES = ["width", "knn", "kde", "committee", "classifier"]
LABEL = {"width": "CV+ width", "knn": "kNN distance", "kde": "KDE dissimilarity",
         "committee": "inner-model SD", "classifier": "domain classifier"}
METRICS = ["coverage", "width_median", "width_median_sd", "mae", "bias", "interval_score"]


def _md(df, digits=3):
    return df.round(digits).to_markdown(index=False) if len(df) else "(none)"


# ----------------------------------------------------------------------------- per-group table
def unit_csvs(ctx, target, design_rep):
    n_out = 2 * ctx.m["outer_folds_per_half"]
    out = []
    for design, rep in design_rep:
        if design == "rep2" and "rep2" not in ctx.m["blocks"][target]:
            continue
        for o in range(n_out):
            d = pd.read_csv(os.path.join(ctx.d["unitsB"], f"{target}_{design}_r{rep}_o{o:02d}.csv"))
            out.append(d.assign(rep=rep, outer=o))
    return pd.concat(out, ignore_index=True)


def group_distances(ctx, target):
    """Median kNN distance (and KDE dissimilarity, repeat 1) of each held-out group, per repeat."""
    rows = []
    for f in sorted(glob.glob(os.path.join(ctx.d8["scores"], f"{target}_r*_o*.npz"))):
        rep = int(os.path.basename(f).split("_r")[1].split("_")[0])
        with np.load(f) as z:
            g = ctx.L["sga01_number"].to_numpy()[z["test_rows"]]
            for u in np.unique(g):
                m = g == u
                rows.append(dict(group=int(u), rep=rep, knn_median=float(np.median(z["knn"][m])),
                                 kde_median=float(np.median(z["kde"][m])) if "kde" in z.files else np.nan))
    return pd.DataFrame(rows)


def group_table(ctx, target):
    a = float(ctx.m["alpha_primary"])
    T = unit_csvs(ctx, target, (("primary", 1), ("rep2", 2)))
    T = T[(T.group != "ALL") & (T.method == "cvplus") & np.isclose(T.alpha, a)].copy()
    T["group"] = T.group.astype(int)
    G = T.groupby("group").agg(**{k: (k, "mean") for k in METRICS}, n=("n", "first"),
                               sd_y=("sd_y", "mean"), n_repeats=("rep", "nunique")).reset_index()
    Dd = group_distances(ctx, target)
    if len(Dd):
        dist = Dd.groupby("group").agg(knn_median=("knn_median", "mean"),
                                       n_distance_repeats=("rep", "nunique")).reset_index()
        kde = Dd[Dd.rep == 1][["group", "kde_median"]]
        G = G.merge(dist, on="group", how="left").merge(kde, on="group", how="left")
    else:
        G["knn_median"], G["kde_median"], G["n_distance_repeats"] = np.nan, np.nan, 0
    G["offset_share"] = G.bias.abs() / G.mae
    G.insert(1, "symbol", G.group.map(lambda g: ctx.sgsym.get(int(g), str(g))))
    # consistency with the primary coverage used by v7 / fewshot8
    P = primary_group_coverage(ctx, target).set_index("group")
    if not np.allclose(G.set_index("group").coverage.loc[P.index], P.coverage, atol=1e-12):
        raise AssertionError("group coverage differs from primary_group_coverage")
    return G.sort_values("coverage").reset_index(drop=True)


def across_groups(G, B, n_perm, seed):
    out = dict(n_groups=int(len(G)), share_below_090=float((G.coverage < 0.90).mean()),
               share_below_080=float((G.coverage < 0.80).mean()),
               median_coverage=float(G.coverage.median()))
    for col in ("knn_median", "kde_median"):
        x = G[col].to_numpy(float)
        ok = np.isfinite(x)
        if ok.sum() < 5:
            continue
        rho, lo, hi, nb = boot_spearman(G.coverage.to_numpy()[ok], x[ok], B, seed)
        _, p = perm_spearman_mc(G.coverage.to_numpy()[ok], x[ok], n_perm, seed + 1)
        out[f"spearman_coverage_{col}"] = dict(rho=rho, ci95=[lo, hi], p_perm=p, n=int(ok.sum()), n_boot=nb)
    x = G.knn_median.to_numpy(float)
    ok = np.isfinite(x)
    if ok.sum() >= 5:
        out["spearman_mae_knn_median"] = spearman(G.mae.to_numpy()[ok], x[ok])
        out["spearman_coverage_offset_share"] = spearman(G.coverage, G.offset_share)
    return out


# ----------------------------------------------------------------------------- H3 and HJ+
def h3_summary(ctx, target):
    files = sorted(glob.glob(os.path.join(ctx.d["unitsH3"], f"{target}_H3_p*.csv")))
    if not files:
        return pd.DataFrame()
    D = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    D["inf"] = ~np.isfinite(D.threshold)
    rows = []
    for (meth, a), g in D.groupby(["method", "alpha"]):
        gw = g.groupby(["partition", "reshuffle"]).coverage.mean()
        rows.append(dict(method=meth, alpha=a, group_weighted_mean=float(gw.mean()), group_weighted_sd=float(gw.std()),
                         share_groups_at_nominal=float((g.coverage >= 1 - a - 1e-12).mean()),
                         share_groups_below_080=float((g.coverage < 0.80).mean()),
                         share_unbounded=float(g.inf.mean()), n_test_group_instances=int(len(g)),
                         n_partitions=int(g.partition.nunique())))
    return pd.DataFrame(rows)


def hjplus_summary(ctx, target):
    n_out = 2 * ctx.m["outer_folds_per_half"]
    files = [os.path.join(ctx.d8["hjplus"], f"{target}_hj_o{o:02d}.csv") for o in range(n_out)]
    have = [f for f in files if os.path.exists(f)]
    if not have:
        return None
    H = pd.concat([pd.read_csv(f) for f in have], ignore_index=True)
    a = float(ctx.m["alpha_primary"])
    V = unit_csvs(ctx, target, (("primary", 1),))
    V = V[(V.group != "ALL") & np.isclose(V.alpha, a) & V.method.isin(["cvplus", "cvplus_gw"])]
    V = V[V.outer.isin([int(os.path.basename(f)[-6:-4]) for f in have])]
    rows, per = [], []
    for al, g in H[H.group != "ALL"].groupby("alpha"):
        rows.append(dict(alpha=al, n_units=len(have), n_groups=len(g),
                         marginal=float(np.average(g.coverage, weights=g.n)),
                         group_weighted=float(g.coverage.mean()),
                         share_groups_at_nominal=float((g.coverage >= 1 - al - 1e-12).mean()),
                         share_groups_below_080=float((g.coverage < 0.80).mean()),
                         median_width_sd=float(g.width_median_sd.median()),
                         n_unbounded=int(g.n_unbounded.sum()), theorem5_bound=1 - 2 * al))
    g = H[(H.group != "ALL") & np.isclose(H.alpha, a)][["group", "coverage", "width_median_sd", "n"]].copy()
    g["group"] = g.group.astype(int)
    for meth in ("cvplus", "cvplus_gw"):
        v = V[V.method == meth][["group", "coverage", "width_median_sd"]].copy()
        v["group"] = v.group.astype(int)
        g = g.merge(v.rename(columns={"coverage": f"cov_{meth}", "width_median_sd": f"wsd_{meth}"}), on="group")
    g.insert(1, "symbol", g.group.map(lambda x: ctx.sgsym.get(int(x), str(x))))
    return dict(summary=pd.DataFrame(rows), per_group=g.rename(columns={"coverage": "cov_hjplus",
                                                                         "width_median_sd": "wsd_hjplus"}),
                complete=len(have) == n_out)


# ----------------------------------------------------------------------------- H4
def structure_table(ctx, target):
    a = float(ctx.m["alpha_primary"])
    n_out = 2 * ctx.m["outer_folds_per_half"]
    frames = []
    for o in range(n_out):
        f = os.path.join(ctx.d8["scores"], f"{target}_r1_o{o:02d}.npz")
        if not os.path.exists(f):
            raise RuntimeError(f"{f} missing: run scores8 first")
        U = ctx.unit(target, "primary", 1, o)
        with np.load(f) as z:
            S = {k: z[k] for k in z.files}
        if not np.array_equal(S["test_rows"], U["test_rows"]):
            raise AssertionError(f"score rows differ from unit rows ({target} o{o})")
        mu = U["mu_te"].astype(float)
        lo, hi = U[f"cvplus__a{a:.2f}__lo"].astype(float), U[f"cvplus__a{a:.2f}__hi"].astype(float)
        y = U["y_test"].astype(float)
        pred = mu.mean(axis=0)
        frames.append(pd.DataFrame(dict(
            target=target, outer=o, row=U["test_rows"], group=U["test_grp"].astype(int),
            formula=ctx.formula[U["test_rows"]], y=y, pred=pred, abs_err=np.abs(y - pred), lo=lo, hi=hi,
            miscovered=(y < lo) | (y > hi), width=hi - lo, committee=mu.std(axis=0, ddof=1),
            knn=S["knn"], kde=S["kde"], classifier=S["classifier"], flag_p=S["flag_p_test"])))
    return pd.concat(frames, ignore_index=True)


def h4_pooled(L, B, seed):
    y = L.miscovered.to_numpy(bool)
    Sc = {s: L[s].to_numpy(float) for s in SCORES}
    point = {s: auroc(y, Sc[s]) for s in SCORES}
    boot = {s: [] for s in SCORES}
    diff = []
    cb = ClusterBootstrap(L.group.to_numpy(), seed)          # resample held-out space groups
    for _ in range(B):
        i = cb.draw()
        a_ = {s: auroc(y[i], Sc[s][i]) for s in SCORES}
        for s in SCORES:
            boot[s].append(a_[s])
        diff.append(max(a_["knn"], a_["kde"]) - a_["width"])
    P = pd.DataFrame([dict(score=s, auroc=point[s], lo=percentile_ci(boot[s])[0], hi=percentile_ci(boot[s])[1],
                           n=int(y.size), n_miscovered=int(y.sum())) for s in SCORES])
    d = max(point["knn"], point["kde"]) - point["width"]
    lo, hi = percentile_ci(diff)
    return P, dict(diff=d, diff_lo=lo, diff_hi=hi, best_distance=max(point["knn"], point["kde"]),
                   best_distance_score="knn" if point["knn"] >= point["kde"] else "kde")


def h4_within(L, B, seed, min_class):
    rows, elig = [], []
    for g, d in L.groupby("group", sort=True):
        yy = d.miscovered.to_numpy(bool)
        ok = yy.sum() >= min_class and (~yy).sum() >= min_class
        r = dict(group=int(g), n=len(d), n_miscovered=int(yy.sum()), eligible=bool(ok))
        for s in SCORES:
            r[s] = auroc(yy, d[s].to_numpy(float)) if ok else np.nan
        rows.append(r)
        if ok:
            elig.append(d)
    per = pd.DataFrame(rows)
    if not elig:
        return per, pd.DataFrame([dict(score=s, weighted_auroc=np.nan, lo=np.nan, hi=np.nan, n_groups=0) for s in SCORES])
    w = np.array([len(d) for d in elig], float)
    ys = [d.miscovered.to_numpy(bool) for d in elig]
    Ss = [{s: d[s].to_numpy(float) for s in SCORES} for d in elig]
    cbs = [ClusterBootstrap(d.formula.to_numpy(), seed + k) for k, d in enumerate(elig)]  # formulas within group
    point = {s: float(np.average([auroc(y, S_[s]) for y, S_ in zip(ys, Ss)], weights=w)) for s in SCORES}
    boot = {s: [] for s in SCORES}
    for _ in range(B):
        idx = [cb.draw() for cb in cbs]
        for s in SCORES:
            v = np.array([auroc(y[i], S_[s][i]) for y, S_, i in zip(ys, Ss, idx)])
            ok = np.isfinite(v)
            boot[s].append(np.average(v[ok], weights=w[ok]) if ok.any() else np.nan)
    summ = pd.DataFrame([dict(score=s, weighted_auroc=point[s], lo=percentile_ci(boot[s])[0],
                              hi=percentile_ci(boot[s])[1], n_groups=len(elig)) for s in SCORES])
    return per, summ


def h4_flags(ctx, target, L, alpha_flag):
    n_out = 2 * ctx.m["outer_folds_per_half"]
    rows = []
    for o in range(n_out):
        with np.load(os.path.join(ctx.d8["scores"], f"{target}_r1_o{o:02d}.npz")) as z:
            p_id, p_te, ncal = z["flag_p_idtest"], z["flag_p_test"], int(z["flag_n_cal"])
        k = int(np.sum(p_id <= alpha_flag))
        lo, hi = clopper_pearson(k, len(p_id))
        rows.append(dict(outer=o, n_cal=ncal, n_idtest=len(p_id), false_flag_rate=k / len(p_id), ffr_lo=lo,
                         ffr_hi=hi, n_heldout=len(p_te), flag_rate_heldout=float(np.mean(p_te <= alpha_flag))))
    F = pd.DataFrame(rows)
    fl = L.flag_p <= alpha_flag
    prec = dict(miscoverage_if_flagged=float(L.miscovered[fl].mean()) if fl.any() else np.nan,
                miscoverage_if_not_flagged=float(L.miscovered[~fl].mean()) if (~fl).any() else np.nan,
                share_flagged=float(fl.mean()))
    return F, prec


def h4_decide(dd, P, Ws, rules):
    h4 = bool(dd["diff"] >= rules["delta"] and dd["diff_lo"] > 0 and dd["best_distance"] >= rules["min_distance_auroc"])
    g4 = bool((P.auroc < rules["g4_threshold"]).all() and (Ws.weighted_auroc.fillna(0) < rules["g4_threshold"]).all())
    return dict(H4J_supported=h4, G4J_triggered=g4, detail={k: (round(float(v), 4) if isinstance(v, float) else v)
                                                            for k, v in dd.items()},
                pooled_auroc=P.set_index("score").auroc.round(3).to_dict(),
                within_weighted_auroc=Ws.set_index("score").weighted_auroc.round(3).to_dict())


# ----------------------------------------------------------------------------- H5
def h5_tables(ctx, alpha, rules, seed):
    from .report5 import _boot_ratio
    files = sorted(glob.glob(os.path.join(ctx.d8["fewshot"], "*_g*.csv.gz")))
    if not files:
        return None
    D = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    W = pd.read_csv(os.path.join(ctx.d8["fewshot"], "worst_groups.csv"))
    keys = ["target", "group", "scheme", "null", "m", "method"]
    S = D.groupby(keys).agg(n_draws=("cov", "size"), mean_cov=("cov", "mean"), sd_cov=("cov", "std"),
                            q05=("cov", lambda x: x.quantile(0.05)), p_below_080=("cov", lambda x: float(np.mean(x < 0.80))),
                            width_median=("width", "median"), frac_infinite=("infinite", "mean"),
                            mae_raw=("mae_raw", "mean"), mae_offset=("mae_offset", "mean")).reset_index()
    S = S.merge(W[["target", "group", "symbol", "sd_y"]], on=["target", "group"], how="left")
    S["width_sd"] = S.width_median / S.sd_y
    S["nominal"] = [nominal_mean(int(m), alpha) for m in S.m]
    rng = np.random.default_rng(seed)
    rho_rows = []
    sel = D[D.scheme.isin(["clustered", "one_per_group"]) & (D.method == "split")]
    for (t, g, sch, m), G in sel.groupby(["target", "group", "scheme", "m"]):
        x, y = G[~G.null]["cov"].to_numpy(), G[G.null]["cov"].to_numpy()
        if len(x) < 10 or len(y) < 10:
            continue
        lo, hi = _boot_ratio(x, y, 1000, rng) if m == rules["m"] else (np.nan, np.nan)
        rho_rows.append(dict(target=t, group=g, scheme=sch, m=m, mean_real=x.mean(), mean_null=y.mean(),
                             nominal=nominal_mean(int(m), alpha), sd_real=np.std(x, ddof=1), sd_null=np.std(y, ddof=1),
                             rho=np.std(x, ddof=1) / np.std(y, ddof=1), rho_lo=lo, rho_hi=hi))
    Rr = pd.DataFrame(rho_rows).merge(W[["target", "group", "symbol"]], on=["target", "group"], how="left")
    dec = {}
    for t, Wt in W.groupby("target"):
        cell = Rr[(Rr.target == t) & (Rr.scheme == "one_per_group") & (Rr.m == rules["m"])].set_index("group")
        groups = Wt.group.tolist()
        means = {int(g): (float(cell.loc[g, "mean_real"]) if g in cell.index else np.nan) for g in groups}
        rhos = {int(g): (float(cell.loc[g, "rho"]) if g in cell.index else np.nan) for g in groups}
        n_mean_ok = sum(1 for v in means.values() if np.isfinite(v) and v >= rules["mean_coverage"])
        n_rho_ok = sum(1 for v in rhos.values() if np.isfinite(v) and v <= rules["dispersion_ratio"])
        ok = n_mean_ok == len(groups) and n_rho_ok >= rules["min_groups"]
        dec[t] = dict(H5J_supported=bool(ok), m=rules["m"], groups={int(g): ctx.sgsym.get(int(g), str(g)) for g in groups},
                      mean_coverage=means, dispersion_ratio=rhos, n_mean_ok=n_mean_ok, n_rho_ok=n_rho_ok,
                      evaluable=bool(all(np.isfinite(v) for v in means.values())))
    return dict(summary=S, rho=Rr, decisions=dec, worst=W)


# ----------------------------------------------------------------------------- H2
def h2_tables(ctx, target):
    a = float(ctx.m["alpha_primary"])
    h = ctx.m8["h2"]
    n_out = 2 * ctx.m["outer_folds_per_half"]
    tail_rows, abl_rows = [], []
    for o in range(n_out):
        U = ctx.unit(target, "primary", 1, o)
        y_cal = ctx.y[target][U["cal_rows"]]
        R = np.abs(y_cal - U["oof"].astype(float))
        rows, tau = tail_enrichment(R, U["cal_grp"], h["tail_quantile"])
        for r in rows:
            tail_rows.append(dict(outer=o, tau=tau, **r))
        G = pd.DataFrame(rows).sort_values("oof_mae", ascending=False)
        n_drop = int(math.ceil(h["drop_share"] * len(G) - 1e-9))
        drop = G.group.to_numpy()[:n_drop]
        keep = ~np.isin(U["cal_grp"], drop)
        mu = U["mu_te"].astype(float)
        K = int(U["n_folds"])
        lo, hi = cf.cvplus(mu, mu, None, R[keep], U["cal_fold"].astype(int)[keep], K, a)
        lo0, hi0 = U[f"cvplus__a{a:.2f}__lo"].astype(float), U[f"cvplus__a{a:.2f}__hi"].astype(float)
        y = U["y_test"].astype(float)
        for g in np.unique(U["test_grp"]):
            m = U["test_grp"] == g
            abl_rows.append(dict(outer=o, group=int(g), n=int(m.sum()),
                                 coverage=float(np.mean((y[m] >= lo0[m]) & (y[m] <= hi0[m]))),
                                 coverage_ablated=float(np.mean((y[m] >= lo[m]) & (y[m] <= hi[m]))),
                                 width_ratio=float(np.median(hi[m] - lo[m]) / np.median(hi0[m] - lo0[m])),
                                 dropped=",".join(ctx.sgsym.get(int(d), str(d)) for d in drop)))
    T = pd.DataFrame(tail_rows)
    Tg = T.groupby("group").agg(n_units=("outer", "nunique"), enrichment=("enrichment", "mean"),
                                share_tail=("share_tail", "mean"), share_cal=("share_cal", "mean"),
                                oof_mae=("oof_mae", "mean")).reset_index().sort_values("enrichment", ascending=False)
    Tg.insert(1, "symbol", Tg.group.map(lambda g: ctx.sgsym.get(int(g), str(g))))
    A = pd.DataFrame(abl_rows)
    A.insert(2, "symbol", A.group.map(lambda g: ctx.sgsym.get(int(g), str(g))))
    return Tg, A


# ----------------------------------------------------------------------------- oxides
def oxide_tables(ctx, target):
    a = float(ctx.m["alpha_primary"])
    ox = ctx.F[ctx.m8["oxide"]["rule"]].to_numpy(bool)
    n_out = 2 * ctx.m["outer_folds_per_half"]
    rows = []
    for design, rep in (("primary", 1), ("rep2", 2)):
        if design == "rep2" and "rep2" not in ctx.m["blocks"][target]:
            continue
        for o in range(n_out):
            U = ctx.unit(target, design, rep, o)
            y = U["y_test"].astype(float)
            lo, hi = U[f"cvplus__a{a:.2f}__lo"].astype(float), U[f"cvplus__a{a:.2f}__hi"].astype(float)
            cov = (y >= lo) & (y <= hi)
            is_ox = ox[U["test_rows"]]
            for g in np.unique(U["test_grp"]):
                m = U["test_grp"] == g
                rows.append(dict(rep=rep, group=int(g), n_oxide=int((m & is_ox).sum()), n_other=int((m & ~is_ox).sum()),
                                 cov_oxide=float(cov[m & is_ox].mean()) if (m & is_ox).any() else np.nan,
                                 cov_other=float(cov[m & ~is_ox].mean()) if (m & ~is_ox).any() else np.nan))
    O = pd.DataFrame(rows).groupby("group").agg(n_oxide=("n_oxide", "first"), n_other=("n_other", "first"),
                                                cov_oxide=("cov_oxide", "mean"), cov_other=("cov_other", "mean")).reset_index()
    O.insert(1, "symbol", O.group.map(lambda g: ctx.sgsym.get(int(g), str(g))))
    sel = O[O.n_oxide >= ctx.m8["oxide"]["min_rows"]].sort_values("cov_oxide").reset_index(drop=True)
    summ = dict(n_groups=int(len(sel)), median_cov_oxide=float(sel.cov_oxide.median()) if len(sel) else np.nan,
                share_below_090=float((sel.cov_oxide < 0.90).mean()) if len(sel) else np.nan,
                share_below_080=float((sel.cov_oxide < 0.80).mean()) if len(sel) else np.nan,
                n_oxide_rows=int(O.n_oxide.sum()))
    return sel, summ


# ----------------------------------------------------------------------------- driver
def run_report8(cfg, allow_deviation=False):
    from .setup8 import require_prereg8
    pre = require_prereg8(cfg, allow_deviation)
    ctx = Ctx8(cfg)
    m8_ = ctx.m8
    a = float(ctx.m["alpha_primary"])
    rules = ctx.m["decisions"]
    B, n_perm, seed = int(m8_["bootstrap"]), int(m8_["n_perm"]), int(m8_["seed"])
    rep = ctx.d8["report"]
    lines = ["# Phase 6 Task 6.5 report (JARVIS-DFT; independent replication, not pooled with NOMAD)", "",
             f"Addendum written {pre['time_utc']}; v7 rules registered {pre['v7_preregistration']['time_utc']}.", ""]
    decisions = {}

    for t in ctx.m["targets"]:
        G = group_table(ctx, t)
        atomic_write_csv(os.path.join(rep, f"{t}_groups.csv"), G)
        acr = across_groups(G, B, n_perm, seed + 11)
        decisions.setdefault(t, {})["across_groups"] = acr
        sp = acr.get("spearman_coverage_knn_median", {})
        lines += [f"## {t}: held-out groups (primary CV+, alpha {a}; mean over repeats)", "",
                  f"{acr['n_groups']} groups; median coverage {acr['median_coverage']:.3f}; below 0.90: "
                  f"{acr['share_below_090']:.1%}; below 0.80: {acr['share_below_080']:.1%}.", ""]
        if sp:
            lines += [f"Spearman(group coverage, median kNN distance) = {sp['rho']:.3f} (95 % bootstrap CI over groups "
                      f"{sp['ci95'][0]:.3f} to {sp['ci95'][1]:.3f}; permutation p {sp['p_perm']:.4f}; n = {sp['n']}).", ""]
        lines += [_md(G[["group", "symbol", "n", "coverage", "width_median_sd", "mae", "bias", "offset_share",
                         "interval_score", "knn_median", "kde_median"]]), ""]

    for t in ctx.m["targets"]:
        if "h3" not in ctx.m["blocks"][t]:
            continue
        H = h3_summary(ctx, t)
        atomic_write_csv(os.path.join(rep, f"{t}_h3_nominal.csv"), H)
        decisions[t]["H3_nominal"] = H.to_dict("records")
        lines += [f"## {t}: H3 valid group-level design (v7 runs)", "",
                  "group_weighted_mean = mean over partitions x reshuffles of the mean within-group coverage; "
                  "share_groups_at_nominal = how often a held-out group is covered at >= 1 - alpha.", "", _md(H), ""]

    lines += ["## Hierarchical jackknife+ (Lee, Barber & Willett, arXiv:2306.06342 v4, App. C.2.2, eq. 20)", "",
              "Theorem 5: if the K + 1 groups are hierarchically exchangeable and the fitting algorithm is "
              "hierarchically symmetric, P{Y_test in C(X_test)} >= 1 - 2 alpha. Quote the theorem verbatim from "
              "the paper in the manuscript.", ""]
    for t in m8_["hjplus"]["targets"]:
        HJ = hjplus_summary(ctx, t)
        if HJ is None:
            lines += [f"* {t}: not run yet (hjplus8).", ""]
            continue
        atomic_write_csv(os.path.join(rep, f"{t}_hjplus_summary.csv"), HJ["summary"])
        atomic_write_csv(os.path.join(rep, f"{t}_hjplus_by_group.csv"), HJ["per_group"])
        decisions[t]["HJplus"] = dict(complete=HJ["complete"], summary=HJ["summary"].to_dict("records"))
        pg = HJ["per_group"]
        lines += [f"### {t}{'' if HJ['complete'] else ' (INCOMPLETE: not all outer folds done)'}", "",
                  _md(HJ["summary"]), "",
                  f"Same held-out structures, alpha {a}: median group coverage HJ+ {pg.cov_hjplus.median():.3f}, "
                  f"CV+ {pg.cov_cvplus.median():.3f}, group-weighted CV+ {pg.cov_cvplus_gw.median():.3f}; "
                  f"groups below 0.80: {(pg.cov_hjplus < 0.8).mean():.1%} / {(pg.cov_cvplus < 0.8).mean():.1%} / "
                  f"{(pg.cov_cvplus_gw < 0.8).mean():.1%}.", "", _md(pg.sort_values('cov_hjplus').head(10)), ""]

    h4res = {}
    for t in DETECT_TARGETS:
        L = structure_table(ctx, t)
        L.to_csv(os.path.join(rep, f"{t}_structures.csv.gz"), index=False, compression="gzip")
        P, dd = h4_pooled(L, B, seed + 21)
        Wp, Ws = h4_within(L, B, seed + 22, int(m8_["min_class"]))
        Wp.insert(1, "symbol", Wp.group.map(lambda g: ctx.sgsym.get(int(g), str(g))))
        F, prec = h4_flags(ctx, t, L, float(m8_["flag"]["alpha"]))
        A = []
        cov = ~L.miscovered.to_numpy(bool)
        for s in SCORES + ["oracle_abs_err"]:
            sc = L.abs_err.to_numpy() if s == "oracle_abs_err" else L[s].to_numpy(float)
            for f_, c_, e_ in abstention_curve(sc, cov, L.abs_err.to_numpy(), m8_["abstention_grid"], seed + 23):
                A.append(dict(score=s, abstained=f_, coverage_retained=c_, mae_retained=e_))
        A = pd.DataFrame(A)
        dec = h4_decide(dd, P, Ws, rules["h4"])
        worst = json.load(open(os.path.join(ctx.d8["phase0"], "PREREGISTRATION_phase6_v8.json")))["worst_groups"][t]
        wg = [w["group"] for w in worst]
        dec["worst_groups_within_auroc"] = Wp[Wp.group.isin(wg)].set_index("symbol")[SCORES + ["n", "n_miscovered"]].round(3).to_dict("index")
        dec["flag"] = dict(mean_false_flag_rate=float(F.false_flag_rate.mean()), **prec)
        decisions[t]["H4J"] = dec
        h4res[t] = dict(P=P, Ws=Ws, A=A)
        for name, df in (("h4_pooled", P), ("h4_within_by_group", Wp), ("h4_within_summary", Ws),
                         ("h4_abstention", A), ("h4_flags", F)):
            atomic_write_csv(os.path.join(rep, f"{t}_{name}.csv"), df)
        piv = A[A.abstained.isin([0.0, 0.1, 0.25, 0.5])].pivot_table(index="score", columns="abstained",
                                                                       values="coverage_retained")
        lines += [f"## {t}: H4-J failure detection (outcome: outside the primary CV+ interval, repeat 1)", "",
                  f"Pooled ({int(L.miscovered.sum())} of {len(L)} structures miscovered); 95 % CIs resample held-out space groups:",
                  "", _md(P), "", "Within-group (n-weighted; CIs resample formulas within groups):", "", _md(Ws), "",
                  f"best distance AUROC - width AUROC = {dd['diff']:.3f} (95 % CI {dd['diff_lo']:.3f} to {dd['diff_hi']:.3f}); "
                  f"**H4-J {'supported' if dec['H4J_supported'] else 'not supported'}; G4-J "
                  f"{'triggered' if dec['G4J_triggered'] else 'not triggered'}**."
                  + (" Both rules are met: a distance score beats the width, but no score reaches "
                     f"{rules['h4']['g4_threshold']} pooled or within group, so it is not a usable detector."
                     if dec["H4J_supported"] and dec["G4J_triggered"] else ""), "",
                  "Three worst-covered groups, within-group AUROC:", "",
                  _md(Wp[Wp.group.isin(wg)][["group", "symbol", "n", "n_miscovered"] + SCORES]), "",
                  "Coverage of retained structures after abstaining on the highest-score fraction:", "",
                  piv.round(3).to_markdown(), "",
                  f"Conformal kNN flag (p <= {m8_['flag']['alpha']}): mean in-distribution false-flag rate "
                  f"{F.false_flag_rate.mean():.3f}; held-out flagged {prec['share_flagged']:.1%}; miscoverage flagged "
                  f"{prec['miscoverage_if_flagged']:.3f} vs not flagged {prec['miscoverage_if_not_flagged']:.3f}.", ""]

    H5 = h5_tables(ctx, a, rules["h5"], seed + 31)
    if H5 is not None:
        atomic_write_csv(os.path.join(rep, "h5_summary.csv"), H5["summary"])
        atomic_write_csv(os.path.join(rep, "h5_dispersion.csv"), H5["rho"])
        for t, d in H5["decisions"].items():
            decisions[t]["H5J"] = d
        S = H5["summary"]
        tab = S[(~S.null) & S.scheme.isin(["one_per_group", "random", "clustered"])][
            ["target", "symbol", "scheme", "m", "method", "mean_cov", "sd_cov", "p_below_080", "nominal", "width_sd",
             "frac_infinite"]]
        lines += ["## H5-J few-shot recalibration of the three worst-covered groups (repeat-1 model)", "",
                  _md(H5["worst"]), "", _md(tab), "", "Dispersion ratio vs permutation null (split):", "",
                  _md(H5["rho"][["target", "symbol", "scheme", "m", "mean_real", "nominal", "rho", "rho_lo", "rho_hi"]]), ""]
        for t, d in H5["decisions"].items():
            lines += [f"* {t}: **H5-J {'supported' if d['H5J_supported'] else 'not supported'}** "
                      f"(m = {d['m']}: mean coverage {json.dumps({d['groups'][g]: round(v, 3) for g, v in d['mean_coverage'].items()})}; "
                      f"dispersion ratio {json.dumps({d['groups'][g]: round(v, 3) for g, v in d['dispersion_ratio'].items()})})"]
        lines.append("")
    else:
        lines += ["## H5-J: fewshot8 not run yet.", ""]

    for t in DETECT_TARGETS:
        Tg, Ab = h2_tables(ctx, t)
        atomic_write_csv(os.path.join(rep, f"{t}_h2_tail.csv"), Tg)
        atomic_write_csv(os.path.join(rep, f"{t}_h2_ablation.csv"), Ab)
        chg = Ab.coverage_ablated - Ab.coverage
        decisions[t]["H2_descriptive"] = dict(median_coverage_change=float(chg.median()),
                                              share_groups_below_080_before=float((Ab.coverage < 0.8).mean()),
                                              share_groups_below_080_after=float((Ab.coverage_ablated < 0.8).mean()),
                                              median_width_ratio=float(Ab.width_ratio.median()))
        lines += [f"## {t}: H2 (descriptive): calibration tail and calibration-only ablation", "",
                  "Groups most over-represented in the calibration tail:", "", _md(Tg.head(8)), "",
                  f"Removing the {m8_['h2']['drop_share']:.0%} highest-MAE calibration groups: median change in held-out "
                  f"group coverage {chg.median():+.3f}; groups below 0.80 {(Ab.coverage < 0.8).mean():.1%} -> "
                  f"{(Ab.coverage_ablated < 0.8).mean():.1%}; median width ratio {Ab.width_ratio.median():.2f}.", ""]

    for t in ctx.m["targets"]:
        O, s = oxide_tables(ctx, t)
        atomic_write_csv(os.path.join(rep, f"{t}_oxides.csv"), O)
        decisions[t]["oxides_descriptive"] = s
        lines += [f"## {t}: oxide subset (descriptive; oxidation-state rule; groups with >= "
                  f"{m8_['oxide']['min_rows']} oxide structures)", "",
                  f"{s['n_groups']} groups; median oxide coverage {s['median_cov_oxide']:.3f}; below 0.90 "
                  f"{s['share_below_090']:.1%}; below 0.80 {s['share_below_080']:.1%}.", "", _md(O), ""]

    # H6: separate results side by side (never pooled)
    v7dec = os.path.join(ctx.d["report"], "phase6_decisions_v7.json")
    v7 = json.load(open(v7dec)) if os.path.exists(v7dec) else {}
    H6 = []
    for comp, nomad in m8_["h6_nomad"].items():
        r = dict(component=comp, NOMAD=nomad)
        for t in DETECT_TARGETS:
            d7, d8 = v7.get(t, {}), decisions.get(t, {})
            val = {"H1": d7.get("H1J", {}).get("verdict"),
                   "G5": ("holds" if d7.get("G5", {}).get("holds") else "fails") if d7.get("G5") else None,
                   "G3": ("positive" if d7.get("G3J", {}).get("positive") else "negative") if d7.get("G3J") else None,
                   "H3": d7.get("H3J", {}).get("verdict"),
                   "H4": (("supported" if d8["H4J"]["H4J_supported"] else "not supported") + "; G4 "
                          + ("triggered" if d8["H4J"]["G4J_triggered"] else "not triggered")) if "H4J" in d8 else None,
                   "H5": ("supported" if d8["H5J"]["H5J_supported"] else "not supported") if "H5J" in d8 else None}.get(comp)
            r[f"JARVIS {t}"] = val if val is not None else "n/a"
        H6.append(r)
    H6 = pd.DataFrame(H6)
    atomic_write_csv(os.path.join(rep, "h6_side_by_side.csv"), H6)
    lines += ["## H6: NOMAD and JARVIS side by side (separate analyses; never pooled)", "", H6.to_markdown(index=False), ""]

    decisions["_meta"] = dict(addendum_time_utc=pre["time_utc"], v7_rules_time_utc=pre["v7_preregistration"]["time_utc"],
                              bootstrap=B, n_perm=n_perm)
    atomic_write_json(os.path.join(rep, "phase6_decisions_v8.json"), json.loads(json.dumps(decisions, default=float)))
    text = "\n".join(lines)
    atomic_write_text(os.path.join(rep, "PHASE6_TASK65_SUMMARY.md"), text)
    print(text)
    print(f"\n[report8] written to {rep}")
