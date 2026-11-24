"""Phase 6 Task 6.5 (v8): computing steps. Everything is read from phase07_v1 (v7 units, designs,
tuned hyperparameters) and phase06_v1 (features, table); everything is written to phase08_v1.

  scores8   per Regime B outer fold (repeat 1, formation and band gap): failure-detection scores
            for every held-out structure, exactly the NOMAD Phase 4 set adapted to JARVIS
              knn         mean distance to the 25 nearest outer-training structures (standardized
                          features, PCA-50; NOMAD Phase 4 distance = the 'distance to training')
              kde         Schultz et al. (2025) KDE dissimilarity (reference set subsampled, an8)
              committee   SD of the 5 inner (space-group fold) models = MatFold's ensemble SD
              width       CV+ width at alpha 0.10 (the primary interval)
              classifier  out-of-fold domain-classifier probability (formula-grouped folds)
            and the conformal kNN flag (Bates et al. 2023) with an in-distribution false-flag
            check. For repeat 2 and for TBmBJ only the kNN distance is computed (group distance).
  fewshot8  H5-J: few-shot recalibration of the 3 worst-covered groups per target with the fixed
            repeat-1 model, NOMAD Phase 5 code (split and offset-corrected full conformal; random,
            clustered and one-per-composition draws; permutation null).
  hjplus8   hierarchical jackknife+ (Lee, Barber & Willett, eq. 20, Theorem 5): for every repeat-1
            outer fold, one model per training space group k fitted without k (rare groups always
            in training; one seed for all K models and rows in a fixed order, so the fitting
            algorithm is symmetric in the groups as Theorem 5 requires). Resumable per model.
"""
from __future__ import annotations

import json
import os
import time
import warnings

import numpy as np
import pandas as pd

from . import cvfast as cf
from .config import atomic_save_npz, atomic_write_csv, atomic_write_json, get_paths
from .design7 import B_rowsets, _seed
from .fit7 import MOCK, Ctx, _lgb_task, _metrics, _run, dirs7

DETECT_TARGETS = ("formation", "bandgap")


# ----------------------------------------------------------------------------- context
def m8(cfg) -> dict:
    m = json.loads(json.dumps(cfg["models8"]))
    if MOCK:
        m["out_dir"] = "phase08_mock"
        m["figures_dir"] = "paper_figures_mock8"
        m["bootstrap"] = 200
        m["n_perm"] = 2000
        m["fewshot"]["draws"] = 60
        m["kde"].update(n_reference=3000, n_bandwidth=1000)
    return m


def dirs8(cfg) -> dict:
    p = get_paths(cfg)
    root = os.path.join(p.drive_root, m8(cfg)["out_dir"])
    d = {k: os.path.join(root, k) for k in ("phase0", "scores", "fewshot", "hjplus", "report")}
    d["hjparts"] = os.path.join(d["hjplus"], "parts")
    for v in d.values():
        os.makedirs(v, exist_ok=True)
    d["root"] = root
    d["figures"] = os.path.join(p.drive_root, m8(cfg)["figures_dir"])
    return d


class Ctx8(Ctx):
    """v7 context (table, labels, features, designs, tuned hyperparameters) plus v8 folders."""

    def __init__(self, cfg):
        super().__init__(cfg)
        self.m8 = m8(cfg)
        self.d8 = dirs8(cfg)
        self.formula = self.F["reduced_formula"].astype(str).to_numpy()
        self.sgsym = dict(zip(self.L["sga01_number"].astype(int), self.L["sga01_symbol"].astype(str)))

    def unit(self, target, design, rep, o) -> dict:
        """A v7 Regime B unit, checked against the design (rows, groups)."""
        stem = os.path.join(self.d["unitsB"], f"{target}_{design}_r{rep}_o{o:02d}")
        if not (os.path.exists(stem + ".npz") and os.path.exists(stem + ".csv")):
            raise FileNotFoundError(f"v7 unit missing: {stem}.npz/.csv (run fit7 first)")
        with np.load(stem + ".npz") as z:
            U = {k: z[k] for k in z.files}
        D = self.design(target)
        R = B_rowsets(D, rep, o, {"rep2": "primary"}.get(design, design))
        rows = D["rows"]
        if not (np.array_equal(U["test_rows"], rows[R["test"]]) and np.array_equal(U["cal_rows"], rows[R["cal"]])):
            raise AssertionError(f"{stem}: rows differ from the v7 design")
        U["rare_rows"] = rows[R["rare"]]
        U["test_grp"] = D["grp"][R["test"]]
        U["cal_grp"] = D["grp"][R["cal"]]
        U["half"] = R["half"]
        U["n_folds"] = R["n_folds"]
        U["csv"] = pd.read_csv(stem + ".csv")
        return U


def _impute(Xfit, *others):
    med = np.nanmedian(Xfit, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    fix = lambda A: np.where(np.isnan(A), med, A)
    return [fix(Xfit)] + [fix(A) for A in others]


# ----------------------------------------------------------------------------- scores8
def score_unit(ctx: Ctx8, target: str, rep: int, o: int, full: bool) -> str:
    from sklearn.model_selection import StratifiedGroupKFold
    from .detect import classifier_oof, conformal_pvalues, knn_distance
    from .an8 import kde_dissimilarity_sub
    f = os.path.join(ctx.d8["scores"], f"{target}_r{rep}_o{o:02d}.npz")
    if os.path.exists(f):
        return "skipped"
    t0 = time.time()
    m = ctx.m8
    U = ctx.unit(target, "primary" if rep == 1 else "rep2", rep, o)
    tr = np.sort(np.concatenate([U["rare_rows"], U["cal_rows"]]))
    te = U["test_rows"]
    Xtr, Xte = _impute(np.asarray(ctx.X[tr], float), np.asarray(ctx.X[te], float))
    seed = _seed(m["seed"], target, "scores", rep, o) % (2 ** 31)
    knn_te, knn = knn_distance(Xtr, Xte, m["distance"]["k"], m["distance"]["pca_components"], seed)
    out = dict(test_rows=te, train_rows=tr, knn=knn_te, knn_explained=np.array(knn.explained_))
    msg = f"knn (PCA-50 explains {knn.explained_:.2f})"
    if full:
        kd = m["kde"]
        kde_te, h, nref, nbw = kde_dissimilarity_sub(Xtr, Xte, kd["bandwidth_quantile"], kd["n_reference"],
                                                     kd["n_bandwidth"], seed + 1)
        c = m["classifier"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            clf_te, _ = classifier_oof(Xtr, Xte, ctx.formula[tr], ctx.formula[te], c["pca_components"],
                                       c["logistic_C"], c["cv_folds"], seed + 2)
        # conformal kNN flag: fit / calibration / in-distribution test parts of the training rows,
        # folds grouped by reduced formula and stratified by space group (as NOMAD Phase 4)
        fl = m["flag"]
        sg_tr = ctx.L["sga01_number"].to_numpy()[tr]
        fold = np.full(len(tr), -1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")      # rare space groups have fewer members than folds
            for k, (_, v) in enumerate(StratifiedGroupKFold(fl["split_folds"], shuffle=True, random_state=seed + 3)
                                       .split(np.zeros(len(tr)), sg_tr, ctx.formula[tr])):
                fold[v] = k
        fit, cal, idt = fold >= 2, fold == 0, fold == 1
        s_te, det = knn_distance(Xtr[fit], Xte, m["distance"]["k"], m["distance"]["pca_components"], seed + 4)
        p_te = conformal_pvalues(det.sigma(Xtr[cal]), s_te)
        p_idt = conformal_pvalues(det.sigma(Xtr[cal]), det.sigma(Xtr[idt]))
        out.update(kde=kde_te, kde_bandwidth=np.array(h), kde_n_reference=np.array(nref),
                   kde_n_bandwidth=np.array(nbw), classifier=clf_te, flag_p_test=p_te, flag_p_idtest=p_idt,
                   flag_idtest_rows=tr[idt], flag_n_cal=np.array(int(cal.sum())))
        msg += f", kde (h {h:.1f}, reference {nref}), classifier, flag (cal {int(cal.sum())})"
    atomic_save_npz(f, **out)
    return f"{(time.time() - t0) / 60:.1f} min | train {len(tr)} | test {len(te)} | {msg}"


def run_scores8(cfg, targets=None, allow_deviation=False):
    from .setup8 import require_prereg8
    require_prereg8(cfg, allow_deviation)
    ctx = Ctx8(cfg)
    n_out = 2 * ctx.m["outer_folds_per_half"]
    targets = targets or ctx.m8["targets"]
    for t in targets:
        reps = [1, 2] if "rep2" in ctx.m["blocks"][t] else [1]
        for rep in reps:
            for o in range(n_out):
                full = rep == 1 and t in DETECT_TARGETS
                print(f"  [scores8] {t} r{rep} o{o}: {score_unit(ctx, t, rep, o, full)}", flush=True)
    print("[scores8] done. Next: fewshot8, hjplus8, report8")


# ----------------------------------------------------------------------------- fewshot8 (H5-J)
def primary_group_coverage(ctx: Ctx8, target: str) -> pd.DataFrame:
    """Primary coverage per group (CV+, alpha_primary), from the v7 unit CSVs, mean over repeats."""
    a = float(ctx.m["alpha_primary"])
    n_out = 2 * ctx.m["outer_folds_per_half"]
    frames = []
    for design, rep in (("primary", 1), ("rep2", 2)):
        if design == "rep2" and "rep2" not in ctx.m["blocks"][target]:
            continue
        for o in range(n_out):
            d = pd.read_csv(os.path.join(ctx.d["unitsB"], f"{target}_{design}_r{rep}_o{o:02d}.csv"))
            d = d[(d.group != "ALL") & (d.method == "cvplus") & np.isclose(d.alpha, a)]
            frames.append(d.assign(rep=rep, outer=o))
    T = pd.concat(frames, ignore_index=True)
    T["group"] = T.group.astype(int)
    P = T.pivot_table(index="group", columns="rep", values="coverage")
    out = pd.DataFrame(dict(group=P.index.to_numpy(), coverage=P.mean(axis=1).to_numpy(),
                            coverage_r1=P[1].to_numpy()))
    out["n_repeats"] = P.notna().sum(axis=1).to_numpy()
    return out


def worst_groups(ctx: Ctx8, target: str, n: int) -> pd.DataFrame:
    """Pre-registered rule (v7): the n groups with the lowest primary coverage (mean over repeats);
    ties broken by lower repeat-1 coverage, then by space-group number."""
    C = primary_group_coverage(ctx, target)
    return C.sort_values(["coverage", "coverage_r1", "group"]).head(n).reset_index(drop=True)


def family_from_unit(ctx: Ctx8, target: str, g: int):
    """Residuals z = y - mu (mu = mean of the 5 repeat-1 inner models, which never saw group g)."""
    from .fewshot import Family
    n_out = 2 * ctx.m["outer_folds_per_half"]
    for o in range(n_out):
        U = ctx.unit(target, "primary", 1, o)
        if g in set(U["test_grp"].tolist()):
            msk = U["test_grp"] == g
            rows = U["test_rows"][msk]
            mu = U["mu_te"].astype(float).mean(axis=0)[msk]
            y = U["y_test"][msk].astype(float)
            return Family(z=y - mu, groups=ctx.formula[rows], ids=rows, y=y), o
    raise AssertionError(f"group {g} not found in any repeat-1 unit of {target}")


def run_family8(fam, m8_, alpha, seed, label) -> pd.DataFrame:
    """NOMAD Phase 5 procedure (phase5.run_family) for one JARVIS family; same functions."""
    from .fewshot import draw_clustered, draw_one_per_group, draw_random, group_weighted
    from .phase5 import _evaluate
    f = m8_["fewshot"]
    D, mtg = int(f["draws"]), int(f["min_test_groups"])
    allidx = np.arange(fam.N)
    out = []
    for scheme, null in (("random", False), ("clustered", False), ("clustered", True),
                         ("one_per_group", False), ("one_per_group", True)):
        for m in f["m"]:
            rng = np.random.default_rng(_seed(seed, label, scheme, m, null) % (2 ** 32 - 1))
            for d in range(D):
                z = fam.z[rng.permutation(fam.N)] if null else fam.z
                new, tg = None, None
                if scheme == "random":
                    cal, rest = draw_random(fam, m, rng)
                else:
                    r = (draw_clustered if scheme == "clustered" else draw_one_per_group)(fam, m, rng, mtg)
                    if r is None:
                        break
                    cal, new, tg = r
                    rest = np.setdiff1d(allidx, cal, assume_unique=True)
                new_fn = ((lambda c, _t=new, _g=tg: group_weighted(c, _t, fam, _g)) if new is not None else None)
                for row in _evaluate(z[cal], (z[new] if new is not None else None), z[rest], alpha, new_fn):
                    out.append(dict(scheme=scheme, null=null, m=int(m), draw=d, n_test=len(rest),
                                    n_test_groups=(len(tg) if tg is not None else np.nan), **row))
    return pd.DataFrame(out)


def run_fewshot8(cfg, allow_deviation=False):
    from .setup8 import require_prereg8
    require_prereg8(cfg, allow_deviation)
    ctx = Ctx8(cfg)
    a = float(ctx.m["alpha_primary"])
    sel = []
    for t in DETECT_TARGETS:
        W = worst_groups(ctx, t, int(ctx.m8["fewshot"]["n_worst"]))
        for _, w in W.iterrows():
            g = int(w.group)
            fam, o = family_from_unit(ctx, t, g)
            sel.append(dict(target=t, group=g, symbol=ctx.sgsym.get(g, str(g)), primary_coverage=w.coverage,
                            coverage_r1=w.coverage_r1, outer_fold=o, N=fam.N, n_formulas=fam.G,
                            bias=float(fam.z.mean()), mae=float(np.abs(fam.z).mean()),
                            sd_y=float(np.std(fam.y, ddof=1))))
            fn = os.path.join(ctx.d8["fewshot"], f"{t}_g{g}.csv.gz")
            if os.path.exists(fn):
                print(f"  [fewshot8] {t} {ctx.sgsym.get(g, g)}: done earlier, skipped", flush=True)
                continue
            t0 = time.time()
            df = run_family8(fam, ctx.m8, a, ctx.m8["seed"], f"{t}|{g}")
            df.insert(0, "group", g)
            df.insert(0, "target", t)
            df.to_csv(fn + ".tmp", index=False, compression="gzip")
            os.replace(fn + ".tmp", fn)
            print(f"  [fewshot8] {t} {ctx.sgsym.get(g, g)} (coverage {w.coverage:.3f}): N={fam.N}, "
                  f"formulas={fam.G}, {len(df)} rows ({time.time() - t0:.0f} s)", flush=True)
    atomic_write_csv(os.path.join(ctx.d8["fewshot"], "worst_groups.csv"), pd.DataFrame(sel))
    print("[fewshot8] done.")


# ----------------------------------------------------------------------------- hjplus8
def hjplus_unit(ctx: Ctx8, target: str, o: int, cfg) -> str:
    from .models import plan_parallel
    from .an8 import hjplus_definition
    stem = os.path.join(ctx.d8["hjplus"], f"{target}_hj_o{o:02d}")
    if os.path.exists(stem + ".npz") and os.path.exists(stem + ".csv"):
        return "skipped"
    t0 = time.time()
    D = ctx.design(target)
    R = B_rowsets(D, 1, o, "primary")
    rows = D["rows"]
    te, cal, rare = rows[R["test"]], rows[R["cal"]], rows[R["rare"]]
    cal_grp = D["grp"][R["cal"]]
    groups = np.unique(cal_grp)
    K = groups.size
    y = ctx.y[target]
    hp_src = f"h{R['half']}"
    hp = ctx.hp(target, hp_src)
    seed = _seed(ctx.m8["seed"], target, "hjplus", o)       # ONE seed for all K models (symmetry)
    todo = []
    for g in groups:
        pf = os.path.join(ctx.d8["hjparts"], f"{target}_o{o:02d}_g{int(g)}.npz")
        if not os.path.exists(pf):
            todo.append((g, pf))
    w, _ = plan_parallel(cfg)
    batch = max(1, 2 * w)
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        tasks = []
        for g, _ in chunk:
            fit = np.sort(np.concatenate([rare, cal[cal_grp != g]]))   # fixed (sorted) order
            tasks.append(_lgb_task(ctx, fit, y[fit], {"oof": cal[cal_grp == g], "test": te}, hp, seed))
        for (g, pf), r in zip(chunk, _run(tasks, cfg)):
            atomic_save_npz(pf, oof=r["oof"], test=r["test"], best_iteration=np.array(r["best_iteration"]),
                            group=np.array(int(g)))
        done = min(i + batch, len(todo))
        if done == len(todo) or (i // batch) % 5 == 4:
            print(f"    [hjplus8] {target} o{o}: {done}/{len(todo)} leave-one-group-out models "
                  f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    mu_te = np.empty((K, len(te)))
    oof = np.empty(len(cal))
    for k, g in enumerate(groups):
        with np.load(os.path.join(ctx.d8["hjparts"], f"{target}_o{o:02d}_g{int(g)}.npz")) as z:
            mu_te[k] = z["test"]
            oof[cal_grp == g] = z["oof"]
    gidx = np.searchsorted(groups, cal_grp)
    Rabs = np.abs(y[cal] - oof)
    alphas = [float(a) for a in ctx.m8["hjplus"]["alphas"]]
    a_p = float(ctx.m["alpha_primary"])
    intervals = {("hjplus", a): cf.cvplus(mu_te, mu_te, None, Rabs, gidx, K, a, group=cal_grp) for a in alphas}
    # check the production routine against eq. (20) written literally, on random test points
    rng = np.random.default_rng(seed % (2 ** 31))
    cols = rng.choice(len(te), size=min(int(ctx.m8["hjplus"]["spot_check"]), len(te)), replace=False)
    blo, bhi = hjplus_definition(mu_te[:, cols], Rabs, gidx, a_p)
    lo, hi = intervals[("hjplus", a_p)]
    for u, v in ((lo[cols], blo), (hi[cols], bhi)):
        same_inf = np.array_equal(np.isinf(u), np.isinf(v)) and np.array_equal(np.sign(u[np.isinf(u)]), np.sign(v[np.isinf(v)]))
        fin = np.isfinite(v)
        if not same_inf or (fin.any() and np.max(np.abs(u[fin] - v[fin]) / (1 + np.abs(v[fin]))) > 1e-9):
            raise AssertionError(f"hierarchical jackknife+ differs from eq. (20) ({target} o{o})")
    meta = dict(target=target, design="hjplus", rep=1, outer=o, half=R["half"], hp=hp_src,
                n_train=int(len(rare) + len(cal)), n_cal=int(len(cal)), n_rare=int(len(rare)), n_folds=K,
                n_cal_groups=K)
    rowsM = _metrics(meta, y[te], mu_te.mean(axis=0), intervals, D["grp"][R["test"]], alphas)
    arrays = dict(test_rows=te, y_test=y[te], mu_te=mu_te.astype(np.float32), cal_rows=cal,
                  cal_group=cal_grp.astype(np.int32), groups=groups.astype(np.int32), oof=oof.astype(np.float32),
                  spot_cols=cols)
    for (meth, a), (l, h) in intervals.items():
        arrays[f"{meth}__a{a:.2f}__lo"], arrays[f"{meth}__a{a:.2f}__hi"] = l.astype(np.float32), h.astype(np.float32)
    atomic_save_npz(stem + ".npz", **arrays)
    atomic_write_csv(stem + ".csv", pd.DataFrame(rowsM))
    c = [r for r in rowsM if r["group"] == "ALL" and abs(r["alpha"] - a_p) < 1e-9][0]
    return (f"{(time.time() - t0) / 60:.1f} min | K = {K} leave-one-group-out models | test {len(te)} rows | "
            f"coverage {c['coverage']:.3f} (unbounded {c['n_unbounded']}) | eq. (20) check on {len(cols)} points passed")


def run_hjplus8(cfg, targets=None, allow_deviation=False):
    from .setup8 import require_prereg8
    require_prereg8(cfg, allow_deviation)
    ctx = Ctx8(cfg)
    n_out = 2 * ctx.m["outer_folds_per_half"]
    for t in targets or ctx.m8["hjplus"]["targets"]:
        for o in range(n_out):
            print(f"  [hjplus8] {t} o{o}: {hjplus_unit(ctx, t, o, cfg)}", flush=True)
    print("[hjplus8] done.")
