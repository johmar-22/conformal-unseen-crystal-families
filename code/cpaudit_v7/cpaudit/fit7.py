"""Phase 6 Task 6.4 (v7): LightGBM fits and conformal intervals on JARVIS-DFT, in resumable units.

Unit types (each saved atomically as <stem>.npz + <stem>.csv; finished units are skipped):
  tune   cross-half nested tuning h0 / h1 (Regime B, H3) and RA (Regime A); Optuna TPE, MAE,
         early stopping, inner folds by space group (h0, h1) or by reduced formula (RA)
  B      Regime B outer fold: inner models (mean; plus 5 % / 95 % quantile models in `primary`),
         kNN sigma per inner fold, intervals
           primary  : CV+, group-weighted CV+, normalized CV+, CQR-CV+ (alpha 0.10), ensemble SD
           rep2     : as primary without CQR, new outer folds within the same halves
           polyex   : as rep2 on the polymorph-excluded training set
           randinner: CV+ with random formula-grouped inner folds (in-distribution calibration)
           locso    : one crystal system held out, inner leave-one-crystal-system-out
  A      Regime A split: CV+ random, split, Mondrian per space group, normalized split, CQR split
  H3     H3 partition: one fit on the training groups, HCP and single subsampling on the pool

The fast CV+ routine is spot-checked against the brute-force definition inside every unit.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import time

import numpy as np
import pandas as pd

from . import conformal as cp
from . import cvfast as cf
from .config import atomic_save_npz, atomic_write_csv, atomic_write_json, get_paths
from .design7 import (TARGETS, A_rowsets, B_rowsets, H3_partition, LOCSO_rowsets, _seed, build_target_design,
                      check_design, content_hash, load_design, save_design, scope_mask, tune_rowsets)

MOCK = os.environ.get("CPAUDIT_MOCK7", "") == "1"     # verification only: synthetic targets, tiny models


# ----------------------------------------------------------------------------- context
def m7(cfg):
    m = dict(cfg["models7"])
    if MOCK:
        m = json.loads(json.dumps(m))
        m["out_dir"] = "phase07_mock"
        m["tuning"].update(n_trials=2, n_startup_trials=1, max_rounds=30, early_stopping_rounds=5)
        m["h3"]["n_reshuffles"] = 20
    return m


def dirs7(cfg) -> dict:
    p = get_paths(cfg)
    root = os.path.join(p.drive_root, m7(cfg)["out_dir"])
    d = {k: os.path.join(root, k) for k in ("phase0", "design", "tune", "unitsB", "unitsA", "unitsH3", "report")}
    for v in d.values():
        os.makedirs(v, exist_ok=True)
    d["root"] = root
    d["phase6"] = os.path.join(p.drive_root, cfg["paths"]["out_dir"])
    d["cache"] = p.cache
    return d


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def content_hashes(d) -> dict:
    """SHA-256 of contents (decompressed table bytes, raw array bytes), not of compressed files,
    because gzip / npz bytes depend on the zlib build (classic zlib in Colab, zlib-ng elsewhere)."""
    p6 = d["phase6"]
    out = {}
    with gzip.open(os.path.join(p6, "data", "jarvis_frozen.csv.gz"), "rb") as f:
        out["frozen_table"] = sha_bytes(f.read())
    with gzip.open(os.path.join(p6, "splits", "labels.csv.gz"), "rb") as f:
        out["labels"] = sha_bytes(f.read())
    with np.load(os.path.join(p6, "features", "features_v1.npz")) as Z:
        out["features_X"] = sha_bytes(np.ascontiguousarray(Z["X"]).tobytes())
        out["features_jid"] = sha_bytes("\n".join(Z["jid"].tolist()).encode())
    with np.load(os.path.join(p6, "splits", "holdouts_v2.npz")) as Z:
        h = hashlib.sha256()
        for k in sorted(Z.files):
            h.update(k.encode())
            h.update(np.ascontiguousarray(Z[k]).tobytes())
        out["holdouts_v2"] = h.hexdigest()
    return out


class Ctx:
    """Loaded inputs: frozen table, labels, features (cached as .npy for memory-mapped workers)."""

    def __init__(self, cfg):
        from .data6 import read_table
        self.cfg, self.m = cfg, m7(cfg)
        self.d = dirs7(cfg)
        p6 = self.d["phase6"]
        self.F = read_table(os.path.join(p6, "data", "jarvis_frozen.csv.gz"))
        self.L = pd.read_csv(os.path.join(p6, "splits", "labels.csv.gz"), keep_default_na=False, na_values=[""])
        if not np.array_equal(self.L.jid.astype(str).to_numpy(), self.F.jid.astype(str).to_numpy()):
            raise AssertionError("labels.csv.gz rows differ from the frozen table")
        self.X_path = os.path.join(self.d["cache"], "jarvis_X_v1.npy")
        if not os.path.exists(self.X_path):
            with np.load(os.path.join(p6, "features", "features_v1.npz")) as Z:
                if not np.array_equal(Z["jid"].astype(str), self.F.jid.astype(str).to_numpy()):
                    raise AssertionError("feature rows differ from the frozen table")
                np.save(self.X_path + ".tmp.npy", np.ascontiguousarray(Z["X"], dtype=np.float64))
            os.replace(self.X_path + ".tmp.npy", self.X_path)
        self.X = np.load(self.X_path, mmap_mode="r")
        self.y = {}
        for t in TARGETS:
            y = pd.to_numeric(self.F[TARGETS[t]["col"]], errors="coerce").to_numpy(float)
            if MOCK:
                y = np.where(scope_mask(self.F, t), _mock_target(self.X, t), np.nan)
            self.y[t] = y
        self._design = {}

    def design(self, target):
        if target not in self._design:
            f = os.path.join(self.d["design"], f"design_{target}.npz")
            if not os.path.exists(f):
                raise FileNotFoundError(f"{f} not found: run setup7 first")
            self._design[target] = load_design(f)
        return self._design[target]

    def hp(self, target, which):
        f = os.path.join(self.d["tune"], f"hp_{target}_{which}.json")
        if not os.path.exists(f):
            raise FileNotFoundError(f"{f} not found: run tune7 first")
        with open(f) as fh:
            hp = json.load(fh)
        if MOCK:
            hp["n_rounds"] = min(int(hp["n_rounds"]), 30)
        return hp


def _mock_target(X, target):
    """Synthetic target for the verification run only (never used in the real analysis)."""
    rng = np.random.default_rng({"formation": 1, "bandgap": 2, "tbmbj": 3}[target])
    Z = np.nan_to_num((X - np.nanmean(X, 0)) / (np.nanstd(X, 0) + 1e-9))
    w = rng.normal(size=X.shape[1]) * (rng.random(X.shape[1]) < 0.1)
    return np.tanh(Z @ w / 4) + 0.2 * np.sin(Z[:, 3] * Z[:, 30]) + 0.1 * rng.normal(size=len(X))


# ----------------------------------------------------------------------------- design (setup7)
def build_designs(cfg, ctx=None) -> dict:
    ctx = ctx or Ctx(cfg)
    m = ctx.m
    out = {}
    for t in m["targets"]:
        f = os.path.join(ctx.d["design"], f"design_{t}.npz")
        if os.path.exists(f):
            D = load_design(f)
            print(f"[design7] {t}: reusing {f}", flush=True)
        else:
            t0 = time.time()
            D = build_target_design(ctx.F, ctx.L, t, m)
            save_design(f, D)
            D = load_design(f)
            print(f"[design7] {t}: built in {(time.time() - t0) / 60:.1f} min", flush=True)
        summ = check_design(D, m)
        if TARGETS[t]["scope"] == "all":
            summ["matches_holdouts_v2"] = _check_against_holdouts(ctx, D)
        summ["content_hash"] = content_hash(D)
        out[t] = summ
        print(f"[design7] {t}: {summ['n_rows']} rows, {summ['n_eligible']} eligible groups "
              f"({summ['rows_in_eligible']} rows), halves {summ['half_sizes_groups']} groups; "
              f"groups per outer fold {summ['outer_groups_per_fold_r1']}; checks passed", flush=True)
    return out


def _check_against_holdouts(ctx, D) -> bool:
    """For the all-chemistries scope, the eligible groups and their rows must be exactly the
    MatFold-label LOSGO test sets of v6.1 (which matfold6 matched to MatFold itself)."""
    with np.load(os.path.join(ctx.d["phase6"], "splits", "holdouts_v2.npz")) as Z:
        keys = [k for k in Z.files if k.startswith("losgo_mf|") and k.endswith("|test")]
        groups = sorted(int(k.split("|")[1]) for k in keys)
        if groups != sorted(int(g) for g in D["eligible"]):
            raise AssertionError("eligible groups differ from the v6.1 MatFold holdouts")
        for g in groups:
            if not np.array_equal(np.sort(Z[f"losgo_mf|{g}|test"]), np.sort(D["rows"][D["grp"] == g])):
                raise AssertionError(f"rows of space group {g} differ from the v6.1 MatFold holdout")
    return True


# ----------------------------------------------------------------------------- fitting helpers
def _lgb_task(ctx, rows_fit, y_fit, predict: dict, hp, seed, quantile=None):
    from .models import lgb_params, plan_parallel
    _, threads = plan_parallel(ctx.cfg)
    return dict(X_path=ctx.X_path, max_bin=ctx.cfg["lgbm"]["max_bin"], n_rounds=int(hp["n_rounds"]),
                rows_train=np.asarray(rows_fit), y_train=np.asarray(y_fit),
                params=lgb_params(hp, ctx.cfg, int(seed) % (2 ** 31), threads, quantile=quantile),
                predict={k: np.asarray(v) for k, v in predict.items()})


def _run(tasks, cfg):
    from .models import plan_parallel, run_tasks
    workers, _ = plan_parallel(cfg)
    return run_tasks(tasks, workers)


def _sigma(ctx, rows_fit, queries: dict, seed):
    """kNN difficulty (as on NOMAD: mean distance to the 25 nearest fit rows in PCA-50 of the
    standardized features), with NaN replaced by the fit rows' column medians."""
    from .difficulty import KNNDifficulty
    nz = ctx.m["normalized"]
    Xf = np.asarray(ctx.X[np.sort(rows_fit)], float)
    med = np.nanmedian(Xf, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    Xf = np.where(np.isnan(Xf), med, Xf)
    knn = KNNDifficulty(nz["knn_k"], nz["pca_components"], nz["beta"], int(seed) % (2 ** 31)).fit(Xf)
    out = {}
    for k, r in queries.items():
        Q = np.asarray(ctx.X[np.asarray(r)], float)
        out[k] = knn.sigma(np.where(np.isnan(Q), med, Q))
    return out, knn.explained_


def _metrics(rows_meta, y, mubar, intervals, groups, alphas_of):
    """One row per (group, method, alpha)."""
    out = []
    ug = np.unique(groups)
    for g in list(ug) + ["ALL"]:
        m = np.ones(len(y), bool) if g == "ALL" else groups == g
        yy, mm = y[m], mubar[m]
        sd = float(np.std(yy, ddof=1)) if m.sum() > 1 else np.nan
        base = dict(rows_meta, group=str(g), n=int(m.sum()), mae=float(np.mean(np.abs(yy - mm))),
                    bias=float(np.mean(yy - mm)), sd_y=sd)
        for (meth, a), (lo, hi) in intervals.items():
            l, h = lo[m], hi[m]
            w = h - l
            cov = (yy >= l) & (yy <= h)
            fin = np.isfinite(w)
            out.append(dict(base, method=meth, alpha=a, coverage=float(cov.mean()),
                            miss_below=float((yy < l).mean()), miss_above=float((yy > h).mean()),
                            width_median=float(np.median(w)), width_mean=float(np.mean(w)) if fin.all() else np.inf,
                            width_median_sd=float(np.median(w) / sd) if sd and np.isfinite(sd) else np.nan,
                            interval_score=float(np.mean(cf.interval_score(yy, l, h, a))) if fin.all() else np.inf,
                            n_unbounded=int((~fin).sum())))
    return out


def _done(stem):
    return os.path.exists(stem + ".npz") and os.path.exists(stem + ".csv")


def _f32(a):
    return np.asarray(a, np.float32)


# ----------------------------------------------------------------------------- tuning
def run_tune(cfg, ctx, target, which):
    from .models import tune
    D = ctx.design(target)
    pos, folds = tune_rowsets(D, which)
    rows = D["rows"][pos]
    out_json = os.path.join(ctx.d["tune"], f"hp_{target}_{which}.json")
    t0 = time.time()
    hp = tune(ctx.X_path, rows, ctx.y[target][rows], folds, ctx.m["tuning"], cfg,
              seed=_seed(ctx.m["seed"], target, "tune", which) % (2 ** 31), out_json=out_json,
              log_json=os.path.join(ctx.d["tune"], f"trials_{target}_{which}.json"), label=f"{target} {which}")
    print(f"[tune7] {target} {which}: n={len(rows)}, {len(folds)} folds, n_rounds {hp['n_rounds']}, "
          f"CV MAE {hp['_cv_mae']:.5f} ({(time.time() - t0) / 60:.1f} min)", flush=True)
    return hp


# ----------------------------------------------------------------------------- Regime B unit
def run_B(cfg, ctx, target, design, rep, o):
    stem = os.path.join(ctx.d["unitsB"], f"{target}_{design}_r{rep}_o{o:02d}")
    if _done(stem):
        return "skipped"
    t0 = time.time()
    m = ctx.m
    D = ctx.design(target)
    alphas = [float(a) for a in m["alphas"]]
    a_p = float(m["alpha_primary"])
    if design == "locso":
        R = LOCSO_rowsets(D, o)
        hp_src = m["locso_hp"]
        grp_key = D["crys"]
    else:
        R = B_rowsets(D, rep, o, {"rep2": "primary"}.get(design, design))
        hp_src = f"h{R['half']}"
        grp_key = D["grp"]
    if len(R["test"]) == 0:
        return "empty"
    hp = ctx.hp(target, hp_src)
    rows, y = D["rows"], ctx.y[target]
    te, cal, rare = rows[R["test"]], rows[R["cal"]], rows[R["rare"]]
    f, K = R["cal_fold"], R["n_folds"]
    y_te, y_cal = y[te], y[cal]
    quant = design == "primary" and target != "tbmbj"
    seed = _seed(m["seed"], target, design, rep, o)
    tasks, names = [], []
    for k in range(K):
        fit = np.concatenate([rare, cal[f != k]])
        pred = {"oof": cal[f == k], "test": te}
        tasks.append(_lgb_task(ctx, fit, y[fit], pred, hp, seed + 11 * k))
        names.append(("mu", k))
        if quant:
            for nm, q in zip(("qlo", "qhi"), m["cqr_quantiles"]):
                tasks.append(_lgb_task(ctx, fit, y[fit], pred, hp, seed + 11 * k, quantile=float(q)))
                names.append((nm, k))
    res = dict(zip(names, _run(tasks, cfg)))

    def stack(nm):
        te_ = np.stack([res[(nm, k)]["test"] for k in range(K)])
        oof = np.empty(len(cal))
        for k in range(K):
            oof[f == k] = res[(nm, k)]["oof"]
        return te_, oof

    mu_te, mu_oof = stack("mu")
    Rabs = np.abs(y_cal - mu_oof)
    mubar = mu_te.mean(axis=0)
    cal_grp = grp_key[R["cal"]]
    intervals = {}
    for a in alphas:
        intervals[("cvplus", a)] = cf.cvplus(mu_te, mu_te, None, Rabs, f, K, a)
    if design != "randinner":
        sig_te = np.empty((K, len(te)))
        sig_oof = np.empty(len(cal))
        expl = []
        for k in range(K):
            fit = np.concatenate([rare, cal[f != k]])
            s, e = _sigma(ctx, fit, {"oof": cal[f == k], "test": te}, seed + 7 * k)
            sig_oof[f == k], sig_te[k] = s["oof"], s["test"]
            expl.append(e)
        Snorm = Rabs / sig_oof
        sd_ens = mu_te.std(axis=0, ddof=1)
        from scipy.stats import norm
        for a in alphas:
            if design != "locso":
                intervals[("cvplus_gw", a)] = cf.cvplus(mu_te, mu_te, None, Rabs, f, K, a, group=cal_grp)
            intervals[("ncvplus", a)] = cf.cvplus(mu_te, mu_te, sig_te, Snorm, f, K, a)
            z = float(norm.ppf(1 - a / 2))
            intervals[("enssd", a)] = (mubar - z * sd_ens, mubar + z * sd_ens)
    if quant:
        ql_te, ql_oof = stack("qlo")
        qh_te, qh_oof = stack("qhi")
        E = np.maximum(ql_oof - y_cal, y_cal - qh_oof)
        intervals[("cqr_cvplus", a_p)] = cf.cvplus(ql_te, qh_te, None, E, f, K, a_p)
    # spot checks of the fast CV+ against the definition
    sc = m["spot_check"]
    cf.spot_check(mu_te, mu_te, None, Rabs, f, K, a_p, *intervals[("cvplus", a_p)], n_check=sc["n_unweighted"], seed=seed)
    if ("cvplus_gw", a_p) in intervals:
        cf.spot_check(mu_te, mu_te, None, Rabs, f, K, a_p, *intervals[("cvplus_gw", a_p)], group=cal_grp,
                      n_check=sc["n_weighted"], seed=seed)
    if ("ncvplus", a_p) in intervals:
        cf.spot_check(mu_te, mu_te, sig_te, Snorm, f, K, a_p, *intervals[("ncvplus", a_p)], n_check=sc["n_unweighted"], seed=seed)
    if quant:
        cf.spot_check(ql_te, qh_te, None, E, f, K, a_p, *intervals[("cqr_cvplus", a_p)], n_check=sc["n_unweighted"], seed=seed)

    meta = dict(target=target, design=design, rep=rep, outer=o, half=R["half"], hp=hp_src,
                n_train=int(len(rare) + len(cal)), n_cal=int(len(cal)), n_rare=int(len(rare)), n_folds=K,
                n_cal_groups=int(len(np.unique(cal_grp))))
    rowsM = _metrics(meta, y_te, mubar, intervals, grp_key[R["test"]], alphas)
    arrays = dict(test_rows=te, y_test=y_te, mu_te=_f32(mu_te), cal_rows=cal, cal_fold=f.astype(np.int8),
                  cal_group=np.asarray(cal_grp).astype(str) if design == "locso" else cal_grp.astype(np.int32),
                  oof=_f32(mu_oof), test_groups=np.asarray(R["test_groups"]).astype(str),
                  best_iters=np.array([res[n]["best_iteration"] for n in names]))
    if design != "randinner":
        arrays.update(sig_te=_f32(sig_te), sig_oof=_f32(sig_oof), pca_explained=np.array(expl))
    if quant:
        arrays.update(ql_te=_f32(ql_te), qh_te=_f32(qh_te), E=_f32(E))
    for (meth, a), (lo, hi) in intervals.items():
        arrays[f"{meth}__a{a:.2f}__lo"], arrays[f"{meth}__a{a:.2f}__hi"] = _f32(lo), _f32(hi)
    atomic_save_npz(stem + ".npz", **arrays)
    atomic_write_json(stem + ".hp.json", hp)
    atomic_write_csv(stem + ".csv", pd.DataFrame(rowsM))
    c = [r for r in rowsM if r["group"] == "ALL" and r["method"] == "cvplus" and abs(r["alpha"] - a_p) < 1e-9][0]
    return (f"{(time.time() - t0) / 60:.1f} min | test {len(te)} rows in {len(R['test_groups'])} groups | "
            f"MAE {c['mae']:.4f} | CV+ coverage (all held-out rows) {c['coverage']:.3f}")


# ----------------------------------------------------------------------------- Regime A unit
def run_A(cfg, ctx, target, s):
    stem = os.path.join(ctx.d["unitsA"], f"{target}_A_s{s:02d}")
    if _done(stem):
        return "skipped"
    t0 = time.time()
    m = ctx.m
    D = ctx.design(target)
    alphas = [float(a) for a in m["alphas"]]
    a_p = float(m["alpha_primary"])
    R = A_rowsets(D, s)
    hp = ctx.hp(target, "RA")
    rows, y = D["rows"], ctx.y[target]
    te, pool, prop, cal = rows[R["test"]], rows[R["pool"]], rows[R["proper"]], rows[R["cal"]]
    fold = R["cv_fold"]
    seed = _seed(m["seed"], target, "A", s)
    tasks, names = [], []
    for k in range(5):
        fit = pool[fold != k]
        tasks.append(_lgb_task(ctx, fit, y[fit], {"oof": pool[fold == k], "test": te}, hp, seed + 11 * k))
        names.append(("cv", k))
    tasks.append(_lgb_task(ctx, prop, y[prop], {"cal": cal, "test": te}, hp, seed + 101))
    names.append(("mu", 0))
    for nm, q in zip(("qlo", "qhi"), m["cqr_quantiles"]):
        tasks.append(_lgb_task(ctx, prop, y[prop], {"cal": cal, "test": te}, hp, seed + 101, quantile=float(q)))
        names.append((nm, 0))
    res = dict(zip(names, _run(tasks, cfg)))
    mu_cv_te = np.stack([res[("cv", k)]["test"] for k in range(5)])
    oof = np.empty(len(pool))
    for k in range(5):
        oof[fold == k] = res[("cv", k)]["oof"]
    Rcv = np.abs(y[pool] - oof)
    mu_cal, mu_te = res[("mu", 0)]["cal"], res[("mu", 0)]["test"]
    ql_cal, ql_te = res[("qlo", 0)]["cal"], res[("qlo", 0)]["test"]
    qh_cal, qh_te = res[("qhi", 0)]["cal"], res[("qhi", 0)]["test"]
    s_, expl = _sigma(ctx, prop, {"cal": cal, "test": te}, seed + 7)
    elig = D["eligible"]
    cat = lambda pos: np.where(np.isin(D["grp"][pos], elig), D["grp"][pos], -1)
    cat_cal, cat_te = cat(R["cal"]), cat(R["test"])
    y_cal, y_te = y[cal], y[te]
    intervals = {}
    for a in alphas:
        intervals[("cvplus_random", a)] = cf.cvplus(mu_cv_te, mu_cv_te, None, Rcv, fold, 5, a)
        lo, hi, _ = cp.split_interval(mu_te, np.abs(y_cal - mu_cal), a)
        intervals[("split", a)] = (lo, hi)
        lo, hi, _ = cp.mondrian_interval(mu_te, cat_te, np.abs(y_cal - mu_cal), cat_cal, a)
        intervals[("mondrian", a)] = (lo, hi)
        lo, hi, _ = cp.normalized_interval(mu_te, s_["test"], y_cal, mu_cal, s_["cal"], a)
        intervals[("nsplit", a)] = (lo, hi)
    lo, hi, _ = cp.cqr_interval(ql_te, qh_te, y_cal, ql_cal, qh_cal, a_p)
    intervals[("cqr_split", a_p)] = (lo, hi)
    cf.spot_check(mu_cv_te, mu_cv_te, None, Rcv, fold, 5, a_p, *intervals[("cvplus_random", a_p)],
                  n_check=m["spot_check"]["n_unweighted"], seed=seed)
    meta = dict(target=target, design="regimeA", rep=s // 5, outer=s, half=-1, hp="RA",
                n_train=int(len(pool)), n_cal=int(len(cal)), n_rare=0, n_folds=5, n_cal_groups=-1)
    # point prediction for MAE: the CV+ ensemble for CV+, the proper model otherwise; report the proper model
    rowsM = _metrics(meta, y_te, mu_te, intervals, cat_te, alphas)
    mae_cv = float(np.mean(np.abs(y_te - mu_cv_te.mean(axis=0))))
    for r in rowsM:
        r["mae_cv_ensemble"] = mae_cv if r["group"] == "ALL" else np.nan
    arrays = dict(test_rows=te, y_test=y_te, mu_cv_te=_f32(mu_cv_te), pool_rows=pool, cv_fold=fold.astype(np.int8),
                  oof=_f32(oof), cal_rows=cal, mu_cal=_f32(mu_cal), mu_te=_f32(mu_te), ql_cal=_f32(ql_cal),
                  qh_cal=_f32(qh_cal), ql_te=_f32(ql_te), qh_te=_f32(qh_te), sig_cal=_f32(s_["cal"]),
                  sig_te=_f32(s_["test"]), pca_explained=np.array([expl]),
                  best_iters=np.array([res[n]["best_iteration"] for n in names]))
    for (meth, a), (lo, hi) in intervals.items():
        arrays[f"{meth}__a{a:.2f}__lo"], arrays[f"{meth}__a{a:.2f}__hi"] = _f32(lo), _f32(hi)
    atomic_save_npz(stem + ".npz", **arrays)
    atomic_write_csv(stem + ".csv", pd.DataFrame(rowsM))
    c = [r for r in rowsM if r["group"] == "ALL" and r["method"] == "cvplus_random" and abs(r["alpha"] - a_p) < 1e-9][0]
    return f"{(time.time() - t0) / 60:.1f} min | test {len(te)} | CV+ random coverage {c['coverage']:.3f}"


# ----------------------------------------------------------------------------- H3 unit
def run_H3(cfg, ctx, target, p):
    stem = os.path.join(ctx.d["unitsH3"], f"{target}_H3_p{p:02d}")
    if _done(stem):
        return "skipped"
    t0 = time.time()
    m = ctx.m
    D = ctx.design(target)
    P = H3_partition(D, p, m)
    hp = ctx.hp(target, f"h{P['half']}")
    rows, y = D["rows"], ctx.y[target]
    tr, pool = rows[P["train"]], rows[P["pool"]]
    res = _run([_lgb_task(ctx, tr, y[tr], {"pool": pool}, hp, _seed(m["seed"], target, "H3fit", p))], cfg)[0]
    mu = res["pool"]
    score = np.abs(y[pool] - mu)
    g = D["grp"][P["pool"]]
    rng = np.random.default_rng(P["seed"])
    out = []
    for r in range(int(m["h3"]["n_reshuffles"])):
        cal_groups = rng.choice(P["pool_groups"], P["n_cal"], replace=False)
        is_cal = np.isin(g, cal_groups)
        test_groups = [x for x in P["pool_groups"] if x not in set(cal_groups)]
        for a in [float(x) for x in m["alphas"]]:
            T = {"hcp": cf.hcp_threshold(score[is_cal], g[is_cal], a),
                 "single_subsample": cf.single_subsample_threshold(score[is_cal], g[is_cal], a, rng)}
            for meth, thr in T.items():
                for tg in test_groups:
                    mm = g == tg
                    out.append(dict(target=target, partition=p, half=P["half"], reshuffle=r, method=meth, alpha=a,
                                    test_group=int(tg), n=int(mm.sum()), threshold=thr,
                                    coverage=float(np.mean(score[mm] <= thr))))
    atomic_save_npz(stem + ".npz", pool_rows=pool, mu=_f32(mu), score=_f32(score), pool_group=g.astype(np.int32),
                    best_iteration=np.array([res["best_iteration"]]))
    df = pd.DataFrame(out)
    atomic_write_csv(stem + ".csv", df)
    gw = df[(df.method == "hcp") & (df.alpha.round(3) == round(float(m["alpha_primary"]), 3))]
    return (f"{(time.time() - t0) / 60:.1f} min | pool {len(P['pool_groups'])} groups | HCP group-weighted "
            f"coverage {gw.groupby('reshuffle').coverage.mean().mean():.3f}")


# ----------------------------------------------------------------------------- orchestration
BLOCK_UNITS = {
    "primary": lambda D, m: [("B", "primary", 1, o) for o in range(2 * m["outer_folds_per_half"])],
    "randinner": lambda D, m: [("B", "randinner", 1, o) for o in range(2 * m["outer_folds_per_half"])],
    "polyex": lambda D, m: [("B", "polyex", 1, o) for o in range(2 * m["outer_folds_per_half"])],
    "rep2": lambda D, m: [("B", "rep2", 2, o) for o in range(2 * m["outer_folds_per_half"])],
    "locso": lambda D, m: [("B", "locso", 0, o) for o in range(len(D["crystal_systems"]))],
    "regimeA": lambda D, m: [("A", None, None, s) for s in range(10)],
    "h3": lambda D, m: [("H3", None, None, p) for p in range(int(m["h3"]["n_partitions"]))],
}
TUNE_RUNS = lambda blocks: ["h0", "h1"] + (["RA"] if "regimeA" in blocks else [])


def run_tune7(cfg, targets=None, allow_deviation=False):
    from .setup7 import require_preregistration
    require_preregistration(cfg, allow_deviation)
    ctx = Ctx(cfg)
    for t in targets or ctx.m["targets"]:
        for which in TUNE_RUNS(ctx.m["blocks"][t]):
            run_tune(cfg, ctx, t, which)


def run_fit7(cfg, targets=None, blocks=None, allow_deviation=False):
    from .setup7 import require_preregistration
    require_preregistration(cfg, allow_deviation)
    ctx = Ctx(cfg)
    m = ctx.m
    order = ["primary", "regimeA", "randinner", "polyex", "h3", "rep2", "locso"]
    for b in order:
        for t in targets or m["targets"]:
            if b not in m["blocks"][t] or (blocks and b not in blocks):
                continue
            D = ctx.design(t)
            units = BLOCK_UNITS[b](D, m)
            print(f"[fit7] {t} {b}: {len(units)} units", flush=True)
            for kind, design, rep, k in units:
                if kind == "B":
                    msg = run_B(cfg, ctx, t, design, rep, k)
                elif kind == "A":
                    msg = run_A(cfg, ctx, t, k)
                else:
                    msg = run_H3(cfg, ctx, t, k)
                print(f"  [fit7] {t} {b} unit {k}: {msg}", flush=True)
