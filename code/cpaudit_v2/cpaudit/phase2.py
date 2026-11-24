"""Phase 2: the leave-one-lattice-out (LOSO) mechanism, Tasks 2.1-2.5 of the revision plan.

For each target and held-out lattice g (the test set is every entry of g):

  HPO       Optuna TPE with leave-one-lattice-out inner folds over the five training lattices
            (nested: g is never used for tuning).
  Fits      5 lattice-fold mean models, 10 lattice-fold quantile models (CQR), 5 composition-grouped
            random-fold mean models, all on the five training lattices only.
  Methods   cvplus_lattice      CV+ (Barber et al. 2021) with the lattice folds
            ncvplus_knn         nested CV+ with |y - mu| / sigma scores, sigma = kNN distance fitted
                                on each fold's training lattices (Gupta et al. 2022; Hu et al. 2022)
            cqr_cvplus          nested CV+ with CQR scores (Romano et al. 2019; Gupta et al. 2022)
            cvplus_random_in    CV+ with random composition-grouped folds: in-distribution residuals
            *_noFd3m_cal        calibration-only ablation: the Fd-3m fold's terms removed (g != Fd-3m)
            group_level         Dunn et al. (2023) family-level threshold (illustrative, Task 2.4)
  Records   every out-of-fold residual with its lattice (Task 2.1 calibration composition).

Leave-two-out (Task 2.2): for X != Fd-3m, the same procedure on the four lattices other than
X and Fd-3m (4 lattice folds, hyperparameters of LOSO-X), tested on X.

All outputs go to <out_dir>/loso and <out_dir>/hpo and are written atomically; finished units
are skipped on restart and HPO resumes from its per-trial log.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

from . import conformal as cp
from .config import atomic_save_npz, atomic_write_csv, atomic_write_json, get_paths
from .data import LATTICE_NAME, SUBFAMILY_NAME, build_dataset, ensure_feature_cache, load_features
from .difficulty import KNNDifficulty
from .metrics import summarize
from .models import lgb_params, plan_parallel, run_tasks, tune
from .setup2 import check_preregistration

REF = 227


def _seed(cfg, target, g, extra=0):
    return int(cfg["seed_base"]) + 100 * int(g) + (0 if target == "bandgap" else 50) + extra


def _sigma_per_fold(X, rows_tr, fold, rows_te, cfg, seed):
    """sigma_{-k}(x): kNN distance to fold k's training lattices; for OOF points and test points."""
    nz = cfg["loso"]["normalized"]
    K = int(fold.max()) + 1
    sig_oof = np.empty(len(rows_tr))
    sig_te = np.empty((K, len(rows_te)))
    Xte = np.asarray(X[rows_te], np.float64)
    for k in range(K):
        knn = KNNDifficulty(nz["knn_k"], nz["pca_components"], nz["beta"], seed + k).fit(
            np.asarray(X[rows_tr[fold != k]], np.float64))
        sig_oof[fold == k] = knn.sigma(np.asarray(X[rows_tr[fold == k]], np.float64))
        sig_te[k] = knn.sigma(Xte)
    return sig_oof, sig_te


def _metric_rows(unit, y_te, masks, intervals, preds):
    rows = []
    for (method, a), (lo, hi) in intervals.items():
        for gname, m in masks:
            if m.sum() == 0:
                continue
            r = dict(unit, method=method, alpha=a, group=gname)
            r.update(summarize(y_te[m], preds[method][m], lo[m], hi[m], a))
            rows.append(r)
    return rows


def _masks(ds, te, g):
    conflict = ds.conflict[te]
    masks = [("heldout", np.ones(len(te), bool)), ("heldout_noconflict", ~conflict)]
    if g == REF:
        masks += [(SUBFAMILY_NAME[s], ds.sub[te] == s) for s in ("227a", "227b")]
    return masks


def _calibration_composition(R, lat, alphas, label):
    rows = []
    n = R.size
    for a in alphas:
        q = cp.conformal_quantile(R, a)
        above = R >= q
        for c in np.unique(lat):
            m = lat == c
            share_pts = m.sum() / n
            share_above = (above & m).sum() / max(above.sum(), 1)
            rows.append(dict(folds=label, alpha=a, cal_lattice=LATTICE_NAME[c], n_cal=int(m.sum()),
                             share_points=share_pts, n_at_or_above_q=int((above & m).sum()),
                             share_at_or_above_q=share_above, ratio=share_above / share_pts,
                             pooled_q=q, lattice_p90=float(np.quantile(R[m], 0.9)),
                             lattice_mean_resid=float(R[m].mean())))
    return rows


# ----------------------------------------------------------------------------- LOSO unit
def run_loso_unit(ds, target, g, cfg, paths, X_path, X):
    name = LATTICE_NAME[g]
    stem = os.path.join(paths.loso, f"{target}_loso{g}")
    if os.path.exists(stem + ".csv") and os.path.exists(stem + ".npz"):
        return "skipped"
    t0 = time.time()
    L = cfg["loso"]
    alphas = [float(a) for a in cfg["alpha"]]
    a_p = float(cfg["alpha_primary"])
    rows_all, y = ds.orig_row, ds.y[target]
    te = np.where(ds.sg == g)[0]
    tr = np.where(ds.sg != g)[0]
    lat_tr = ds.sg[tr]
    lats = np.unique(lat_tr)
    fold_l = np.searchsorted(lats, lat_tr).astype(np.int64)
    rows_tr, rows_te, y_tr, y_te = rows_all[tr], rows_all[te], y[tr], y[te]
    seed = _seed(cfg, target, g)

    # --- nested HPO: inner folds = leave one training lattice out ---
    folds = [(np.where(fold_l != k)[0], np.where(fold_l == k)[0]) for k in range(len(lats))]
    hp_json = os.path.join(paths.hpo, f"hp_{target}_loso{g}.json")
    if not os.path.exists(hp_json):
        print(f"  [HPO] {target} LOSO-{name}: {L['hpo']['n_trials']} trials x {len(lats)} lattice folds "
              f"(n={len(tr)})", flush=True)
    hp = tune(X_path, rows_tr, y_tr, folds, L["hpo"], cfg, seed=seed + 1,
              out_json=hp_json, log_json=os.path.join(paths.hpo, f"trials_{target}_loso{g}.json"),
              label=f"{target} LOSO-{name}")

    # --- fits ---
    workers, threads = plan_parallel(cfg)
    base = dict(X_path=X_path, max_bin=cfg["lgbm"]["max_bin"], n_rounds=hp["n_rounds"])
    fold_r = np.full(len(tr), -1, np.int64)
    sgkf = StratifiedGroupKFold(L["random_cv_k"], shuffle=True, random_state=seed + 2)
    for k, (_, v) in enumerate(sgkf.split(np.zeros(len(tr)), lat_tr, ds.group[tr])):
        fold_r[v] = k
    assert (fold_r >= 0).all()
    tasks, names = [], []
    for k in range(len(lats)):
        fit, oof = fold_l != k, fold_l == k
        tasks.append(dict(base, rows_train=rows_tr[fit], y_train=y_tr[fit],
                          params=lgb_params(hp, cfg, seed + 10 + k, threads),
                          predict={"oof": rows_tr[oof], "test": rows_te}))
        names.append(f"l{k}")
        if L["cqr"]:
            for nm, q in (("qlo", a_p / 2), ("qhi", 1 - a_p / 2)):
                tasks.append(dict(base, rows_train=rows_tr[fit], y_train=y_tr[fit],
                                  params=lgb_params(hp, cfg, seed + 10 + k, threads, quantile=q),
                                  predict={"oof": rows_tr[oof], "test": rows_te}))
                names.append(f"{nm}{k}")
    for k in range(L["random_cv_k"]):
        fit, oof = fold_r != k, fold_r == k
        tasks.append(dict(base, rows_train=rows_tr[fit], y_train=y_tr[fit],
                          params=lgb_params(hp, cfg, seed + 30 + k, threads),
                          predict={"oof": rows_tr[oof], "test": rows_te}))
        names.append(f"r{k}")
    res = dict(zip(names, run_tasks(tasks, workers)))

    K = len(lats)
    mu_te_l = np.stack([res[f"l{k}"]["test"] for k in range(K)])
    oof_l = np.empty(len(tr))
    for k in range(K):
        oof_l[fold_l == k] = res[f"l{k}"]["oof"]
    R_l = np.abs(y_tr - oof_l)
    Kr = L["random_cv_k"]
    mu_te_r = np.stack([res[f"r{k}"]["test"] for k in range(Kr)])
    oof_r = np.empty(len(tr))
    for k in range(Kr):
        oof_r[fold_r == k] = res[f"r{k}"]["oof"]
    R_r = np.abs(y_tr - oof_r)
    sig_oof, sig_te = _sigma_per_fold(X, rows_tr, fold_l, rows_te, cfg, seed + 40)
    if L["cqr"]:
        ql_te = np.stack([res[f"qlo{k}"]["test"] for k in range(K)])
        qh_te = np.stack([res[f"qhi{k}"]["test"] for k in range(K)])
        ql_oof, qh_oof = np.empty(len(tr)), np.empty(len(tr))
        for k in range(K):
            ql_oof[fold_l == k] = res[f"qlo{k}"]["oof"]
            qh_oof[fold_l == k] = res[f"qhi{k}"]["oof"]
        E = np.maximum(ql_oof - y_tr, y_tr - qh_oof)

    # --- intervals ---
    intervals = {}
    keep_sets = [("", np.ones(len(tr), bool))]
    if g != REF:
        keep_sets.append(("_noFd3m_cal", lat_tr != REF))
    for suffix, keep in keep_sets:
        f, Rk = fold_l[keep], R_l[keep]
        mu = mu_te_l[f]
        nsc = (R_l / sig_oof)[keep]
        for a in alphas:
            intervals[(f"cvplus_lattice{suffix}", a)] = cp.cvplus_interval(mu_te_l, f, Rk, a)
            intervals[(f"ncvplus_knn{suffix}", a)] = cp.cvplus_from_terms(
                mu - nsc[:, None] * sig_te[f], mu + nsc[:, None] * sig_te[f], a)
        if L["cqr"]:
            intervals[(f"cqr_cvplus{suffix}", a_p)] = cp.cvplus_from_terms(
                ql_te[f] - E[keep][:, None], qh_te[f] + E[keep][:, None], a_p)
    for a in alphas:
        intervals[("cvplus_random_in", a)] = cp.cvplus_interval(mu_te_r, fold_r, R_r, a)
    tau = float(L["group_conformal_tau"])
    fam = np.array([cp.higher_quantile(R_l[fold_l == k], tau) for k in range(K)])
    mubar_l, mubar_r = mu_te_l.mean(axis=0), mu_te_r.mean(axis=0)
    group_rows = []
    for a in alphas:
        T = cp.group_level_threshold(fam, a)
        intervals[("group_level", a)] = (mubar_l - T, mubar_l + T)
        cov = np.mean(np.abs(y_te - mubar_l) <= T)
        group_rows.append(dict(target=target, heldout=name, alpha=a, K=K, tau=tau,
                               k_needed=cp.k_upper(K, a), threshold=T, width=2 * T, coverage_of_heldout=cov,
                               heldout_tau_quantile=cp.higher_quantile(np.abs(y_te - mubar_l), tau),
                               family_event_covered=bool(cp.higher_quantile(np.abs(y_te - mubar_l), tau) <= T),
                               **{f"S_{LATTICE_NAME[l]}": s for l, s in zip(lats, fam)}))
    preds = {m: (mubar_r if m == "cvplus_random_in" else mubar_l) for m, _ in intervals}

    unit = dict(target=target, design="loso", heldout=name, sg=int(g), n_train=len(tr), n_test=len(te),
                mae_lattice_ensemble=float(np.mean(np.abs(y_te - mubar_l))),
                mae_random_ensemble=float(np.mean(np.abs(y_te - mubar_r))))
    df = pd.DataFrame(_metric_rows(unit, y_te, _masks(ds, te, g), intervals, preds))
    comp = pd.DataFrame(_calibration_composition(R_l, lat_tr, alphas, "lattice")
                        + _calibration_composition(R_r, lat_tr, alphas, "random"))
    comp.insert(0, "heldout", name)
    comp.insert(0, "target", target)

    arrays = dict(test_id=ds.ids[te], train_id=ds.ids[tr], y_test=y_te, y_train=y_tr, lat_train=lat_tr,
                  sub_test=ds.sub[te].astype("U4"), conflict_test=ds.conflict[te],
                  fold_l=fold_l, R_l=R_l, oof_l=oof_l, mu_te_l=mu_te_l, sig_oof=sig_oof, sig_te=sig_te,
                  fold_r=fold_r, R_r=R_r, mu_te_r=mu_te_r, family_scores=fam, lattices=lats)
    if L["cqr"]:
        arrays.update(ql_te=ql_te, qh_te=qh_te, E=E)
    for (m, a), (lo, hi) in intervals.items():
        arrays[f"{m}__a{a:.2f}__lo"], arrays[f"{m}__a{a:.2f}__hi"] = lo, hi
    atomic_save_npz(stem + ".npz", **arrays)
    atomic_write_csv(stem + "_calcomp.csv", comp)
    atomic_write_csv(stem + "_grouplevel.csv", pd.DataFrame(group_rows))
    atomic_write_json(stem + ".hp.json", hp)
    atomic_write_csv(stem + ".csv", df)
    c = df[(df.method == "cvplus_lattice") & (df.alpha.round(4) == round(a_p, 4)) & (df.group == "heldout")]
    return (f"done in {(time.time() - t0) / 60:.1f} min | MAE {unit['mae_lattice_ensemble']:.4f} | "
            f"CV+ coverage {float(c.coverage.iloc[0]):.3f}")


# ----------------------------------------------------------------------------- leave-two-out unit
def run_l2o_unit(ds, target, x, cfg, paths, X_path, X):
    name = LATTICE_NAME[x]
    stem = os.path.join(paths.loso, f"{target}_l2o{x}")
    if os.path.exists(stem + ".csv") and os.path.exists(stem + ".npz"):
        return "skipped"
    hp_json = os.path.join(paths.hpo, f"hp_{target}_loso{x}.json")
    if not os.path.exists(hp_json):
        raise RuntimeError(f"LOSO-{name} hyperparameters missing ({hp_json}); run LOSO-{name} first.")
    hp = json.load(open(hp_json))
    t0 = time.time()
    alphas = [float(a) for a in cfg["alpha"]]
    rows_all, y = ds.orig_row, ds.y[target]
    te = np.where(ds.sg == x)[0]
    tr = np.where((ds.sg != x) & (ds.sg != REF))[0]
    lat_tr = ds.sg[tr]
    lats = np.unique(lat_tr)
    fold = np.searchsorted(lats, lat_tr).astype(np.int64)
    rows_tr, rows_te, y_tr, y_te = rows_all[tr], rows_all[te], y[tr], y[te]
    seed = _seed(cfg, target, x, extra=60)
    workers, threads = plan_parallel(cfg)
    base = dict(X_path=X_path, max_bin=cfg["lgbm"]["max_bin"], n_rounds=hp["n_rounds"])
    tasks = [dict(base, rows_train=rows_tr[fold != k], y_train=y_tr[fold != k],
                  params=lgb_params(hp, cfg, seed + k, threads),
                  predict={"oof": rows_tr[fold == k], "test": rows_te}) for k in range(len(lats))]
    res = run_tasks(tasks, workers)
    mu_te = np.stack([r["test"] for r in res])
    oof = np.empty(len(tr))
    for k, r in enumerate(res):
        oof[fold == k] = r["oof"]
    R = np.abs(y_tr - oof)
    sig_oof, sig_te = _sigma_per_fold(X, rows_tr, fold, rows_te, cfg, seed + 40)
    mu = mu_te[fold]
    nsc = R / sig_oof
    intervals = {}
    for a in alphas:
        intervals[("cvplus_lattice", a)] = cp.cvplus_interval(mu_te, fold, R, a)
        intervals[("ncvplus_knn", a)] = cp.cvplus_from_terms(mu - nsc[:, None] * sig_te[fold],
                                                             mu + nsc[:, None] * sig_te[fold], a)
    mubar = mu_te.mean(axis=0)
    preds = {m: mubar for m, _ in intervals}
    unit = dict(target=target, design="leave_two_out", heldout=name, sg=int(x), n_train=len(tr), n_test=len(te),
                mae_lattice_ensemble=float(np.mean(np.abs(y_te - mubar))), mae_random_ensemble=np.nan)
    df = pd.DataFrame(_metric_rows(unit, y_te, _masks(ds, te, x), intervals, preds))
    arrays = dict(test_id=ds.ids[te], train_id=ds.ids[tr], y_test=y_te, lat_train=lat_tr, fold=fold, R=R,
                  mu_te=mu_te, sig_oof=sig_oof, sig_te=sig_te)
    for (m, a), (lo, hi) in intervals.items():
        arrays[f"{m}__a{a:.2f}__lo"], arrays[f"{m}__a{a:.2f}__hi"] = lo, hi
    atomic_save_npz(stem + ".npz", **arrays)
    atomic_write_csv(stem + ".csv", df)
    return f"done in {(time.time() - t0) / 60:.1f} min | MAE {unit['mae_lattice_ensemble']:.4f}"


# ----------------------------------------------------------------------------- driver
def run_phase2(cfg, targets=("bandgap", "formation"), allow_deviation=False):
    paths = get_paths(cfg)
    check_preregistration(cfg, paths, allow_deviation)
    if not os.path.exists(os.path.join(paths.overlap, "overlap_diagnostics.csv")):
        raise RuntimeError("Run the `overlap` step (Phase 3 diagnostics) before Phase 2, as pre-registered.")
    X_path = ensure_feature_cache(cfg, paths)
    X = load_features(cfg, paths)
    ds = build_dataset(cfg, paths)
    w, t = plan_parallel(cfg)
    order = [int(g) for g in cfg["loso"]["order"]]
    print(f"Phase 2 | n={ds.n} | targets={list(targets)} | holdouts={[LATTICE_NAME[g] for g in order]} | "
          f"{w} parallel fits x {t} threads", flush=True)
    t_start = time.time()
    units = [(tg, g) for tg in targets for g in order]
    for i, (target, g) in enumerate(units, 1):
        st = run_loso_unit(ds, target, g, cfg, paths, X_path, X)
        print(f"  [{i}/{len(units)}] {target:<9} LOSO-{LATTICE_NAME[g]:<8} {st}   "
              f"(elapsed {(time.time() - t_start) / 60:.0f} min)", flush=True)
        if cfg["loso"]["leave_two_out"] and g != REF:
            st = run_l2o_unit(ds, target, g, cfg, paths, X_path, X)
            print(f"        {target:<9} L2O-{LATTICE_NAME[g]:<8} {st}", flush=True)
    print("Phase 2 finished. Next: report2")
