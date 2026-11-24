"""Phase 4: failure-detection scores for every held-out structure (no model is refitted).

Configurations: the 6 nested-LOSO holdouts and the 5 leave-two-out holdouts of Phase 2. For each
configuration the training set is the set of structures the Phase 2 models were trained on, and
the screened set is the held-out lattice. Feature-based scores (knn, kde, classifier) depend only
on these two sets and are shared by both targets; width, committee and subsample come from the
saved Phase 2 predictions of each target.

Outcomes per structure:
  miscovered    y outside the reference interval (cvplus_lattice at alpha_primary, the design of the
                submitted manuscript); split into miss_below / miss_above
  large_error   |error| above the in-distribution 90th percentile, i.e. the 0.9-quantile of the
                random-fold out-of-fold |residuals| on the same training lattices (LOSO only)

Flag (pre-registered, LOSO and leave-two-out): the training structures are split by composition
group into fit / calibration / in-distribution test (3 : 1 : 1 folds); the kNN distance is fitted
on the fit part; conformal p-values (Bates et al. 2023) are computed from the calibration part;
structures with p <= alpha_flag are flagged. The in-distribution test part measures the
false-flag rate, which is <= alpha_flag in expectation when calibration and test structures are
exchangeable.
"""
from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

from .config import atomic_save_npz, atomic_write_csv, get_paths
from .data import LATTICE_NAME, id_table, load_features
from .detect import classifier_oof, conformal_pvalues, kde_dissimilarity, knn_distance
from .setup4 import check_preregistration

REF = 227
TARGETS = ("bandgap", "formation")


def _configs(paths):
    out = []
    for g in (227, 12, 33, 167, 194, 206):
        out.append(("loso", g))
    for g in (12, 33, 167, 194, 206):
        out.append(("l2o", g))
    for design, g in out:
        for t in TARGETS:
            f = os.path.join(paths.phase2_loso, f"{t}_{design}{g}.npz")
            if not os.path.exists(f):
                raise FileNotFoundError(f"Phase 2 unit missing: {f}")
    return out


def _feature_scores(design, g, paths, cfg, X, tab):
    """kNN, KDE and classifier scores for one configuration (cached)."""
    f = os.path.join(paths.scores, f"features_{design}{g}.npz")
    if os.path.exists(f):
        return dict(np.load(f))
    d = cfg["detect"]
    z = np.load(os.path.join(paths.phase2_loso, f"bandgap_{design}{g}.npz"))
    tr_id, te_id = z["train_id"], z["test_id"]
    zf = np.load(os.path.join(paths.phase2_loso, f"formation_{design}{g}.npz"))
    if not (np.array_equal(tr_id, zf["train_id"]) and np.array_equal(te_id, zf["test_id"])):
        raise AssertionError(f"band gap and formation units of {design}{g} use different structures")
    seed = int(cfg["seed_base"]) + 100 * g + (0 if design == "loso" else 7)
    Xtr = np.asarray(X[tr_id - 1], np.float64)
    Xte = np.asarray(X[te_id - 1], np.float64)
    t0 = time.time()
    knn_te, knn = knn_distance(Xtr, Xte, d["knn"]["k"], d["knn"]["pca_components"], seed)
    kde_te, kde_tr, h = kde_dissimilarity(Xtr, Xte, d["kde"]["bandwidth_quantile"])
    c = d["classifier"]
    clf_te, clf_tr = classifier_oof(Xtr, Xte, tab.loc[tr_id, "group"].to_numpy(),
                                    tab.loc[te_id, "group"].to_numpy(), c["pca_components"],
                                    c["logistic_C"], c["cv_folds"], seed)
    # conformal flag: fit / calibration / in-distribution test split of the training structures
    fl = d["flag"]
    lat_tr = tab.loc[tr_id, "sg"].to_numpy()
    fold = np.full(len(tr_id), -1)
    for k, (_, v) in enumerate(StratifiedGroupKFold(fl["split_folds"], shuffle=True, random_state=seed + 1)
                               .split(np.zeros(len(tr_id)), lat_tr, tab.loc[tr_id, "group"].to_numpy())):
        fold[v] = k
    fit, cal, idt = fold >= 2, fold == 0, fold == 1
    s_te, det = knn_distance(Xtr[fit], Xte, d["knn"]["k"], d["knn"]["pca_components"], seed + 2)
    s_cal = det.sigma(Xtr[cal])
    s_idt = det.sigma(Xtr[idt])
    p_te = conformal_pvalues(s_cal, s_te)
    p_idt = conformal_pvalues(s_cal, s_idt)
    out = dict(train_id=tr_id, test_id=te_id, knn=knn_te, kde=kde_te, classifier=clf_te,
               kde_train=kde_tr, classifier_train=clf_tr, kde_bandwidth=np.array(h),
               knn_explained=np.array(knn.explained_), flag_p_test=p_te, flag_p_idtest=p_idt,
               flag_idtest_id=tr_id[idt], flag_n_cal=np.array(int(cal.sum())))
    atomic_save_npz(f, **out)
    print(f"    features {design}-{LATTICE_NAME[g]}: n_train={len(tr_id)} n_test={len(te_id)} "
          f"kde bandwidth={h:.1f} PCA explained={knn.explained_:.3f} ({time.time() - t0:.0f} s)", flush=True)
    return out


def _structure_table(design, g, target, feat, paths, cfg, tab):
    z = np.load(os.path.join(paths.phase2_loso, f"{target}_{design}{g}.npz"))
    if not np.array_equal(z["test_id"], feat["test_id"]):
        raise AssertionError(f"structure order differs between Phase 2 unit and features ({design}{g})")
    a = float(cfg["alpha_primary"])
    ref = cfg["reference_interval"]
    lo, hi = z[f"{ref}__a{a:.2f}__lo"], z[f"{ref}__a{a:.2f}__hi"]
    y = z["y_test"]
    mu = z["mu_te_l"] if design == "loso" else z["mu_te"]
    pred = mu.mean(axis=0)
    err = pred - y
    te = z["test_id"]
    df = pd.DataFrame(dict(design=design, target=target, heldout=LATTICE_NAME[g], sg=g, id=te,
                           group=tab.loc[te, "group"].to_numpy(), sub=tab.loc[te, "sub"].to_numpy(),
                           y=y, pred=pred, err=err, abs_err=np.abs(err), lo=lo, hi=hi,
                           miss_below=y < lo, miss_above=y > hi))
    df["miscovered"] = df.miss_below | df.miss_above
    df["width"] = hi - lo
    df["committee"] = mu.std(axis=0)
    if design == "loso":
        df["subsample"] = z["mu_te_r"].std(axis=0)
        tau = float(np.quantile(z["R_r"], cfg["detect"]["large_error_quantile"]))
        df["large_error"] = df.abs_err > tau
        df["large_error_threshold"] = tau
    else:
        df["subsample"] = np.nan
        df["large_error"] = np.nan
        df["large_error_threshold"] = np.nan
    for s in ("knn", "kde", "classifier"):
        df[s] = feat[s]
    df["flag_p"] = feat["flag_p_test"]
    df["flagged"] = df.flag_p <= cfg["detect"]["flag"]["alpha"]
    return df


def run_phase4(cfg, allow_deviation=False):
    paths = get_paths(cfg)
    check_preregistration(cfg, paths, allow_deviation)
    tab = id_table(cfg, paths)
    X = load_features(cfg, paths)
    configs = _configs(paths)
    print(f"Phase 4 | {len(configs)} configurations x {len(TARGETS)} targets | features from "
          f"{cfg['paths']['soap_npz']}", flush=True)
    frames, flags = [], []
    t0 = time.time()
    for i, (design, g) in enumerate(configs, 1):
        feat = _feature_scores(design, g, paths, cfg, X, tab)
        for t in TARGETS:
            frames.append(_structure_table(design, g, t, feat, paths, cfg, tab))
        a = cfg["detect"]["flag"]["alpha"]
        flags.append(dict(design=design, heldout=LATTICE_NAME[g], n_cal=int(feat["flag_n_cal"]),
                          n_idtest=len(feat["flag_p_idtest"]),
                          false_flag_rate=float(np.mean(feat["flag_p_idtest"] <= a)),
                          n_test=len(feat["flag_p_test"]),
                          flag_rate_heldout=float(np.mean(feat["flag_p_test"] <= a))))
        print(f"  [{i}/{len(configs)}] {design}-{LATTICE_NAME[g]} done (elapsed {(time.time() - t0) / 60:.1f} min)",
              flush=True)
    table = pd.concat(frames, ignore_index=True)
    atomic_write_csv(os.path.join(paths.scores, "phase4_structures.csv"), table)
    atomic_write_csv(os.path.join(paths.scores, "phase4_flags_by_config.csv"), pd.DataFrame(flags))
    print(pd.DataFrame(flags).round(3).to_string(index=False))
    print("Phase 4 scores written. Next: report4")
