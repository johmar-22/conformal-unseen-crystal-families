"""Phase 3: covariate-overlap diagnostics for weighted conformal prediction under LOSO.

Weighted conformal prediction (Tibshirani et al., NeurIPS 2019) and its jackknife+ form JAW
(Prinster et al., NeurIPS 2022) re-weight calibration points by the likelihood ratio
w(x) = dP_test(x) / dP_cal(x). With a probabilistic classifier separating the held-out lattice g
(label 1) from the five calibration lattices (label 0),
    w(x) = p(x) / (1 - p(x)) * n_rest / n_g.
Two things can make weighted CP unusable. (i) A few calibration points carry all the weight:
the Kish effective sample size ESS = (sum w)^2 / sum w^2 collapses. (ii) The test point's own
weight dominates: weighted CP places mass p_test = w(x) / (sum_i w_i + w(x)) at +inf, so the
interval is unbounded whenever p_test > alpha. With near-perfect separation (AUROC ~ 1) the
calibration weights can all be uniformly tiny, so (i) may look harmless while (ii) makes every
interval infinite; both are therefore reported and both enter the gate.

Out-of-fold probabilities come from composition-grouped folds (no composition group in both the
fitting and the scoring part). Two classifiers are used, logistic regression (linear) and kNN
(local), on PCA components of the standardized SOAP vector; the gate uses the larger AUROC.
The same machinery checks whether the two Fd-3m cell settings (227a, 227b) are separable.
Everything here is written once, before any LOSO result exists (see the pre-registration).
"""
from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .config import atomic_save_npz, atomic_write_csv, atomic_write_json, get_paths
from .conformal import kish_ess
from .data import LATTICE_NAME, build_dataset, load_features
from .setup2 import check_preregistration


def pca_scores(X, n_components: int, seed: int) -> tuple[np.ndarray, float]:
    Z = StandardScaler().fit_transform(np.asarray(X, np.float64))
    p = PCA(n_components=n_components, svd_solver="randomized", random_state=seed).fit(Z)
    return p.transform(Z), float(p.explained_variance_ratio_.sum())


def oof_proba(Z, y, groups, model, folds: int, seed: int) -> np.ndarray:
    p = np.full(len(y), np.nan)
    for tr, te in StratifiedGroupKFold(folds, shuffle=True, random_state=seed).split(Z, y, groups):
        p[te] = model().fit(Z[tr], y[tr]).predict_proba(Z[te])[:, 1]
    assert np.isfinite(p).all()
    return p


def run_overlap(cfg, allow_deviation=False):
    paths = get_paths(cfg)
    check_preregistration(cfg, paths, allow_deviation)
    f_csv = os.path.join(paths.overlap, "overlap_diagnostics.csv")
    if os.path.exists(f_csv):
        print(f"  overlap diagnostics already exist (kept unchanged): {f_csv}")
        print(pd.read_csv(f_csv).to_string(index=False))
        return
    o = cfg["overlap"]
    seed = int(cfg["seed_base"]) + 7000
    ds = build_dataset(cfg, paths)
    X = load_features(cfg, paths)
    Z, expl = pca_scores(X[ds.orig_row], o["pca_components"], seed)
    print(f"  PCA({o['pca_components']}) on standardized SOAP: explained variance {expl:.3f}", flush=True)
    models = {
        "logistic": lambda: make_pipeline(StandardScaler(), LogisticRegression(C=o["logistic_C"], max_iter=5000)),
        "knn": lambda: make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=o["knn_k"])),
    }
    clip = float(o["prob_clip"])
    a_p = float(cfg["alpha_primary"])
    warnings.filterwarnings("ignore", message="Unknown solver options")
    rows, arrays = [], {"id": ds.ids}
    for g in cfg["loso"]["order"]:
        y = (ds.sg == g).astype(int)
        res = dict(heldout=LATTICE_NAME[g], sg=int(g), n_g=int(y.sum()), n_rest=int((1 - y).sum()))
        probs = {}
        for name, mk in models.items():
            p = oof_proba(Z, y, ds.group, mk, o["cv_folds"], seed + g)
            probs[name] = p
            res[f"auroc_{name}"] = float(roc_auc_score(y, p))
        res["auroc_max"] = max(res["auroc_logistic"], res["auroc_knn"])
        p = np.clip(probs["logistic"], clip, 1 - clip)   # smooth density-ratio estimate for the weights
        w = p / (1 - p) * res["n_rest"] / res["n_g"]
        w_rest = w[y == 0]
        res["ess_rest"] = kish_ess(w_rest)
        res["ess_frac"] = res["ess_rest"] / res["n_rest"]
        res["max_weight_share"] = float(w_rest.max() / w_rest.sum())
        share = w[y == 1] / (w_rest.sum() + w[y == 1])
        res["median_test_weight_share"] = float(np.median(share))
        res["frac_test_unbounded"] = float(np.mean(share > a_p))   # weighted interval infinite at alpha_primary
        res["overlap_ok"] = bool(res["auroc_max"] <= o["auroc_max"] and res["ess_rest"] >= o["ess_min"]
                                 and res["median_test_weight_share"] < a_p)
        rows.append(res)
        arrays[f"p_logistic_{g}"] = probs["logistic"]
        arrays[f"p_knn_{g}"] = probs["knn"]
        arrays[f"w_{g}"] = w
        print(f"  holdout {LATTICE_NAME[g]:<8} AUROC logistic {res['auroc_logistic']:.3f}  kNN {res['auroc_knn']:.3f}  "
              f"ESS {res['ess_rest']:.1f} of {res['n_rest']}  test-weight share {res['median_test_weight_share']:.3f}  "
              f"overlap_ok={res['overlap_ok']}", flush=True)
    # Fd-3m cell settings: are 227a and 227b distinguishable in SOAP space?
    m = ds.sg == 227
    y = (ds.sub[m] == "227b").astype(int)
    sub = dict(heldout="227a vs 227b (within Fd-3m)", sg=227, n_g=int(y.sum()), n_rest=int((1 - y).sum()))
    for name, mk in models.items():
        sub[f"auroc_{name}"] = float(roc_auc_score(y, oof_proba(Z[m], y, ds.group[m], mk, o["cv_folds"], seed)))
    sub["auroc_max"] = max(sub["auroc_logistic"], sub["auroc_knn"])
    rows.append(sub)
    print(f"  Fd-3m 227a vs 227b: AUROC logistic {sub['auroc_logistic']:.3f}  kNN {sub['auroc_knn']:.3f}")
    df = pd.DataFrame(rows)
    atomic_save_npz(os.path.join(paths.overlap, "overlap_weights.npz"), **arrays)
    atomic_write_json(os.path.join(paths.overlap, "overlap_decision.json"),
                      dict(pca_explained=expl, gate=dict(auroc_max=o["auroc_max"], ess_min=o["ess_min"],
                                                         median_test_weight_share_below=a_p),
                           weighted_cvplus_for=[r["heldout"] for r in rows if r.get("overlap_ok")]))
    atomic_write_csv(f_csv, df)
    print(df.round(3).to_string(index=False))
