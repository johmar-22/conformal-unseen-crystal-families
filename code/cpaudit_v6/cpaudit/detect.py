"""Failure-detection scores and their evaluation (Phase 4).

Scores (higher = more doubt about the prediction):
  width       CV+ interval width of the reference method (what the submitted manuscript reports)
  knn         mean Euclidean distance to the k nearest training structures, standardized SOAP,
              PCA-50 (the same distance as the Phase 1-2 normalized score; cf. Janet et al.,
              Chem. Sci. 10, 7913 (2019) for distance-to-data as an uncertainty metric)
  kde         kernel-density dissimilarity d = 1 - KDE(x) / max_train KDE, Epanechnikov kernel,
              standardized features, bandwidth from the nearest-neighbour rule of scikit-learn's
              estimate_bandwidth (Schultz, Wang, Jacobs & Morgan, npj Comput. Mater. 2025)
  committee   standard deviation of the lattice-fold models' predictions (leave-one-family-out
              committee; disagreement between models that each miss one training lattice)
  subsample   standard deviation of the five random-fold CV+ models' predictions (an ensemble of
              80 % grouped subsamples; used instead of a bootstrap ensemble so that no model is
              refitted; ensemble spread as in Lakshminarayanan et al., NeurIPS 2017)
  classifier  out-of-fold probability that a structure belongs to the screened (held-out) set,
              from a domain classifier trained on training-vs-screened structures (transductive:
              it uses the unlabelled screening set, as a practitioner could)

Evaluation: AUROC (Mann-Whitney, ties = 1/2), average precision, cluster bootstrap over
composition groups (the exchangeability unit established in Phase 1), selective-prediction
(abstention) curves (Geifman & El-Yaniv, NeurIPS 2017), and conformal p-values for outlier
flagging (Bates, Candes, Lei, Romano & Sesia, Ann. Stat. 51, 149, 2023).
"""
from __future__ import annotations

import math
import warnings

import numpy as np
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .difficulty import KNNDifficulty

SCORES = ["width", "knn", "kde", "committee", "subsample", "classifier"]
DISTANCE_SCORES = ["knn", "kde"]


# ----------------------------------------------------------------------------- metrics
def auroc(y, s) -> float:
    """P(score of a random positive > score of a random negative), ties counted 1/2."""
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def avg_precision(y, s) -> float:
    y = np.asarray(y, bool)
    return float(average_precision_score(y, s)) if 0 < y.sum() < y.size else np.nan


class ClusterBootstrap:
    """Resample composition groups with replacement (Field & Welsh 2007; cluster bootstrap)."""

    def __init__(self, groups, seed: int):
        g = np.asarray(groups)
        self.idx = [np.flatnonzero(g == u) for u in np.unique(g)]
        self.rng = np.random.default_rng(seed)

    def draw(self) -> np.ndarray:
        pick = self.rng.integers(0, len(self.idx), len(self.idx))
        return np.concatenate([self.idx[i] for i in pick])


def percentile_ci(x, level=0.95):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (np.nan, np.nan)
    a = (1 - level) / 2
    return float(np.quantile(x, a)), float(np.quantile(x, 1 - a))


# ----------------------------------------------------------------------------- scores
def _sqdist(A: np.ndarray, B: np.ndarray, chunk: int = 512) -> np.ndarray:
    """Squared Euclidean distances, float64, computed in row chunks of A."""
    bb = np.einsum("ij,ij->i", B, B)
    out = np.empty((A.shape[0], B.shape[0]))
    for i in range(0, A.shape[0], chunk):
        a = A[i:i + chunk]
        out[i:i + chunk] = np.einsum("ij,ij->i", a, a)[:, None] + bb[None, :] - 2.0 * a @ B.T
    np.maximum(out, 0.0, out=out)
    return out


def estimate_bandwidth_nn(D2_train: np.ndarray, quantile: float) -> float:
    """scikit-learn estimate_bandwidth: mean over training points of the distance to their
    n_neighbors-th nearest neighbour (self included as the first), n_neighbors = int(n*quantile)."""
    n = D2_train.shape[0]
    k = max(1, int(n * quantile))
    kth = np.partition(D2_train, k - 1, axis=1)[:, k - 1]
    return float(np.sqrt(kth).mean())


def kde_dissimilarity(X_train, X_test, quantile: float = 0.3):
    """Schultz et al. (2025): d = 1 - KDE(x)/max_i KDE(x_i), Epanechnikov kernel on standardized
    features. Returns (d_test, d_train, bandwidth)."""
    sc = StandardScaler().fit(X_train)
    A = sc.transform(X_train).astype(np.float64)
    B = sc.transform(X_test).astype(np.float64)
    D2tr = _sqdist(A, A)
    h = estimate_bandwidth_nn(D2tr, quantile)
    dens_tr = np.clip(1.0 - D2tr / h ** 2, 0, None).sum(axis=1)
    del D2tr
    dens_te = np.clip(1.0 - _sqdist(B, A) / h ** 2, 0, None).sum(axis=1)
    m = dens_tr.max()
    return 1.0 - dens_te / m, 1.0 - dens_tr / m, h


def knn_distance(X_train, X_test, k: int, n_components: int, seed: int):
    knn = KNNDifficulty(k=k, n_components=n_components, beta=0.0, seed=seed).fit(X_train)
    return knn.sigma(X_test), knn


def classifier_oof(X_train, X_test, groups_train, groups_test, n_components, C, folds, seed):
    """Out-of-fold P(screened | x) for the screened structures, grouped folds."""
    warnings.filterwarnings("ignore", message="Unknown solver options")
    X = np.vstack([X_train, X_test]).astype(np.float64)
    y = np.r_[np.zeros(len(X_train), int), np.ones(len(X_test), int)]
    g = np.r_[np.asarray(groups_train), np.asarray(groups_test)]
    Z = PCA(n_components=n_components, svd_solver="randomized", random_state=seed).fit_transform(
        StandardScaler().fit_transform(X))
    p = np.full(len(y), np.nan)
    for tr, te in StratifiedGroupKFold(folds, shuffle=True, random_state=seed).split(Z, y, g):
        m = make_pipeline(StandardScaler(), LogisticRegression(C=C, max_iter=5000)).fit(Z[tr], y[tr])
        p[te] = m.predict_proba(Z[te])[:, 1]
    return p[len(X_train):], p[:len(X_train)]


def conformal_pvalues(cal_scores, test_scores) -> np.ndarray:
    """Marginal conformal p-values (Bates et al. 2023, eq. 1): p = (1 + #{cal >= s}) / (n + 1)."""
    c = np.sort(np.asarray(cal_scores, float))
    s = np.asarray(test_scores, float)
    n_ge = c.size - np.searchsorted(c, s, side="left")
    return (1.0 + n_ge) / (c.size + 1.0)


# ----------------------------------------------------------------------------- abstention
def abstention_curve(score, covered, abs_err, grid, seed: int):
    """Abstain on the highest-score fraction f; report interval coverage and MAE of the rest.
    Ties are broken at random (fixed seed) so a saturated score is not favoured by order."""
    rng = np.random.default_rng(seed)
    order = np.lexsort((rng.random(len(score)), np.asarray(score, float)))   # ascending score
    cov, err = np.asarray(covered, float)[order], np.asarray(abs_err, float)[order]
    n = len(order)
    out = []
    for f in grid:
        keep = max(1, int(math.ceil((1 - f) * n)))
        out.append((float(f), float(cov[:keep].mean()), float(err[:keep].mean())))
    return out
