"""Conformal prediction primitives, written out explicitly so every formula can be checked.

References
----------
Split conformal ........ Lei, G'Sell, Rinaldo, Tibshirani & Wasserman, JASA 113, 1094 (2018).
Mondrian (per group) ... Vovk, Mach. Learn. 92, 349 (2013); Bostrom et al., PMLR 179, 24 (2022).
CV+ .................... Barber, Candes, Ramdas & Tibshirani, Ann. Stat. 49, 486 (2021), Sec. 3-4.
Normalized scores ...... Papadopoulos, Vovk & Gammerman, JAIR 40, 815 (2011); Lei et al. (2018).
CQR .................... Romano, Patterson & Candes, NeurIPS 32 (2019).
CV+ for general scores . Gupta, Kuchibhotla & Ramdas, Pattern Recognit. 127, 108496 (2022)
                         (nested conformal: CV+/jackknife+ for normalized and CQR scores).
Weighted (covariate .... Tibshirani, Barber, Candes & Ramdas, NeurIPS 32 (2019); jackknife+ form
 shift)                  (JAW): Prinster, Liu & Saria, NeurIPS 35 (2022).
Group-level ............ Dunn, Wasserman & Ramdas, JASA 118, 2491 (2023).

Conventions (alpha = miscoverage, n = number of calibration scores):
  k_hi = ceil((1 - alpha)(n + 1));  if k_hi > n the upper bound is +inf.
  k_lo = floor(alpha (n + 1));      if k_lo < 1 the lower bound is -inf.
"""
from __future__ import annotations

import math

import numpy as np


def k_upper(n: int, alpha: float) -> int:
    # small tolerance so (n + 1)(1 - alpha) that is an integer in exact arithmetic is not
    # pushed up by binary rounding (e.g. 10 * 0.9 = 9.000000000000002)
    return int(math.ceil((n + 1) * (1.0 - alpha) - 1e-9))


def k_lower(n: int, alpha: float) -> int:
    return int(math.floor(alpha * (n + 1) + 1e-9))


def conformal_quantile(scores, alpha: float) -> float:
    """k_hi-th smallest score (1-indexed); +inf if there are too few scores for this alpha."""
    s = np.sort(np.asarray(scores, dtype=float))
    n = s.size
    if n == 0:
        return np.inf
    k = k_upper(n, alpha)
    return np.inf if k > n else float(s[k - 1])


# ----------------------------------------------------------------------------- split / Mondrian
def split_interval(pred_test, cal_scores, alpha):
    q = conformal_quantile(cal_scores, alpha)
    return pred_test - q, pred_test + q, q


def mondrian_interval(pred_test, cat_test, cal_scores, cat_cal, alpha):
    """Separate quantile per category; a category absent from calibration gets +inf."""
    pred_test = np.asarray(pred_test, float)
    q = np.empty_like(pred_test)
    qmap = {}
    for c in np.unique(cat_test):
        qmap[c] = conformal_quantile(np.asarray(cal_scores)[np.asarray(cat_cal) == c], alpha)
        q[np.asarray(cat_test) == c] = qmap[c]
    return pred_test - q, pred_test + q, qmap


def normalized_interval(pred_test, sigma_test, y_cal, pred_cal, sigma_cal, alpha):
    scores = np.abs(np.asarray(y_cal) - np.asarray(pred_cal)) / np.asarray(sigma_cal)
    q = conformal_quantile(scores, alpha)
    return pred_test - q * sigma_test, pred_test + q * sigma_test, q


def cqr_interval(qlo_test, qhi_test, y_cal, qlo_cal, qhi_cal, alpha):
    scores = np.maximum(np.asarray(qlo_cal) - y_cal, y_cal - np.asarray(qhi_cal))
    q = conformal_quantile(scores, alpha)
    return qlo_test - q, qhi_test + q, q


# ----------------------------------------------------------------------------- CV+
def cvplus_interval(fold_pred_test: np.ndarray, fold_of_train: np.ndarray,
                    resid_train: np.ndarray, alpha: float):
    """CV+ interval as defined by Barber et al. (2021), Section 4 (K-fold version of jackknife+).

    fold_pred_test : (K, n_test) prediction of the model trained without fold k
    fold_of_train  : (n,) fold label of each training point
    resid_train    : (n,) out-of-fold absolute residual R_i = |y_i - mu_{-k(i)}(x_i)|

    lower = floor(alpha(n+1))-th smallest of  mu_{-k(i)}(x) - R_i
    upper = ceil((1-alpha)(n+1))-th smallest of mu_{-k(i)}(x) + R_i
    """
    fold_of_train = np.asarray(fold_of_train)
    R = np.asarray(resid_train, float)
    n = R.size
    mu = np.asarray(fold_pred_test, float)[fold_of_train, :]          # (n, n_test)
    lo_vals = mu - R[:, None]
    hi_vals = mu + R[:, None]
    kl, kh = k_lower(n, alpha), k_upper(n, alpha)
    lo = np.full(mu.shape[1], -np.inf) if kl < 1 else np.partition(lo_vals, kl - 1, axis=0)[kl - 1]
    hi = np.full(mu.shape[1], np.inf) if kh > n else np.partition(hi_vals, kh - 1, axis=0)[kh - 1]
    return lo, hi


def cvplus_from_terms(lo_terms: np.ndarray, hi_terms: np.ndarray, alpha: float):
    """CV+ for a general nested score (Gupta et al. 2022, Sec. 3).

    lo_terms[i, j] = lower end of the set that just contains calibration point i, evaluated at
    test point j with the model of fold k(i); hi_terms likewise. For absolute residuals these are
    mu_{-k(i)}(x_j) -+ R_i; for normalized scores mu -+ R_i sigma_{-k(i)}(x_j); for CQR
    q_lo,-k(i)(x_j) - E_i and q_hi,-k(i)(x_j) + E_i.
    """
    n = lo_terms.shape[0]
    kl, kh = k_lower(n, alpha), k_upper(n, alpha)
    lo = np.full(lo_terms.shape[1], -np.inf) if kl < 1 else np.partition(lo_terms, kl - 1, axis=0)[kl - 1]
    hi = np.full(hi_terms.shape[1], np.inf) if kh > n else np.partition(hi_terms, kh - 1, axis=0)[kh - 1]
    return lo, hi


def _weighted_upper(values: np.ndarray, p: np.ndarray, p_inf: float, level: float) -> float:
    """Smallest v with sum_{values <= v} p >= level, where an extra atom of mass p_inf sits at +inf."""
    o = np.argsort(values, kind="mergesort")
    c = np.cumsum(p[o])
    idx = np.searchsorted(c, level - 1e-12)
    return float(values[o][idx]) if idx < values.size else np.inf


def weighted_cvplus_interval(fold_pred_test, fold_of_train, resid_train, w_train, w_test, alpha):
    """CV+ with likelihood-ratio weights (JAW form, Prinster et al. 2022, applied to K folds).

    For test point j the calibration terms get mass w_i / (sum_i w_i + w_j) and a point mass
    w_j / (...) is placed at +inf (upper end) and -inf (lower end), as in weighted conformal
    prediction (Tibshirani et al. 2019). With all weights equal this is exactly cvplus_interval.
    The finite-sample guarantee of JAW is proven for the leave-one-out version; for K folds the
    weighting is used here as a diagnostic, not as a guarantee.
    """
    R = np.asarray(resid_train, float)
    mu = np.asarray(fold_pred_test, float)[np.asarray(fold_of_train), :]
    w = np.asarray(w_train, float)
    m = mu.shape[1]
    lo, hi = np.empty(m), np.empty(m)
    for j in range(m):
        tot = w.sum() + float(w_test[j])
        p, p_inf = w / tot, float(w_test[j]) / tot
        hi[j] = _weighted_upper(mu[:, j] + R, p, p_inf, 1 - alpha)
        lo[j] = -_weighted_upper(-(mu[:, j] - R), p, p_inf, 1 - alpha)
    return lo, hi


def group_level_threshold(family_scores, alpha: float) -> float:
    """k-th smallest of K family-level scores, k = ceil((K+1)(1-alpha)); +inf if k > K.
    Valid for a new family when the K + 1 families are exchangeable (Dunn et al. 2023)."""
    s = np.sort(np.asarray(family_scores, float))
    k = k_upper(s.size, alpha)
    return np.inf if k > s.size else float(s[k - 1])


def kish_ess(w) -> float:
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / np.sum(w ** 2)) if w.size and np.sum(w ** 2) > 0 else 0.0


def higher_quantile(x, tau: float) -> float:
    """ceil(tau * n)-th smallest value (a conservative empirical tau-quantile)."""
    s = np.sort(np.asarray(x, float))
    k = max(1, int(math.ceil(tau * s.size - 1e-9)))
    return float(s[k - 1])
