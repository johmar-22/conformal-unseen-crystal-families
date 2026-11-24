"""Conformal prediction primitives, written out explicitly so every formula can be checked.

References
----------
Split conformal ........ Lei, G'Sell, Rinaldo, Tibshirani & Wasserman, JASA 113, 1094 (2018).
Mondrian (per group) ... Vovk, Mach. Learn. 92, 349 (2013); Bostrom et al., PMLR 179, 24 (2022).
CV+ .................... Barber, Candes, Ramdas & Tibshirani, Ann. Stat. 49, 486 (2021), Sec. 3-4.
Normalized scores ...... Papadopoulos, Vovk & Gammerman, JAIR 40, 815 (2011); Lei et al. (2018).
CQR .................... Romano, Patterson & Candes, NeurIPS 32 (2019).

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
