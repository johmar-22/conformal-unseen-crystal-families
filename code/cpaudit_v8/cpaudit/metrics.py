"""Interval metrics: coverage with exact CIs, width, interval score, width-error association,
and the sampling distribution of coverage expected from a valid split-conformal procedure.

Interval score (Gneiting & Raftery, JASA 102, 359, 2007), lower is better:
    IS = (u - l) + (2/alpha)(l - y) 1{y < l} + (2/alpha)(y - u) 1{y > u}
Training-conditional coverage of split conformal (Vovk 2013; Angelopoulos & Bates 2023):
    C | calibration ~ Beta(k, n + 1 - k),  k = ceil((1 - alpha)(n + 1)),
so test coverage on m points is Beta-binomial(m, k, n + 1 - k) / m.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import beta as beta_dist
from scipy.stats import betabinom, spearmanr

from .conformal import k_upper


def clopper_pearson(k: int, n: int, level: float = 0.95):
    if n == 0:
        return np.nan, np.nan
    a = 1 - level
    lo = 0.0 if k == 0 else float(beta_dist.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta_dist.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def interval_score(y, lo, hi, alpha):
    y, lo, hi = map(lambda a: np.asarray(a, float), (y, lo, hi))
    return (hi - lo) + (2 / alpha) * np.clip(lo - y, 0, None) + (2 / alpha) * np.clip(y - hi, 0, None)


def summarize(y, pred, lo, hi, alpha) -> dict:
    y, pred, lo, hi = map(lambda a: np.asarray(a, float), (y, pred, lo, hi))
    n = y.size
    cov = (y >= lo) & (y <= hi)
    k = int(cov.sum())
    w = hi - lo
    finite = np.isfinite(w)
    ci = clopper_pearson(k, n)
    err = np.abs(y - pred)
    if finite.sum() >= 10 and np.nanstd(w[finite]) > 1e-12:
        rho = float(spearmanr(w[finite], err[finite]).statistic)
    else:
        rho = np.nan  # constant width (e.g. global split conformal): association undefined
    return dict(
        n=n, covered=k, coverage=k / n if n else np.nan, cov_ci_lo=ci[0], cov_ci_hi=ci[1],
        width_mean=float(np.mean(w)) if finite.all() else np.inf,
        width_median=float(np.median(w)),
        frac_infinite=float(1 - finite.mean()) if n else np.nan,
        interval_score=float(np.mean(interval_score(y, lo, hi, alpha))) if finite.all() else np.inf,
        mae=float(np.mean(err)), spearman_width_abserr=rho,
    )


def betabinom_prob_below(n_cal: int, n_test: int, alpha: float, threshold: float | None = None) -> float:
    """P(test coverage < threshold) for a valid split-conformal interval (threshold = 1 - alpha)."""
    if n_cal <= 0 or n_test <= 0:
        return np.nan
    k = k_upper(n_cal, alpha)
    if k > n_cal:
        return 0.0  # infinite interval always covers
    thr = 1 - alpha if threshold is None else threshold
    m_crit = int(np.ceil(thr * n_test - 1e-9))   # covering fewer than m_crit points -> below
    return float(betabinom.cdf(m_crit - 1, n_test, k, n_cal + 1 - k))
