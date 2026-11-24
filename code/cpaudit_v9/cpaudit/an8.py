"""Phase 6 Task 6.5 (v8): analysis tools that do not fit models. Each has a direct twin or a
known-answer test in selftest.py (t_phase8).

Contents
--------
kde_dissimilarity_sub   Schultz et al. (npj Comput. Mater. 2025) KDE dissimilarity, as in NOMAD Phase 4
                        (detect.kde_dissimilarity), with the reference set and the bandwidth sample
                        drawn at random when the training set is larger than n_ref / n_bw. With
                        n_train <= both limits it is identical to the Phase 4 function.
hjplus_definition       Hierarchical jackknife+ written directly from Lee, Barber & Willett
                        (arXiv:2306.06342 v4, Appendix C.2.2, eq. 20); the production code uses
                        cvfast.cvplus with one fold per group, which this twin checks.
spearman / boot_spearman / perm_spearman_mc
                        group-level rank correlation, bootstrap over groups (percentile CI) and a
                        Monte Carlo permutation p-value (Phipson & Smyth 2010 correction).
tail_enrichment         share of the calibration tail contributed by each group (H2, descriptive).
sg_tex                  space-group symbol in mathtext (F-43m -> $F\\bar{4}3m$).
"""
from __future__ import annotations

import re

import numpy as np
from scipy.stats import rankdata
from sklearn.preprocessing import StandardScaler

from .detect import _sqdist, estimate_bandwidth_nn

TOL = 1e-12


# ----------------------------------------------------------------------------- KDE
def kde_dissimilarity_sub(X_train, X_test, quantile=0.3, n_ref=20000, n_bw=5000, seed=0, chunk=512):
    """d = 1 - KDE(x) / max_r KDE(x_r), Epanechnikov kernel on standardized features.

    Standardization uses all training rows. The density of a query is a sum over the reference
    set R (all training rows, or a random subset of n_ref of them); the normalising maximum is
    taken over R (each reference point's density includes itself, as in Phase 4). The bandwidth is
    scikit-learn's estimate_bandwidth rule (mean distance to the int(n*quantile)-th neighbour,
    self included) on R or on a random subset of n_bw rows of R (scikit-learn's n_samples
    option). Returns (d_test, bandwidth, n_reference, n_bandwidth)."""
    rng = np.random.default_rng(seed)
    sc = StandardScaler().fit(X_train)
    A = sc.transform(X_train).astype(np.float64)
    B = sc.transform(X_test).astype(np.float64)
    ref = A if len(A) <= n_ref else A[np.sort(rng.choice(len(A), n_ref, replace=False))]
    bws = ref if len(ref) <= n_bw else ref[np.sort(rng.choice(len(ref), n_bw, replace=False))]
    h = estimate_bandwidth_nn(_sqdist(bws, bws), quantile)

    def dens(Q):
        out = np.empty(len(Q))
        for i in range(0, len(Q), chunk):
            out[i:i + chunk] = np.clip(1.0 - _sqdist(Q[i:i + chunk], ref) / h ** 2, 0, None).sum(axis=1)
        return out

    m = dens(ref).max()
    return 1.0 - dens(B) / m, float(h), int(len(ref)), int(len(bws))


# ----------------------------------------------------------------------------- hierarchical jackknife+
def hjplus_definition(mu_te, R, gidx, alpha):
    """Lee, Barber & Willett (2023, v4 2025), Appendix C.2.2, eq. (20), computed literally.

    mu_te[k, j] = prediction for test point j of the model fitted without group k (k = 0..K-1);
    R[i] = |Y_i - mu_{-k(i)}(X_i)| for calibration point i of group k(i) = gidx[i].
    Each calibration point of group k has weight 1/((K+1) N_k); a point mass 1/(K+1) sits at
    -inf (lower end) and at +inf (upper end).
      lower = Q'_alpha = sup{t : P(T <= t) <= alpha}  of  sum w_i delta(mu_{-k(i)}(x) - R_i) + delta(-inf)/(K+1)
      upper = Q_{1-alpha} = inf{t : P(T <= t) >= 1 - alpha}  of  sum w_i delta(mu_{-k(i)}(x) + R_i) + delta(+inf)/(K+1)
    Ties: P(T <= t) counts every atom at t. Comparisons use a 1e-12 tolerance so that exact rational
    boundaries (e.g. 1/(K+1) = alpha) resolve as in exact arithmetic."""
    mu_te = np.asarray(mu_te, float)
    R = np.asarray(R, float)
    gidx = np.asarray(gidx, int)
    K, m = mu_te.shape
    N = np.bincount(gidx, minlength=K)
    if (N == 0).any():
        raise ValueError("every group needs at least one calibration point")
    w = 1.0 / ((K + 1) * N[gidx])
    w_inf = 1.0 / (K + 1)
    lo, hi = np.empty(m), np.empty(m)
    for j in range(m):
        L = mu_te[gidx, j] - R
        v, inv = np.unique(L, return_inverse=True)
        F = w_inf + np.cumsum(np.bincount(inv, weights=w, minlength=v.size))
        if w_inf > alpha + TOL:
            lo[j] = -np.inf
        else:
            lo[j] = v[np.argmax(F > alpha + TOL)]
        U = mu_te[gidx, j] + R
        v, inv = np.unique(U, return_inverse=True)
        F = np.cumsum(np.bincount(inv, weights=w, minlength=v.size))
        hit = F >= 1 - alpha - TOL
        hi[j] = v[np.argmax(hit)] if hit.any() else np.inf
    return lo, hi


# ----------------------------------------------------------------------------- correlation
def spearman(x, y) -> float:
    rx, ry = rankdata(x), rankdata(y)
    if rx.std() == 0 or ry.std() == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def boot_spearman(x, y, B: int, seed: int, level: float = 0.95):
    """Spearman rho and a percentile bootstrap CI, resampling groups (pairs) with replacement."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    rng = np.random.default_rng(seed)
    n = x.size
    r = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        v = spearman(x[i], y[i])
        if np.isfinite(v):
            r.append(v)
    a = (1 - level) / 2
    return spearman(x, y), float(np.quantile(r, a)), float(np.quantile(r, 1 - a)), len(r)


def perm_spearman_mc(x, y, n_perm: int, seed: int):
    """Two-sided Monte Carlo permutation p-value for Spearman rho, (1 + count) / (n_perm + 1)."""
    rng = np.random.default_rng(seed)
    rx, ry = rankdata(x), rankdata(y)
    rx, ry = (rx - rx.mean()) / rx.std(), (ry - ry.mean()) / ry.std()
    rho = float(np.mean(rx * ry))
    cnt = 0
    for i in range(0, n_perm, 10000):
        b = min(10000, n_perm - i)
        P = np.argsort(rng.random((b, rx.size)), axis=1)
        cnt += int(np.sum(np.abs((ry[P] * rx).mean(axis=1)) >= abs(rho) - 1e-12))
    return rho, (1 + cnt) / (n_perm + 1)


# ----------------------------------------------------------------------------- H2 (descriptive)
def tail_enrichment(R, groups, q: float):
    """For calibration residuals R with group labels: threshold tau = q-quantile of R; per group,
    share of the tail (R > tau) it contributes, its share of calibration points, and their ratio."""
    R = np.asarray(R, float)
    g = np.asarray(groups)
    tau = float(np.quantile(R, q))
    tail = R > tau
    out = []
    n_tail = max(int(tail.sum()), 1)
    for u in np.unique(g):
        m = g == u
        st, sc = tail[m].sum() / n_tail, m.mean()
        out.append(dict(group=u, n_cal=int(m.sum()), share_cal=float(sc), share_tail=float(st),
                        enrichment=float(st / sc) if sc > 0 else np.nan, oof_mae=float(R[m].mean())))
    return out, tau


# ----------------------------------------------------------------------------- labels
def sg_tex(symbol: str) -> str:
    """Hermann-Mauguin symbol to mathtext: '-n' -> overbar n, '_n' -> subscript n."""
    s = re.sub(r"-(\d)", r"\\bar{\1}", str(symbol))
    s = re.sub(r"_(\d)", r"_{\1}", s)
    return f"${s}$"
