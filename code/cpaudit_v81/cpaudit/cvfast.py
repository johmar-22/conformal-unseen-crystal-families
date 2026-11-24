"""Exact CV+ intervals for large calibration sets, with optional group weights, and the group-level
split methods used for H3. Every function has a brute-force twin used by the self-tests and by a
runtime spot check inside each unit.

CV+ for a general nested score (Barber et al. 2021; Gupta, Kuchibhotla & Ramdas 2022). Calibration
point i lies in inner fold k(i) and has score S_i. For test point x and the model of fold k:
    lower term  L_i(x) = c_lo[k(i)](x) - s[k(i)](x) * S_i
    upper term  U_i(x) = c_hi[k(i)](x) + s[k(i)](x) * S_i         (s > 0)
absolute residual: c_lo = c_hi = mu_k(x), s = 1, S = |y - mu_k(x_i)|
normalized       : c = mu_k(x), s = sigma_k(x), S = |y - mu_k| / sigma_k(x_i)
CQR              : c_lo = q_lo,k(x), c_hi = q_hi,k(x), s = 1, S = E_i (may be negative)

Unweighted (Barber et al. 2021): lower = floor(alpha (n+1))-th smallest L_i, upper =
ceil((1-alpha)(n+1))-th smallest U_i (+-inf when the index is out of range).

Group weighted: calibration point i of group g gets weight 1 / ((G + 1) n_g) and an atom of mass
1 / (G + 1) sits at -inf (lower) and +inf (upper), the weighting of Lee, Barber & Willett
(arXiv:2306.06342 v4, eq. 6 and eq. 20). lower = sup{t: P(L <= t) <= alpha}, i.e. the smallest
L_i at which the atom plus the weight of all L_j <= L_i exceeds alpha; upper = smallest U_i whose
cumulative weight reaches 1 - alpha. Lee et al. prove this weighting for split conformal
(Theorem 1) and leave-one-group-out jackknife+ (Theorem 5); the K-fold version is our extension.
"""
from __future__ import annotations

import math

import numpy as np

from .conformal import k_lower, k_upper

TOL = 1e-12


def _as2d(a, K, n):
    a = np.asarray(a, float)
    return np.broadcast_to(a, (K, n)) if a.ndim < 2 else a


class _Folds:
    """Calibration scores sorted within each inner fold, with cumulative weights."""

    def __init__(self, S, fold, K, w=None):
        S = np.asarray(S, float)
        fold = np.asarray(fold)
        self.K = K
        self.S, self.cw, self.tot = [], [], []
        for k in range(K):
            m = fold == k
            o = np.argsort(S[m], kind="mergesort")
            self.S.append(S[m][o])
            if w is not None:
                wk = np.asarray(w, float)[m][o]
                self.cw.append(np.concatenate([[0.0], np.cumsum(wk)]))
                self.tot.append(float(wk.sum()))
        self.n = sum(len(s) for s in self.S)


def _idx_upper(Sk, c, s, t):
    """Exact number of upper terms c + s*S_i <= t (terms evaluated exactly as in the definition).
    searchsorted gives the index up to rounding; it is then corrected with the exact terms, which
    are non-decreasing in i because S is sorted and s > 0."""
    n = len(Sk)
    idx = np.searchsorted(Sk, (t - c) / s, side="right")
    for _ in range(n + 1):
        i1 = np.minimum(idx, n - 1)
        up = (idx < n) & (c + s * Sk[i1] <= t)
        i0 = np.maximum(idx - 1, 0)
        down = (idx > 0) & (c + s * Sk[i0] > t)
        if not (up.any() or down.any()):
            break
        idx = idx + up - down
    return idx


def _idx_lower(Sk, c, s, t):
    """Exact p = number of lower terms c - s*S_i > t (they form a prefix, terms non-increasing)."""
    n = len(Sk)
    p = np.searchsorted(Sk, (c - t) / s, side="left")
    for _ in range(n + 1):
        i1 = np.minimum(p, n - 1)
        up = (p < n) & (c - s * Sk[i1] > t)
        i0 = np.maximum(p - 1, 0)
        down = (p > 0) & (c - s * Sk[i0] <= t)
        if not (up.any() or down.any()):
            break
        p = p + up - down
    return p


def _count(F, c, s, t, side, weighted):
    """Number (or weight) of terms <= t, per test point."""
    out = np.zeros(t.shape)
    for k in range(F.K):
        nk = len(F.S[k])
        if nk == 0:
            continue
        if side == "upper":
            idx = _idx_upper(F.S[k], c[k], s[k], t)
            out += F.cw[k][idx] if weighted else idx
        else:
            p = _idx_lower(F.S[k], c[k], s[k], t)
            out += (F.tot[k] - F.cw[k][p]) if weighted else (nk - p)
    return out


def _smallest_term_above(F, c, s, t, side):
    """Smallest term strictly above t, over all folds (exact term values)."""
    best = np.full(t.shape, np.inf)
    for k in range(F.K):
        Sk = F.S[k]
        nk = len(Sk)
        if nk == 0:
            continue
        if side == "upper":
            idx = _idx_upper(Sk, c[k], s[k], t)
            val = np.where(idx < nk, c[k] + s[k] * Sk[np.minimum(idx, nk - 1)], np.inf)
        else:
            p = _idx_lower(Sk, c[k], s[k], t)
            val = np.where(p > 0, c[k] - s[k] * Sk[np.maximum(p - 1, 0)], np.inf)
        best = np.minimum(best, val)
    return best


def _select(F, c, s, side, target, weighted, strict, n_iter=200):
    """Smallest term t (exact) with count(t) >= target (or > target if strict), per test point."""
    allS = [Sk for Sk in F.S if len(Sk)]
    smin, smax = min(a[0] for a in allS), max(a[-1] for a in allS)
    if side == "upper":
        lo = (c + s * smin).min(axis=0) - 1.0 - np.abs((c + s * smin).min(axis=0))
        hi = (c + s * smax).max(axis=0) + 1.0
    else:
        lo = (c - s * smax).min(axis=0) - 1.0 - np.abs((c - s * smax).min(axis=0))
        hi = (c - s * smin).max(axis=0) + 1.0
    if strict:
        ok = lambda v: v > target + TOL
    else:
        ok = lambda v: v >= target - (TOL if weighted else 0)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        done = (mid <= lo) | (mid >= hi)
        if done.all():
            break
        good = ok(_count(F, c, s, mid, side, weighted))
        hi = np.where(good & ~done, mid, hi)
        lo = np.where(~good & ~done, mid, lo)
    ans = _smallest_term_above(F, c, s, lo, side)
    for _ in range(64):          # if the bisection stopped early, walk to the next term
        fail = ~ok(_count(F, c, s, ans, side, weighted))
        if not fail.any():
            return ans
        lo = np.where(fail, ans, lo)
        ans = np.where(fail, _smallest_term_above(F, c, s, lo, side), ans)
    raise RuntimeError("CV+ selection did not converge")


def cvplus(c_lo, c_hi, scale, S, fold, K, alpha, group=None):
    """CV+ interval (lo, hi) for every test point; group=None gives the unweighted version."""
    S = np.asarray(S, float)
    n_te = np.asarray(c_lo).shape[-1]
    c_lo, c_hi = _as2d(c_lo, K, n_te), _as2d(c_hi, K, n_te)
    sc = _as2d(1.0 if scale is None else scale, K, n_te)
    n = S.size
    if group is None:
        F = _Folds(S, fold, K)
        kl, kh = k_lower(n, alpha), k_upper(n, alpha)
        lo = np.full(n_te, -np.inf) if kl < 1 else _select(F, c_lo, sc, "lower", kl, False, False)
        hi = np.full(n_te, np.inf) if kh > n else _select(F, c_hi, sc, "upper", kh, False, False)
        return lo, hi
    g = np.asarray(group)
    ug, inv, cnt = np.unique(g, return_inverse=True, return_counts=True)
    G = ug.size
    w = 1.0 / ((G + 1) * cnt[inv])
    w_inf = 1.0 / (G + 1)
    F = _Folds(S, fold, K, w)
    lo = np.full(n_te, -np.inf) if w_inf > alpha + TOL else _select(F, c_lo, sc, "lower", alpha - w_inf, True, True)
    hi = np.full(n_te, np.inf) if 1.0 - w_inf < 1.0 - alpha - TOL else _select(F, c_hi, sc, "upper", 1.0 - alpha, True, False)
    return lo, hi


# ----------------------------------------------------------------------------- brute force twins
def cvplus_brute(c_lo, c_hi, scale, S, fold, K, alpha, group=None, cols=None):
    """Direct definition; builds all n x n_test terms (use on small problems or a few columns)."""
    S = np.asarray(S, float)
    fold = np.asarray(fold)
    n_te = np.asarray(c_lo).shape[-1]
    c_lo, c_hi = _as2d(c_lo, K, n_te), _as2d(c_hi, K, n_te)
    sc = _as2d(1.0 if scale is None else scale, K, n_te)
    cols = np.arange(n_te) if cols is None else np.asarray(cols)
    # select the test columns first: c[fold] for all columns would be n x n_test (GBs at JARVIS size)
    L = c_lo[:, cols][fold] - sc[:, cols][fold] * S[:, None]
    U = c_hi[:, cols][fold] + sc[:, cols][fold] * S[:, None]
    n = S.size
    if group is None:
        kl, kh = k_lower(n, alpha), k_upper(n, alpha)
        lo = np.full(len(cols), -np.inf) if kl < 1 else np.partition(L, kl - 1, axis=0)[kl - 1]
        hi = np.full(len(cols), np.inf) if kh > n else np.partition(U, kh - 1, axis=0)[kh - 1]
        return lo, hi
    ug, inv, cnt = np.unique(np.asarray(group), return_inverse=True, return_counts=True)
    G = ug.size
    w = 1.0 / ((G + 1) * cnt[inv])
    w_inf = 1.0 / (G + 1)
    lo, hi = np.empty(len(cols)), np.empty(len(cols))
    for j in range(len(cols)):
        o = np.argsort(L[:, j], kind="mergesort")
        cum = w_inf + np.cumsum(w[o])
        i = np.argmax(cum > alpha + TOL) if (cum > alpha + TOL).any() else None
        lo[j] = -np.inf if w_inf > alpha + TOL else L[o[i], j]
        o = np.argsort(U[:, j], kind="mergesort")
        cum = np.cumsum(w[o])
        hit = cum >= 1 - alpha - TOL
        hi[j] = U[o[np.argmax(hit)], j] if hit.any() else np.inf
    return lo, hi


def spot_check(c_lo, c_hi, scale, S, fold, K, alpha, lo, hi, group=None, n_check=48, seed=0):
    """Compare the fast result with the brute-force definition on a random subset of test points."""
    n_te = np.asarray(lo).size
    rng = np.random.default_rng(seed)
    cols = rng.choice(n_te, size=min(n_check, n_te), replace=False)
    blo, bhi = cvplus_brute(c_lo, c_hi, scale, S, fold, K, alpha, group=group, cols=cols)
    for a, b in ((np.asarray(lo)[cols], blo), (np.asarray(hi)[cols], bhi)):
        if not np.array_equal(np.isposinf(a), np.isposinf(b)) or not np.array_equal(np.isneginf(a), np.isneginf(b)):
            raise AssertionError("fast CV+ differs from the definition: different unbounded ends")
        f = np.isfinite(a)
        err = float(np.max(np.abs(a[f] - b[f]) / (1 + np.abs(b[f])))) if f.any() else 0.0
        if err > 1e-9:
            raise AssertionError(f"fast CV+ differs from the definition (max rel err {err:.2e})")
    return True


# ----------------------------------------------------------------------------- split / group level
def hcp_threshold(scores, groups, alpha):
    """Hierarchical split conformal threshold (Lee et al. eq. 6): Q_{1-alpha} of the distribution that
    puts 1/((K1+1) N_k) on each calibration score of group k and 1/(K1+1) on +inf."""
    s = np.asarray(scores, float)
    ug, inv, cnt = np.unique(np.asarray(groups), return_inverse=True, return_counts=True)
    K1 = ug.size
    w = 1.0 / ((K1 + 1) * cnt[inv])
    o = np.argsort(s, kind="mergesort")
    cum = np.cumsum(w[o])
    hit = cum >= 1 - alpha - TOL
    return float(s[o][np.argmax(hit)]) if hit.any() else np.inf


def single_subsample_threshold(scores, groups, alpha, rng):
    """One random calibration score per group (Dunn et al. 2023 single subsampling), then the
    ceil((K+1)(1-alpha))-th smallest of the K drawn scores."""
    s = np.asarray(scores, float)
    g = np.asarray(groups)
    drawn = []
    for u in np.unique(g):
        idx = np.flatnonzero(g == u)
        drawn.append(s[rng.choice(idx)])
    drawn = np.sort(drawn)
    k = k_upper(drawn.size, alpha)
    return float(drawn[k - 1]) if k <= drawn.size else np.inf


def interval_score(y, lo, hi, alpha):
    """Gneiting & Raftery (2007) interval score; inf for an unbounded interval."""
    y, lo, hi = (np.asarray(v, float) for v in (y, lo, hi))
    with np.errstate(invalid="ignore"):
        s = (hi - lo) + (2 / alpha) * np.maximum(lo - y, 0) + (2 / alpha) * np.maximum(y - hi, 0)
    return s
