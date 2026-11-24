"""Few-shot, family-specific conformal calibration (Phase 5): methods, sampling schemes and the
exact reference law. Everything here works on residuals z = y - mu(x) of ONE structure family,
where mu is a fixed model trained without that family. All three intervals have the form
mu(x) + Z with a set Z that does not depend on x, so a draw is summarised by Z alone.

Methods (alpha = miscoverage, m = number of labelled structures of the family)
---------------------------------------------------------------------------------
split ............ split conformal with |z| scores (Vovk et al. 2005; Lei et al., JASA 113,
                   1094, 2018): Z = [-q, q], q = k-th smallest |z_i|, k = ceil((1 - alpha)(m + 1)).
offset_full ...... full conformal for the location model y = mu(x) + b + e, with b fitted by the
                   mean residual of the augmented data (Vovk, Gammerman & Shafer 2005, Ch. 2;
                   Lei et al. 2018). Exact marginal coverage 1 - floor(alpha(m+1))/(m+1) under
                   exchangeability, like split conformal, but the family offset is removed.
                   Closed form (derived in the docstring of `offset_full_set`).
offset_jackplus .. jackknife+ (Barber, Candes, Ramdas & Tibshirani, Ann. Stat. 49, 486, 2021)
                   with the leave-one-out mean offset b_-i. For this location model it is
                   IDENTICAL to offset_full: its per-point intervals [b_-i - R_i, b_-i + R_i] are
                   exactly the full-conformal intervals between z_i and c_i, and every one of them
                   contains the mean residual, so the full-conformal set is the single interval
                   [t-th smallest lower end, t-th largest upper end], t = floor(alpha (m + 1)),
                   which is the jackknife+ interval. Hence the offset interval has the exact
                   full-conformal guarantee, not only the generic jackknife+ 1 - 2 alpha. Kept to
                   verify the identity in the self-test; not run separately in Phase 5.

Exact reference
---------------
For split conformal with calibration drawn without replacement from a finite family of N
structures and coverage measured on the other N - m, the coverage is (R - k)/(N - m), where the
population rank R of the threshold has the negative-hypergeometric law
    P(R = r) = C(r-1, k-1) C(N-r, m-k) / C(N, m),  r = k .. N - m + k   (no ties).
Its mean is exactly k/(m+1); as N -> inf it tends to Beta(k, m + 1 - k) (Vovk 2012, Angelopoulos
& Bates 2023, Sec. 3.2).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import gammaln
from scipy.stats import beta as beta_dist

from .conformal import k_lower, k_upper

METHODS = ("split", "offset_full")
INF = float("inf")


# ----------------------------------------------------------------------------- methods
def split_set(z_cal: np.ndarray, alpha: float):
    """[(-q, q)] with q the k-th smallest |z|; [(-inf, inf)] if k > m."""
    m = len(z_cal)
    k = k_upper(m, alpha)
    if k > m:
        return [(-INF, INF)]
    q = float(np.partition(np.abs(z_cal), k - 1)[k - 1])
    return [(-q, q)]


def offset_full_set(z_cal: np.ndarray, alpha: float):
    """Full conformal prediction set for the residual z of a new structure under y = mu + b + e.

    With candidate z, b(z) = (S + z)/(m + 1), S = sum z_i. Scores s_i = |z_i - b(z)|,
    s_0 = |z - b(z)|. The p-value is p(z) = (1 + N(z))/(m + 1), N(z) = #{i : s_i >= s_0}, and z
    is kept iff p(z) > alpha. Since s_i^2 - s_0^2 = (z_i - z)(z_i + z - 2 b(z)) and
    z_i + z - 2 b(z) = (m - 1)/(m + 1) * (z - c_i) with c_i = (2 S - (m + 1) z_i)/(m - 1),
    the event s_i >= s_0 is exactly z in the closed interval between z_i and c_i (m >= 2).
    The set {z : N(z) >= t}, t = smallest integer > alpha (m + 1) - 1, is found by a sweep.
    Returns a list of disjoint closed intervals (lo, hi)."""
    z = np.asarray(z_cal, dtype=float)
    m = z.size
    if m < 2:
        return [(-INF, INF)]
    t = int(math.floor(alpha * (m + 1) - 1 + 1e-9)) + 1
    if t <= 0:
        return [(-INF, INF)]
    S = z.sum()
    c = (2.0 * S - (m + 1) * z) / (m - 1)
    lo, hi = np.minimum(z, c), np.maximum(z, c)
    # sweep: +1 at lo, -1 just after hi (closed intervals: starts before ends at equal coordinates)
    pts = np.concatenate([lo, hi])
    kind = np.concatenate([np.zeros(m, int), np.ones(m, int)])  # 0 = start, 1 = end
    order = np.lexsort((kind, pts))
    count, out, cur = 0, [], None
    for i in order:
        x = pts[i]
        if kind[i] == 0:
            count += 1
            if count == t and cur is None:
                cur = x
        else:
            if count == t and cur is not None:
                out.append((cur, x))
                cur = None
            count -= 1
    merged = []
    for a, b in out:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged if merged else [(float("nan"), float("nan"))]   # empty set (not expected)


def offset_jackplus_set(z_cal: np.ndarray, alpha: float):
    """Jackknife+ with a leave-one-out mean offset: Z = [q_lo{b_-i - R_i}, q_hi{b_-i + R_i}],
    b_-i = mean of the other residuals, R_i = |z_i - b_-i|; q_lo = floor(alpha(m+1))-th smallest,
    q_hi = ceil((1-alpha)(m+1))-th smallest (Barber et al. 2021, eq. 6)."""
    z = np.asarray(z_cal, dtype=float)
    m = z.size
    if m < 2:
        return [(-INF, INF)]
    b = (z.sum() - z) / (m - 1)
    R = np.abs(z - b)
    kl, ku = k_lower(m, alpha), k_upper(m, alpha)
    lo = -INF if kl < 1 else float(np.sort(b - R)[kl - 1])
    hi = INF if ku > m else float(np.sort(b + R)[ku - 1])
    return [(lo, hi)]


SET_FUNCS = {"split": split_set, "offset_full": offset_full_set, "offset_jackplus": offset_jackplus_set}


def in_set(z: np.ndarray, intervals) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    ok = np.zeros(z.shape, bool)
    for a, b in intervals:
        ok |= (z >= a) & (z <= b)
    return ok


def hull_width(intervals) -> float:
    return float(intervals[-1][1] - intervals[0][0])


# ----------------------------------------------------------------------------- exact reference
def _log_comb(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def exact_split_law(N: int, m: int, alpha: float):
    """Exact law of split-conformal coverage on the remaining N - m structures when m are drawn
    at random without replacement (no ties). Returns (coverage values, probabilities)."""
    k = k_upper(m, alpha)
    if k > m:
        return np.array([1.0]), np.array([1.0])
    r = np.arange(k, N - m + k + 1)
    logp = _log_comb(r - 1, k - 1) + _log_comb(N - r, m - k) - _log_comb(N, m)
    p = np.exp(logp)
    p /= p.sum()
    return (r - k) / (N - m), p


def exact_summary(N: int, m: int, alpha: float, floor: float):
    v, p = exact_split_law(N, m, alpha)
    mean = float(np.sum(v * p))
    sd = float(np.sqrt(np.sum(p * (v - mean) ** 2)))
    return dict(exact_mean=mean, exact_sd=sd, exact_p_below=float(p[v < floor].sum()))


def beta_reference(m: int, alpha: float, floor: float):
    k = k_upper(m, alpha)
    if k > m:
        return dict(k=k, beta_mean=1.0, beta_sd=0.0, beta_p_below=0.0)
    b = beta_dist(k, m + 1 - k)
    return dict(k=k, beta_mean=float(b.mean()), beta_sd=float(b.std()), beta_p_below=float(b.cdf(floor)))


def nominal_mean(m: int, alpha: float) -> float:
    """Exact marginal coverage of split and offset_full: 1 - floor(alpha (m+1))/(m+1) = k/(m+1)."""
    return min(1.0, k_upper(m, alpha) / (m + 1))


# ----------------------------------------------------------------------------- sampling schemes
class Family:
    """Residuals of one structure family with composition groups and (Fd-3m) sub-families."""

    def __init__(self, z, groups, sub=None, ids=None, y=None):
        self.z = np.asarray(z, float)
        self.groups = np.asarray(groups)
        self.sub = None if sub is None else np.asarray(sub)
        self.ids = ids
        self.y = y
        self.ug, self.ginv = np.unique(self.groups, return_inverse=True)
        self.members = [np.where(self.ginv == g)[0] for g in range(len(self.ug))]

    @property
    def N(self):
        return self.z.size

    @property
    def G(self):
        return len(self.ug)


def draw_random(fam: Family, m: int, rng):
    p = rng.permutation(fam.N)
    return p[:m], p[m:]


def draw_clustered(fam: Family, m: int, rng, min_test_groups: int):
    """Whole composition groups in random order until m structures (the last group subsampled);
    test = structures of groups not touched. None if fewer than min_test_groups remain."""
    order = rng.permutation(fam.G)
    cal, used = [], 0
    for j, g in enumerate(order):
        mem = fam.members[g]
        need = m - len(cal)
        if len(mem) <= need:
            cal.extend(mem.tolist())
        else:
            cal.extend(rng.choice(mem, need, replace=False).tolist())
        used = j + 1
        if len(cal) >= m:
            break
    if len(cal) < m or fam.G - used < min_test_groups:
        return None
    test_groups = order[used:]
    test = np.concatenate([fam.members[g] for g in test_groups])
    return np.asarray(cal), test, test_groups


def draw_one_per_group(fam: Family, m: int, rng, min_test_groups: int):
    """m groups at random, one random structure from each (single subsampling, Dunn et al. 2023);
    test = all structures of the other groups. None if fewer than min_test_groups remain."""
    if fam.G - m < min_test_groups:
        return None
    order = rng.permutation(fam.G)
    cal = np.array([rng.choice(fam.members[g]) for g in order[:m]])
    test_groups = order[m:]
    test = np.concatenate([fam.members[g] for g in test_groups])
    return cal, test, test_groups


def draw_subfamily(fam: Family, m: int, rng, src: str, dst: str):
    a = np.where(fam.sub == src)[0]
    b = np.where(fam.sub == dst)[0]
    if len(a) < m + 1 or len(b) == 0:
        return None
    return rng.choice(a, m, replace=False), b


def group_weighted(covered: np.ndarray, test_idx: np.ndarray, fam: Family, test_groups) -> float:
    """Mean over test groups of the within-group covered fraction."""
    pos = {int(i): j for j, i in enumerate(test_idx)}
    vals = [covered[[pos[int(i)] for i in fam.members[g]]].mean() for g in test_groups]
    return float(np.mean(vals))
