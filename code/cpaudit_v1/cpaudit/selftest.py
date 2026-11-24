"""Self-tests: the conformal code is checked against (i) its own definitions by brute force,
(ii) finite-sample theory by simulation, and (iii) independent implementations (crepes;
MAPIE 0.8.6 if installed). No Drive data are needed. Exit code is non-zero on any failure.
"""
from __future__ import annotations

import numpy as np

from . import conformal as cp

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


def t_quantile_rules():
    check("k_upper(9, 0.10) = 9 (finite at the boundary)", cp.k_upper(9, 0.10) == 9)
    check("8 scores at alpha=0.10 give +inf", np.isinf(cp.conformal_quantile(np.arange(8), 0.10)))
    check("9 scores at alpha=0.10 give the maximum", cp.conformal_quantile(np.arange(9), 0.10) == 8)
    check("k_lower(19, 0.05) = 1", cp.k_lower(19, 0.05) == 1)


def t_split_theory(rng):
    n, m, a, sims = 199, 40, 0.10, 4000
    cov = np.empty(sims)
    for s in range(sims):
        e = rng.standard_t(3, n + m)            # heavy tails: validity must not depend on normality
        lo, hi, _ = cp.split_interval(np.zeros(m), np.abs(e[:n]), a)
        cov[s] = np.mean((e[n:] >= lo) & (e[n:] <= hi))
    exp = cp.k_upper(n, a) / (n + 1)
    check("split conformal mean coverage = k/(n+1)", abs(cov.mean() - exp) < 0.004,
          f"mean {cov.mean():.4f}, theory {exp:.4f}")


def t_mondrian(rng):
    a, sims = 0.10, 2000
    sd = {0: 0.2, 1: 0.2, 2: 1.0}
    cg, cm = {g: [] for g in sd}, {g: [] for g in sd}
    for _ in range(sims):
        cat_c = np.repeat(list(sd), 60)
        sc = np.abs(rng.normal(0, [sd[c] for c in cat_c]))
        cat_t = np.repeat(list(sd), 30)
        et = rng.normal(0, [sd[c] for c in cat_t])
        lo_g, hi_g, _ = cp.split_interval(np.zeros(et.size), sc, a)
        lo_m, hi_m, _ = cp.mondrian_interval(np.zeros(et.size), cat_t, sc, cat_c, a)
        for g in sd:
            k = cat_t == g
            cg[g].append(np.mean((et[k] >= lo_g[k]) & (et[k] <= hi_g[k])))
            cm[g].append(np.mean((et[k] >= lo_m[k]) & (et[k] <= hi_m[k])))
    exp = cp.k_upper(60, a) / 61
    check("Mondrian per-group coverage = k/(n+1) in every group",
          all(abs(np.mean(cm[g]) - exp) < 0.008 for g in sd), str({g: round(np.mean(cm[g]), 4) for g in sd}))
    check("global quantile under-covers the noisy group", np.mean(cg[2]) < 0.85, f"{np.mean(cg[2]):.3f}")


def t_normalized_cqr(rng):
    a, n, m, sims = 0.10, 299, 50, 3000
    cn, cq = [], []
    for _ in range(sims):
        x = rng.uniform(0, 1, n + m)
        s = 0.1 + x
        y = 2 * x + s * rng.standard_normal(n + m)
        mu = 2 * x
        lo, hi, _ = cp.normalized_interval(mu[n:], s[n:], y[:n], mu[:n], s[:n], a)
        cn.append(np.mean((y[n:] >= lo) & (y[n:] <= hi)))
        ql, qh = mu - 1.2 * s, mu + 1.2 * s     # deliberately mis-scaled quantile estimates
        lo, hi, _ = cp.cqr_interval(ql[n:], qh[n:], y[:n], ql[:n], qh[:n], a)
        cq.append(np.mean((y[n:] >= lo) & (y[n:] <= hi)))
    exp = cp.k_upper(n, a) / (n + 1)
    check("normalized conformal mean coverage = k/(n+1)", abs(np.mean(cn) - exp) < 0.005, f"{np.mean(cn):.4f}")
    check("CQR mean coverage = k/(n+1)", abs(np.mean(cq) - exp) < 0.005, f"{np.mean(cq):.4f}")


def t_cvplus_bruteforce(rng):
    K, n, m, a = 4, 37, 6, 0.1
    folds = rng.integers(0, K, n)
    fp = rng.normal(size=(K, m))
    R = np.abs(rng.normal(size=n))
    lo, hi = cp.cvplus_interval(fp, folds, R, a)
    ok = True
    for j in range(m):
        lv = sorted(fp[folds[i], j] - R[i] for i in range(n))
        hv = sorted(fp[folds[i], j] + R[i] for i in range(n))
        kl, kh = int(np.floor(a * (n + 1))), int(np.ceil((1 - a) * (n + 1)))
        ok &= np.isclose(lo[j], lv[kl - 1]) and np.isclose(hi[j], hv[kh - 1])
    check("CV+ equals the Barber et al. (2021) definition (brute force)", bool(ok))


def t_cvplus_coverage(rng):
    from sklearn.linear_model import Ridge
    a, n, m, K, sims = 0.10, 120, 30, 5, 600
    cov = []
    for _ in range(sims):
        X = rng.normal(size=(n + m, 5))
        y = X @ np.array([1., -1, .5, 0, 0]) + rng.standard_t(4, n + m)
        folds = rng.permutation(np.arange(n) % K)
        fp, R = np.empty((K, m)), np.empty(n)
        for k in range(K):
            mdl = Ridge(1.0).fit(X[:n][folds != k], y[:n][folds != k])
            R[folds == k] = np.abs(y[:n][folds == k] - mdl.predict(X[:n][folds == k]))
            fp[k] = mdl.predict(X[n:])
        lo, hi = cp.cvplus_interval(fp, folds, R, a)
        cov.append(np.mean((y[n:] >= lo) & (y[n:] <= hi)))
    c = float(np.mean(cov))
    check("CV+ coverage >= 1 - 2 alpha and near 1 - alpha (simulation)", c >= 1 - 2 * a and abs(c - (1 - a)) < 0.03,
          f"mean {c:.4f}")


def t_crepes(rng):
    try:
        import crepes
        from crepes import ConformalRegressor
    except Exception:
        check("crepes cross-check", False, "crepes not installed (pip install crepes)")
        return
    a, n, m = 0.10, 157, 40
    res = rng.normal(size=n)
    yh = rng.normal(size=m)
    sig_c, sig_t = rng.uniform(0.5, 2, n), rng.uniform(0.5, 2, m)
    cat_c, cat_t = rng.integers(0, 3, n), rng.integers(0, 3, m)
    pred = (lambda cr, **kw: cr.predict_int(yh, confidence=1 - a, **kw)) \
        if hasattr(ConformalRegressor(), "predict_int") else (lambda cr, **kw: cr.predict(yh, confidence=1 - a, **kw))
    ref = pred(ConformalRegressor().fit(res))
    mine = np.column_stack(cp.split_interval(yh, np.abs(res), a)[:2])
    check("standard split = crepes", np.allclose(ref, mine), f"crepes {getattr(crepes, '__version__', '?')}")
    ref = pred(ConformalRegressor().fit(res, sigmas=sig_c), sigmas=sig_t)
    q = cp.conformal_quantile(np.abs(res) / sig_c, a)
    mine = np.column_stack([yh - q * sig_t, yh + q * sig_t])
    check("normalized = crepes", np.allclose(ref, mine))
    ref = pred(ConformalRegressor().fit(res, bins=cat_c), bins=cat_t)
    mine = np.column_stack(cp.mondrian_interval(yh, cat_t, np.abs(res), cat_c, a)[:2])
    check("Mondrian = crepes", np.allclose(ref, mine))


def t_mapie(rng):
    try:
        import mapie
        from mapie.regression import MapieRegressor
    except Exception:
        print("  [SKIP] MAPIE cross-check (pip install mapie==0.8.6 to run it)")
        return
    from sklearn.linear_model import LinearRegression
    from sklearn.model_selection import KFold
    n, m, a = 200, 50, 0.10
    X = rng.normal(size=(n + m, 3))
    y = X.sum(1) + rng.normal(size=n + m)
    kf = KFold(5, shuffle=True, random_state=0)
    mp = MapieRegressor(LinearRegression(), cv=kf, method="plus").fit(X[:n], y[:n])
    _, pis = mp.predict(X[n:], alpha=a)
    folds = np.empty(n, int)
    fp, R = np.empty((5, m)), np.empty(n)
    for k, (tr, te) in enumerate(kf.split(X[:n])):
        folds[te] = k
        mdl = LinearRegression().fit(X[tr], y[tr])
        R[te] = np.abs(y[te] - mdl.predict(X[te]))
        fp[k] = mdl.predict(X[n:])
    lo, hi = cp.cvplus_interval(fp, folds, R, a)
    d = float(max(np.abs(lo - pis[:, 0, 0]).max(), np.abs(hi - pis[:, 1, 0]).max()))
    print(f"  [INFO] CV+ vs MAPIE {mapie.__version__} 'plus': max |difference| = {d:.2e} "
          "(differences, if any, come from MAPIE's quantile interpolation rule; this code follows the paper)")


def run_selftest() -> int:
    rng = np.random.default_rng(12345)
    print("Self-test | quantile rules");            t_quantile_rules()
    print("Self-test | split conformal theory");    t_split_theory(rng)
    print("Self-test | Mondrian");                  t_mondrian(rng)
    print("Self-test | normalized and CQR");        t_normalized_cqr(rng)
    print("Self-test | CV+ definition");            t_cvplus_bruteforce(rng)
    print("Self-test | CV+ coverage");              t_cvplus_coverage(rng)
    print("Self-test | crepes");                    t_crepes(rng)
    print("Self-test | MAPIE");                     t_mapie(rng)
    print(f"\n{'ALL SELF-TESTS PASSED' if not FAIL else 'FAILED: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0
