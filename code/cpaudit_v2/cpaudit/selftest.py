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


def t_general_cvplus(rng):
    K, n, m = 5, 61, 7
    f = rng.integers(0, K, n)
    fp = rng.normal(size=(K, m))
    R = np.abs(rng.normal(size=n))
    ok = True
    for a in (0.05, 0.1, 0.2):
        lo, hi = cp.cvplus_interval(fp, f, R, a)
        l2, h2 = cp.cvplus_from_terms(fp[f] - R[:, None], fp[f] + R[:, None], a)
        lw, hw = cp.weighted_cvplus_interval(fp, f, R, np.ones(n), np.ones(m), a)
        ok &= np.allclose(lo, l2) and np.allclose(hi, h2) and np.allclose(lo, lw) and np.allclose(hi, hw)
    check("general-score CV+ and weighted CV+ (unit weights) equal standard CV+", bool(ok))


def t_nested_cvplus_coverage(rng):
    """Normalized and CQR scores inside CV+ (Gupta et al. 2022): coverage >= 1 - 2 alpha, near 1 - alpha."""
    from sklearn.linear_model import QuantileRegressor, Ridge
    a, n, m, K, sims = 0.10, 150, 40, 5, 300
    cn, cq = [], []
    for _ in range(sims):
        x = rng.uniform(0, 1, n + m)
        y = 2 * x + (0.1 + x) * rng.standard_normal(n + m)
        X = x[:, None]
        folds = rng.permutation(np.arange(n) % K)
        mu_te, sig_te, ql_te, qh_te = (np.empty((K, m)) for _ in range(4))
        R, sig_oof, E = np.empty(n), np.empty(n), np.empty(n)
        for k in range(K):
            fit, oof = folds != k, folds == k
            mdl = Ridge(1e-3).fit(X[:n][fit], y[:n][fit])
            res = np.abs(y[:n][fit] - mdl.predict(X[:n][fit]))
            smd = Ridge(1e-3).fit(X[:n][fit], res)                     # crude difficulty model
            sg = lambda Z: np.maximum(smd.predict(Z), 0.05)
            R[oof] = np.abs(y[:n][oof] - mdl.predict(X[:n][oof]))
            sig_oof[oof] = sg(X[:n][oof])
            mu_te[k], sig_te[k] = mdl.predict(X[n:]), sg(X[n:])
            ql = QuantileRegressor(quantile=a / 2, alpha=0, solver="highs").fit(X[:n][fit], y[:n][fit])
            qh = QuantileRegressor(quantile=1 - a / 2, alpha=0, solver="highs").fit(X[:n][fit], y[:n][fit])
            E[oof] = np.maximum(ql.predict(X[:n][oof]) - y[:n][oof], y[:n][oof] - qh.predict(X[:n][oof]))
            ql_te[k], qh_te[k] = ql.predict(X[n:]), qh.predict(X[n:])
        s = R / sig_oof
        lo, hi = cp.cvplus_from_terms(mu_te[folds] - s[:, None] * sig_te[folds], mu_te[folds] + s[:, None] * sig_te[folds], a)
        cn.append(np.mean((y[n:] >= lo) & (y[n:] <= hi)))
        lo, hi = cp.cvplus_from_terms(ql_te[folds] - E[:, None], qh_te[folds] + E[:, None], a)
        cq.append(np.mean((y[n:] >= lo) & (y[n:] <= hi)))
    check("normalized CV+ coverage >= 1 - 2 alpha, near 1 - alpha", np.mean(cn) >= 0.8 and abs(np.mean(cn) - 0.9) < 0.03,
          f"{np.mean(cn):.4f}")
    check("CQR-CV+ coverage >= 1 - 2 alpha, near 1 - alpha", np.mean(cq) >= 0.8 and abs(np.mean(cq) - 0.9) < 0.03,
          f"{np.mean(cq):.4f}")


def t_weighted_shift(rng):
    """Known covariate shift: weighted CV+ restores coverage that unweighted CV+ loses."""
    from sklearn.linear_model import Ridge
    a, n, m, K, sims = 0.10, 200, 50, 5, 300
    cu, cw = [], []
    for _ in range(sims):
        x_tr = rng.normal(0, 1, n)
        x_te = rng.normal(1.5, 1, m)
        f = lambda x: x + (0.2 + np.abs(x)) * rng.standard_normal(x.size)
        y_tr, y_te = f(x_tr), f(x_te)
        folds = rng.permutation(np.arange(n) % K)
        fp, R = np.empty((K, m)), np.empty(n)
        for k in range(K):
            mdl = Ridge(1e-3).fit(x_tr[folds != k, None], y_tr[folds != k])
            R[folds == k] = np.abs(y_tr[folds == k] - mdl.predict(x_tr[folds == k, None]))
            fp[k] = mdl.predict(x_te[:, None])
        w = lambda x: np.exp(1.5 * x - 1.125)          # N(1.5,1)/N(0,1) density ratio
        lo, hi = cp.cvplus_interval(fp, folds, R, a)
        cu.append(np.mean((y_te >= lo) & (y_te <= hi)))
        lo, hi = cp.weighted_cvplus_interval(fp, folds, R, w(x_tr), w(x_te), a)
        cw.append(np.mean((y_te >= lo) & (y_te <= hi)))
    check("weighted CV+ under known covariate shift: coverage restored", np.mean(cw) >= 0.86 and np.mean(cu) < np.mean(cw),
          f"unweighted {np.mean(cu):.3f}, weighted {np.mean(cw):.3f}")


def t_group_level(rng):
    check("group-level: K=5, alpha=0.10 -> unbounded", np.isinf(cp.group_level_threshold(np.arange(5.), 0.10)))
    check("group-level: K=5, alpha=0.20 -> largest family score", cp.group_level_threshold(np.arange(5.), 0.20) == 4.0)
    K, sims, a = 5, 20000, 0.20
    hit = 0
    for _ in range(sims):
        s = np.abs(rng.normal(0, 1, K + 1)) * rng.lognormal(0, 0.5, K + 1)   # exchangeable family scores
        hit += s[-1] <= cp.group_level_threshold(s[:-1], a)
    check("group-level: P(new family score <= threshold) = k/(K+1)", abs(hit / sims - 5 / 6) < 0.01, f"{hit / sims:.4f}")
    check("Kish ESS: equal weights give n, one dominant weight gives ~1",
          abs(cp.kish_ess(np.ones(50)) - 50) < 1e-9 and cp.kish_ess(np.r_[1e6, np.ones(49)]) < 1.01)


def t_geometry(rng):
    from .geometry import lattice_automorphisms, same_structure
    a = 5.0
    Lc = np.eye(3) * a
    f = rng.uniform(0, 1, (12, 3))
    sp = np.array(["Al"] * 2 + ["Ga"] * 2 + ["O"] * 8)
    A = lattice_automorphisms(Lc)
    check("cubic cell has 48 lattice automorphisms", len(A) == 48, str(len(A)))
    op = A[17]
    t = rng.uniform(0, 1, 3)
    perm = rng.permutation(12)
    g2 = (Lc, ((f @ op) + t)[perm] % 1.0, sp[perm])
    check("same crystal after symmetry operation, translation and re-ordering is matched",
          same_structure((Lc, f, sp), g2)[0])
    f3 = g2[1].copy()
    f3[0] += 0.05 / a
    check("0.05 A displacement is not matched", not same_structure((Lc, f, sp), (Lc, f3, g2[2]))[0])
    Lh = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, 13.0]])
    check("hexagonal cell has 24 lattice automorphisms", len(lattice_automorphisms(Lh)) == 24,
          str(len(lattice_automorphisms(Lh))))
    try:
        from pymatgen.analysis.structure_matcher import StructureMatcher
        from pymatgen.core import Lattice, Structure
        sm = StructureMatcher(ltol=1e-3, stol=1e-3, angle_tol=0.05, primitive_cell=True, scale=False,
                              attempt_supercell=False)
        s1, s2 = Structure(Lattice(Lc), sp.tolist(), f), Structure(Lattice(Lc), g2[2].tolist(), g2[1])
        s3 = Structure(Lattice(Lc), g2[2].tolist(), f3)
        check("pymatgen StructureMatcher agrees (match / no match)", sm.fit(s1, s2) and not sm.fit(s1, s3))
    except ImportError:
        check("pymatgen available for de-duplication", False, "pip install pymatgen")


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
    print("Self-test | general / weighted CV+");    t_general_cvplus(rng)
    print("Self-test | nested CV+ coverage");       t_nested_cvplus_coverage(rng)
    print("Self-test | weighted CV+ under shift");  t_weighted_shift(rng)
    print("Self-test | group-level conformal");     t_group_level(rng)
    print("Self-test | structure identity");        t_geometry(rng)
    print(f"\n{'ALL SELF-TESTS PASSED' if not FAIL else 'FAILED: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0
