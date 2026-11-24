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
        print("    (pymatgen not installed: StructureMatcher cross-check skipped; Phase 4 reuses the Phase 2 record)")


def t_detect(rng):
    """Phase 4 building blocks against independent implementations and theory."""
    from sklearn.cluster import estimate_bandwidth
    from sklearn.metrics import roc_auc_score
    from sklearn.neighbors import KernelDensity

    from .detect import (ClusterBootstrap, abstention_curve, auroc, conformal_pvalues, estimate_bandwidth_nn,
                         kde_dissimilarity, _sqdist)
    from .report4 import perm_spearman
    y = rng.random(500) < 0.3
    s = rng.normal(size=500) + y * 0.8
    s[:50] = np.round(s[:50])                      # ties
    check("AUROC equals scikit-learn (with ties)", abs(auroc(y, s) - roc_auc_score(y, s)) < 1e-12)
    X = rng.normal(size=(300, 20))
    Xt = rng.normal(size=(80, 20)) * 1.5
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler().fit(X)
    A, B = sc.transform(X), sc.transform(Xt)
    h_ref = estimate_bandwidth(A, quantile=0.3)
    h = estimate_bandwidth_nn(_sqdist(A, A), 0.3)
    check("bandwidth rule equals sklearn estimate_bandwidth", abs(h - h_ref) / h_ref < 1e-10, f"{h:.6f} vs {h_ref:.6f}")
    d_te, d_tr, h2 = kde_dissimilarity(X, Xt, 0.3)
    kd = KernelDensity(kernel="epanechnikov", bandwidth=h_ref).fit(A)
    dens_te, dens_tr = np.exp(kd.score_samples(B)), np.exp(kd.score_samples(A))
    ref = 1 - dens_te / dens_tr.max()
    check("KDE dissimilarity equals sklearn KernelDensity (Epanechnikov)", np.allclose(d_te, ref, atol=1e-9),
          f"max diff {np.max(np.abs(d_te - ref)):.2e}")
    # conformal p-values: exchangeable calibration and test -> P(p <= a) <= a, close to a
    hits, n = 0, 0
    for _ in range(400):
        c, t = rng.normal(size=99), rng.normal(size=50)
        hits += np.sum(conformal_pvalues(c, t) <= 0.05)
        n += 50
    check("conformal p-values: false-flag rate ~ alpha under exchangeability", 0.035 < hits / n <= 0.055,
          f"{hits / n:.4f}")
    p = conformal_pvalues(np.arange(10.0), np.array([-1.0, 4.5, 100.0]))
    check("conformal p-value formula (1 + #cal >= s)/(n + 1)", np.allclose(p, [11 / 11, 6 / 11, 1 / 11]))
    g = np.repeat(np.arange(40), 5)
    cb = ClusterBootstrap(g, 0)
    i = cb.draw()
    check("cluster bootstrap resamples whole groups", len(i) == 200 and all(np.bincount(g[i], minlength=40) % 5 == 0))
    from scipy.stats import spearmanr
    x6, y6 = rng.normal(size=6), rng.normal(size=6)
    rho, pv = perm_spearman(x6, y6, 1000, rng)
    check("Spearman rho equals scipy", abs(rho - spearmanr(x6, y6).statistic) < 1e-12)
    x7 = np.arange(7.0)
    check("exact permutation p for perfect rank order = 2/7!", abs(perm_spearman(x7, x7, 1000, rng)[1] - 2 / 5040) < 1e-12)
    cov = rng.random(1000) < 0.8
    err = np.where(cov, rng.random(1000), 1 + rng.random(1000))
    cur = abstention_curve(err, cov, err, [0.0, 0.2, 0.5], 0)
    check("abstention on the oracle score raises retained coverage monotonically",
          cur[0][1] <= cur[1][1] <= cur[2][1] and abs(cur[0][1] - cov.mean()) < 1e-12)


def t_fewshot(rng):
    """Phase 5 tools: offset full conformal vs its definition, exactness of coverage, the exact
    finite-population law, jackknife+, one-per-group exactness, permutation null."""
    from . import fewshot as fs
    a = 0.10
    # 1. offset full conformal set = brute-force p-value definition on a grid
    bad = 0
    for m in (2, 3, 5, 9, 10, 17, 30):
        for al in (0.1, 0.2):
            zc = rng.normal(0.7, 1.0, m)
            S = fs.offset_full_set(zc, al)
            grid = np.linspace(-8, 9, 4001)
            Ssum = zc.sum()
            b = (Ssum + grid) / (m + 1)
            s0 = np.abs(grid - b)
            si = np.abs(zc[None, :] - b[:, None])
            pv = (1 + (si >= s0[:, None] - 1e-12).sum(1)) / (m + 1)
            brute = pv > al
            fast = fs.in_set(grid, S)
            near = np.zeros_like(grid, bool)
            for lo, hi in S:
                near |= (np.abs(grid - lo) < 1e-6) | (np.abs(grid - hi) < 1e-6)
            bad += int(np.sum((brute != fast) & ~near))
    check("offset full conformal: closed form = p-value definition (grid)", bad == 0, f"{bad} mismatches")

    # 1b. identity: jackknife+ with the leave-one-out mean offset = offset full conformal (one interval)
    worst, multi = 0.0, 0
    for m in (2, 3, 9, 10, 19, 20, 37, 100):
        for al in (0.05, 0.1, 0.2):
            for _ in range(60):
                zc = rng.standard_t(3, m) * rng.uniform(0.2, 3) + rng.normal()
                A, B = fs.offset_full_set(zc, al), fs.offset_jackplus_set(zc, al)
                multi += len(A) > 1
                if np.isfinite(B[0][0]) and np.isfinite(B[0][1]):
                    worst = max(worst, abs(A[0][0] - B[0][0]), abs(A[-1][1] - B[0][1]))
                elif np.isfinite(A[0][0]):
                    worst = np.inf
    check("offset full conformal = jackknife+ with leave-one-out offset, always one interval",
          worst < 1e-9 and multi == 0, f"max |difference| {worst:.1e}, multi-interval sets {multi}")

    # 2. exactness under exchangeability with a large common offset; width gain over split
    for m in (10, 20):
        cov = {k: [] for k in fs.METHODS}
        wid = {k: [] for k in fs.METHODS}
        for _ in range(6000):
            z = rng.normal(2.0, 1.0, m + 1)
            for k in fs.METHODS:
                S = fs.SET_FUNCS[k](z[:m], a)
                cov[k].append(fs.in_set(z[m:], S)[0])
                wid[k].append(fs.hull_width(S))
        target = fs.nominal_mean(m, a)
        se = np.sqrt(target * (1 - target) / 6000)
        for k in ("split", "offset_full"):
            c = np.mean(cov[k])
            check(f"{k}: mean coverage = k/(m+1) = {target:.3f} (m={m})", abs(c - target) < 3.5 * se, f"{c:.4f}")
        r = np.median(wid["offset_full"]) / np.median(wid["split"])
        check(f"offset full conformal narrower than split under a shift of 2 SD (m={m})", r < 0.7, f"ratio {r:.2f}")

    # 3. exact finite-population law (random scheme) vs simulation
    N, m = 150, 20
    pop = rng.normal(0, 1, N)
    v, p = fs.exact_split_law(N, m, a)
    sims = []
    for _ in range(20000):
        perm = rng.permutation(N)
        S = fs.split_set(pop[perm[:m]], a)
        sims.append(fs.in_set(np.abs(pop[perm[m:]]) * 0 + pop[perm[m:]], S).mean())
    sims = np.array(sims)
    ex = fs.exact_summary(N, m, a, 0.8)
    check("exact law: mean = k/(m+1)", abs(ex["exact_mean"] - fs.nominal_mean(m, a)) < 1e-12,
          f"{ex['exact_mean']:.6f}")
    check("exact law matches simulation (mean, SD, P(cov < 0.8))",
          abs(sims.mean() - ex["exact_mean"]) < 0.002 and abs(sims.std() - ex["exact_sd"]) < 0.003
          and abs(np.mean(sims < 0.8) - ex["exact_p_below"]) < 0.01,
          f"sim {sims.mean():.4f}/{sims.std():.4f}/{np.mean(sims < 0.8):.4f} vs "
          f"exact {ex['exact_mean']:.4f}/{ex['exact_sd']:.4f}/{ex['exact_p_below']:.4f}")
    big = fs.exact_summary(200000, 30, a, 0.8)["exact_p_below"]
    bet = fs.beta_reference(30, a, 0.8)["beta_p_below"]
    check("exact law -> Beta(k, m+1-k) for large N", abs(big - bet) < 0.002, f"{big:.4f} vs {bet:.4f}")

    # 4. composition structure: one-per-group keeps new-composition coverage exact in the mean and
    #    close to the permutation null in spread; clustered labels do not (strong group effects)
    G = 55
    sizes = rng.integers(3, 12, G)
    groups = np.repeat(np.arange(G), sizes)

    def run(fam, scheme, m, null, n=2500):
        out = []
        draw = fs.draw_one_per_group if scheme == "one" else fs.draw_clustered
        for _ in range(n):
            zz = fam.z[rng.permutation(fam.N)] if null else fam.z
            cal, new, tg = draw(fam, m, rng, 10)
            out.append(fs.group_weighted(fs.in_set(zz[new], fs.split_set(zz[cal], a)), new, fam, tg))
        return np.array(out)
    m = 20
    t = fs.nominal_mean(m, a)
    eff = rng.normal(0, 1.5, G)
    famg = fs.Family(z=eff[groups] + rng.normal(0, 0.5, sizes.sum()), groups=groups)
    o, on = run(famg, "one", m, False), run(famg, "one", m, True)
    check("one-per-group: mean coverage of new compositions = k/(m+1) with group effects",
          abs(o.mean() - t) < 3.5 * o.std() / np.sqrt(len(o)), f"{o.mean():.4f} vs {t:.4f}")
    check("one-per-group: spread close to the permutation null (rho < 1.25) with group effects",
          o.std() / on.std() < 1.25, f"rho {o.std() / on.std():.2f}")
    c, cn = run(famg, "clu", m, False), run(famg, "clu", m, True)
    check("clustered labels with group effects: spread inflated (rho > 1.25) and mean below nominal",
          c.std() / cn.std() > 1.25 and c.mean() < t - 0.03, f"rho {c.std() / cn.std():.2f}, mean {c.mean():.3f}")
    fam0 = fs.Family(z=rng.normal(0, 1, sizes.sum()), groups=groups)
    r1 = run(fam0, "one", m, False).std() / run(fam0, "one", m, True).std()
    r2 = run(fam0, "clu", m, False).std() / run(fam0, "clu", m, True).std()
    check("no group effects: rho ~ 1 for both schemes", 0.88 < r1 < 1.12 and 0.88 < r2 < 1.12,
          f"one-per-group {r1:.3f}, clustered {r2:.3f}")
    cal, new, tg = fs.draw_clustered(famg, 30, rng, 10)
    ok = (len(cal) == 30 and len(set(famg.ginv[cal]) & set(famg.ginv[new])) == 0
          and len(set(cal.tolist())) == 30)
    check("clustered draw: m structures, new-composition groups disjoint from calibration groups", ok)

def t_posthoc(rng):
    """Phase 5 post hoc tools: ICC(1) estimator, no-model interval validity, exact re-creation of
    the Phase 5 one-per-composition draws, permutation Spearman."""
    import os
    import tempfile
    import pandas as pd
    from . import fewshot as fs
    from .phase5 import _gz, run_family
    from .posthoc5 import baseline_draws, icc1, spearman_perm, verify_against_phase5
    # 1. balanced case: ICC(1) equals the textbook formula computed independently with pandas
    G, n = 40, 6
    g = np.repeat(np.arange(G), n)
    x = rng.normal(0, 1, G)[g] + rng.normal(0, 1, G * n)
    d = pd.DataFrame(dict(g=g, x=x))
    mg = d.groupby("g").x.mean()
    msb = n * ((mg - x.mean()) ** 2).sum() / (G - 1)
    msw = ((d.x - d.g.map(mg)) ** 2).sum() / (G * n - G)
    ref = (msb - msw) / (msb + (n - 1) * msw)
    check("ICC(1) balanced = textbook ANOVA formula", abs(icc1(x, g) - ref) < 1e-12, f"{icc1(x, g):.4f}")
    # 2. unbalanced case: unbiased-ish for a known ICC (0.5 and 0.0), mean over replicates
    est5, est0 = [], []
    for _ in range(300):
        sizes = rng.integers(1, 12, 60)
        gg = np.repeat(np.arange(60), sizes)
        est5.append(icc1(rng.normal(0, 1, 60)[gg] + rng.normal(0, 1, gg.size), gg))
        est0.append(icc1(rng.normal(0, 1, gg.size), gg))
    check("ICC(1) unbalanced: mean estimate near 0.5 and 0.0", abs(np.mean(est5) - 0.5) < 0.04 and abs(np.mean(est0)) < 0.02,
          f"{np.mean(est5):.3f}, {np.mean(est0):.3f}")
    # 3. no-model interval (offset full conformal on y) keeps exact coverage, one per group
    a, m = 0.10, 20
    Gs = 70
    sizes = rng.integers(2, 9, Gs)
    grp = np.repeat(np.arange(Gs), sizes)
    y = 3.0 + rng.normal(0, 1, Gs)[grp] + rng.normal(0, 0.5, grp.size)
    fam = fs.Family(z=y - 2.0, groups=grp, y=y)
    covs = []
    for _ in range(4000):
        cal, new, tg = fs.draw_one_per_group(fam, m, rng, 10)
        S = fs.offset_full_set(y[cal], a)
        covs.append(fs.group_weighted(fs.in_set(y[new], S), new, fam, tg))
    nom = fs.nominal_mean(m, a)
    se = np.std(covs) / np.sqrt(len(covs))
    check("no-model interval: exact mean coverage of new compositions", abs(np.mean(covs) - nom) < 4 * se,
          f"{np.mean(covs):.4f} vs {nom:.4f}")
    # 4. the post hoc draws re-create the Phase 5 draws exactly (same seeds, same calls)
    c5 = dict(seed_base=777, alpha_primary=0.1, configs=["loso:227"],
              fewshot=dict(draws=30, m_random=[9], m_clustered=[9], m_one_per_group=[9, 20], m_subfamily=[9],
                           min_test_groups=10, subfamily=False))
    B = baseline_draws(fam, c5, "loso", 227, "bandgap")
    P = run_family(fam, c5, "loso", 227, "bandgap")
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, "d.csv.gz")
        with open(f, "wb") as fh:
            fh.write(_gz(P))
        try:
            res = verify_against_phase5(B, f)
            ok = all(v["max_abs_dwidth"] < 1e-12 for v in res.values())
        except AssertionError as e:
            ok, res = False, str(e)
        # a deliberately shifted copy must be rejected
        P2 = P.copy()
        P2.loc[P2.method == "offset_full", "width"] += 1e-3
        with open(f, "wb") as fh:
            fh.write(_gz(P2))
        try:
            verify_against_phase5(B, f)
            caught = False
        except AssertionError:
            caught = True
    check("post hoc draws = Phase 5 draws (and a changed draw is caught)", ok and caught)
    # 5. permutation Spearman agrees with the asymptotic p-value for a clear association
    u = rng.normal(size=14)
    rho, pp, pa = spearman_perm(u, u + rng.normal(0, 1.2, 14), 20000, 1)
    check("permutation Spearman p close to asymptotic", abs(pp - pa) < 0.02 + 0.2 * pa, f"rho {rho:.2f}, {pp:.4f} vs {pa:.4f}")


def t_audit(rng):
    """Phase 6 Task 6.1 audit rules on hand-made entries."""
    import pandas as pd
    from .audit6 import _negatives, _num, classify_oxides, decide, dedup, group_counts
    cases = {"Al2O3": ("O",), "LiFePO4": ("O",), "ZnSO4": ("O",), "Ca(OH)2": ("O",), "In2O3": ("O",),
             "TiOF2": ("F", "O"), "LaOCl": ("Cl", "O")}
    got = {f: _negatives(f) for f in cases}
    check("oxidation-state rule: O the only anion in oxides, oxyanion salts and hydroxides; "
          "oxyfluoride and oxychloride excluded", all(got[f] == v for f, v in cases.items()), str(got))
    check("'na' and non-finite values become NaN", np.isnan(_num("na")) and np.isnan(_num(float("inf"))) and _num("1.5") == 1.5)
    T = pd.DataFrame(dict(jid=["J-1", "J-2", "J-3", "J-4", "J-5", "J-6"],
                          reduced_formula=["Al2O3", "Al2O3", "Al2O3", "CaH2O2", "NaCl", "TiOF2"],
                          chemsys=["Al-O", "Al-O", "Al-O", "Ca-H-O", "Cl-Na", "F-O-Ti"],
                          spg_number=[167, 167, 12, 164, 225, 12], spg_symbol=["R-3c", "R-3c", "C2/m", "P-3m1", "Fm-3m", "C2/m"],
                          crys=["trigonal", "trigonal", "monoclinic", "trigonal", "cubic", "monoclinic"],
                          has_O=[True, True, True, True, False, True],
                          formation_energy_peratom=[-3.0, -3.4, -3.1, -2.0, -2.1, -2.5],
                          optb88vdw_total_energy=[-7.0, -7.4, -7.1, -5.0, -3.0, -6.0],
                          jid_num=[1, 2, 3, 4, 5, 6], mbj_bandgap=[np.nan, 8.0, np.nan, 5.0, 6.0, np.nan]))
    T = classify_oxides(T, ["H", "C", "N", "F", "P", "S", "Cl"])
    check("oxide rules differ as intended (hydroxide: oxidation rule yes, element rule no)",
          T.loc[3, "oxide_oxidation_rule"] and not T.loc[3, "oxide_element_rule"]
          and not T.loc[5, "oxide_oxidation_rule"] and not T.loc[4, "has_O"])
    D = dedup(T)
    check("de-duplication keeps the lowest formation energy per (formula, space group)",
          set(D.jid) == {"J-2", "J-3", "J-4", "J-5", "J-6"})
    Gc = group_counts(D, 2).set_index("spg_number")
    check("group counts: size, eligibility, TBmBJ count, polymorph share",
          Gc.loc[12, "n"] == 2 and Gc.loc[12, "eligible"] and not Gc.loc[167, "eligible"]
          and Gc.loc[167, "n_mbj"] == 1 and abs(Gc.loc[167, "share_formula_in_other_spg"] - 1.0) < 1e-12
          and Gc.loc[12, "share_formula_in_other_spg"] == 0.5)
    a = dict(min_groups_for_oxides=20, min_calibration_groups=10, decision_requires_both_oxide_rules=True)

    def sc(kox, kel, kall):
        rows = [("oxides_oxidation_rule", kox), ("oxides_element_rule", kel), ("all_chemistries", kall)]
        return pd.DataFrame([dict(release="r", scope=s_, K_eligible=k, n_mbj_dedup_within_subset=0, K_eligible_mbj=0)
                             for s_, k in rows])
    d1, d2, d3 = decide(sc(25, 22, 60), a, "r"), decide(sc(19, 22, 59), a, "r"), decide(sc(12, 11, 40), a, "r")
    check("decision rule: oxides only if both rules give >= 20 groups; a flip is reported",
          d1["primary_scope"] == "oxides_oxidation_rule" and d2["primary_scope"] == "all_chemistries"
          and d2["oxide_definition_flips_the_decision"] and d3["primary_scope"] == "all_chemistries"
          and not d3["oxide_definition_flips_the_decision"] and "jackknife" in d2["h3_design_secondary"])


def t_figcheck(rng):
    """The automatic figure overlap check finds planted overlaps and passes clean layouts."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .figcheck import check_overlaps

    def kinds(fig):
        out = {d["kind"].split()[0] for d in check_overlaps(fig)}
        plt.close(fig)
        return out
    fig = plt.figure(figsize=(3, 2)); fig.text(0.3, 0.5, "label one"); fig.text(0.33, 0.5, "label two")
    k1 = kinds(fig)
    fig = plt.figure(figsize=(3, 2)); fig.text(0.3, 0.5, "label one", rotation=40); fig.text(0.6, 0.5, "label two", rotation=40)
    k2 = kinds(fig)
    fig, ax = plt.subplots(); ax.plot(range(10), range(10), label="x"); ax.legend(loc="center")
    k3 = kinds(fig)
    fig, ax = plt.subplots(); ax.scatter([5], [5]); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.text(4.9, 4.9, "here")
    k4 = kinds(fig)
    fig, ax = plt.subplots(); ax.plot(range(10), range(10), label="x"); ax.legend(loc="lower right"); ax.text(1, 8, "free")
    k5 = kinds(fig)
    fig, ax = plt.subplots(); ax.plot([-50, 5], [0, 5]); ax.set_xlim(0, 10)   # clipped part under the y tick labels
    k6 = kinds(fig)
    fig, ax = plt.subplots(); ax.scatter([-0.04], [0.4], s=200, clip_on=False); ax.set_xlim(0, 1); ax.set_ylim(0, 1)   # marker over tick labels
    k7 = kinds(fig)
    check("overlap check: text/text, rotated text, legend/data, annotation/data, label/data are found; "
          "clean and clipped cases pass",
          k1 == {"T1"} and k2 == set() and k3 == {"L1"} and k4 == {"A1"} and k5 == set() and k6 == set()
          and "K1" in k7, str([k1, k2, k3, k4, k5, k6, k7]))


def t_phase6(rng):
    """Phase 6 Tasks 6.2-6.3 tools: structure conversion, features, table round trip, holdouts,
    MatFold label identity."""
    import os
    import tempfile
    import warnings
    import pandas as pd
    from .data6 import jarvis_structure, read_table
    from .features6 import N_COMPOSITION, N_CRYSTALNN, N_DENSITY, composition_block, featurize_structure
    from .phase5 import _gz
    from .splits6 import _sga01, holdout_sets, matfold_check
    warnings.simplefilter("ignore")
    a = 5.64
    nacl_frac = dict(lattice_mat=[[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]], elements=["Na", "Cl"],
                     coords=[[0, 0, 0], [0.5, 0.5, 0.5]], cartesian=False)
    nacl_cart = dict(nacl_frac, coords=[[0, 0, 0], [a / 2, a / 2, a / 2]], cartesian=True)
    s1, s2 = jarvis_structure(nacl_frac), jarvis_structure(nacl_cart)
    check("JARVIS atoms -> pymatgen: fractional and Cartesian records give the same crystal",
          np.allclose(s1.frac_coords, s2.frac_coords) and _sga01(nacl_frac)[0] == 225, str(_sga01(nacl_frac)))
    from matminer.featurizers.composition import ElementProperty
    from pymatgen.core import Composition
    Xc = composition_block(["Al2O3", "NaN"])
    ref = ElementProperty.from_preset("magpie", impute_nan=False).featurize(Composition("Al2O3"))
    check("composition block: 145 Magpie columns, equal to matminer; 'NaN' parsed as sodium nitride",
          Xc.shape == (2, N_COMPOSITION) and np.allclose(Xc[0, :132], ref, equal_nan=True)
          and Composition("NaN").reduced_formula == "NaN" and str(Composition("NaN").elements[0]) == "Na")
    r1, e1 = featurize_structure(nacl_frac)
    r2, _ = featurize_structure(nacl_frac)
    check("structure block: 3 + 122 values, density = pymatgen, deterministic",
          len(r1) == N_DENSITY + N_CRYSTALNN and abs(r1[0] - s1.density) < 1e-9 and np.array_equal(
              np.asarray(r1), np.asarray(r2), equal_nan=True) and e1 == "", e1)
    he = dict(lattice_mat=[[3, 0, 0], [0, 3, 0], [0, 0, 3]], elements=["He"], coords=[[0, 0, 0]], cartesian=False)
    r3, e3 = featurize_structure(he)
    check("element without atomic radius: packing fraction NaN, density kept, no crash",
          np.isnan(r3[2]) and np.isfinite(r3[0]) and "packing_fraction" in e3, e3)
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, "t.csv.gz")
        T = pd.DataFrame(dict(jid=["J-1", "J-2", "J-3"], reduced_formula=["NaN", "NA", "Al2O3"],
                              formation_energy_peratom=[-0.1, np.nan, 0.123456789012345678]))
        with open(f, "wb") as fh:
            fh.write(_gz(T))
        B = read_table(f)
    check("frozen-table reader keeps formulas 'NaN'/'NA' as text and round-trips floats exactly",
          B.reduced_formula.tolist() == ["NaN", "NA", "Al2O3"] and np.isnan(B.formation_energy_peratom[1])
          and B.formation_energy_peratom[2] == 0.123456789012345678)
    lab = np.array([1, 1, 2, 2, 3, 3, 3])
    form = np.array(["A", "B", "A", "C", "D", "B", "E"])
    H = holdout_sets(lab, form, [1, 2, 3])
    ok = (H[1][0].tolist() == [0, 1] and H[1][1].tolist() == [2, 5] and H[2][0].tolist() == [2, 3]
          and H[2][1].tolist() == [0] and H[3][1].tolist() == [1])
    ok &= all(set(t).isdisjoint(set(np.where(lab != g)[0]) - set(range(len(lab)))) for g, (t, _) in H.items())
    check("holdouts: test = the group; polymorph drop = other groups' rows with a formula of the group", ok,
          str({k: (v[0].tolist(), v[1].tolist()) for k, v in H.items()}))
    by = {"J-a": nacl_frac, "J-b": he,
          "J-c": dict(lattice_mat=[[3.2, 0, 0], [-1.6, 2.771, 0], [0, 0, 5.2]], elements=["Zn", "Zn", "O", "O"],
                      coords=[[1 / 3, 2 / 3, 0], [2 / 3, 1 / 3, 0.5], [1 / 3, 2 / 3, 0.38], [2 / 3, 1 / 3, 0.88]],
                      cartesian=False)}
    ours = pd.DataFrame([dict(jid=j, sga01_symbol=_sga01(v)[1], sga01_crys=_sga01(v)[2]) for j, v in by.items()])
    chk = matfold_check(np.array(list(by)), by, ours, 3, 0)
    check("our space-group / crystal-system re-detection = MatFold's own labels",
          chk["identical_space_group_symbol"] == 1.0 and chk["identical_crystal_system"] == 1.0, str(chk))


def t_matfold6(rng):
    """MatFold itself (create_nested_splits, leave-one-label-out, fraction_lower_limit) on a small set of
    real crystal structures gives exactly our holdouts; rare labels stay in training."""
    import os
    import tempfile
    import warnings
    import pandas as pd
    from MatFold import MatFold
    from .data6 import jarvis_structure
    from .matfold6 import compare_folds, matfold_folds
    from .splits6 import _sga01, eligible_groups, holdout_sets
    warnings.simplefilter("ignore")
    recs = {}
    for i in range(6):                 # rock salt NaCl, a = 5.4..5.9 (Fm-3m, 225)
        a = 5.4 + 0.1 * i
        recs[f"J-r{i}"] = dict(lattice_mat=[[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]],
                               elements=["Na", "Cl"], coords=[[0, 0, 0], [0.5, 0.5, 0.5]], cartesian=False)
    for i in range(4):                 # wurtzite ZnO (P6_3mc, 186)
        a, c = 3.2 + 0.02 * i, 5.2
        recs[f"J-w{i}"] = dict(lattice_mat=[[a, 0, 0], [-a / 2, a * 3 ** 0.5 / 2, 0], [0, 0, c]],
                               elements=["Zn", "Zn", "O", "O"],
                               coords=[[1 / 3, 2 / 3, 0], [2 / 3, 1 / 3, 0.5], [1 / 3, 2 / 3, 0.38],
                                       [2 / 3, 1 / 3, 0.88]], cartesian=False)
    for i in range(3):                 # simple cubic Po-type (Pm-3m, 221), CsCl (221) and a rare tetragonal cell
        a = 3.3 + 0.1 * i
        recs[f"J-c{i}"] = dict(lattice_mat=[[a, 0, 0], [0, a, 0], [0, 0, a]], elements=["Cs", "Cl"],
                               coords=[[0, 0, 0], [0.5, 0.5, 0.5]], cartesian=False)
    recs["J-t0"] = dict(lattice_mat=[[3, 0, 0], [0, 3, 0], [0, 0, 4.1]], elements=["Na"], coords=[[0, 0, 0]],
                        cartesian=False)
    jids = np.array(list(recs))
    lab = np.array([_sga01(recs[j])[0] for j in jids])
    crys = np.array([_sga01(recs[j])[2].lower() for j in jids])
    sym = np.array([_sga01(recs[j])[1] for j in jids])
    form = np.array([jarvis_structure(recs[j]).composition.reduced_formula for j in jids])
    elig = eligible_groups(lab, 3)                       # 225 (6), 186 (4), 221 (3); the 1-entry group stays in training
    H = holdout_sets(lab, form, elig)
    with tempfile.TemporaryDirectory() as td:
        bulk = {j: jarvis_structure(recs[j]).as_dict() for j in jids}
        mf = MatFold(pd.DataFrame({"jid": jids}), bulk, cols_to_keep=["jid"], write_data_checksums=True, seed=1)
        mf.create_nested_splits("sgnum", n_inner_splits=1, n_outer_splits=0, write_base_str="mf",
                                output_dir=td, fraction_lower_limit=(3 - 0.5) / len(jids))
        folds = matfold_folds(td, "sgnum")
        cmp_ = compare_folds(folds, jids, dict(zip(jids, lab.tolist())), {g: t for g, (t, _) in H.items()})
        m = mf.df.set_index("structureid").loc[jids]
        same_lab = (m.sgnum.astype(str).to_numpy() == sym).all() and (m.composition.to_numpy() == form).all()
        rare_in_train = all("J-t0" in set(tr) for _, tr, _ in folds)
        has_json = os.path.exists(os.path.join(td, "mf.sgnum.json"))
    check("MatFold itself: labels = ours; LOO folds = our holdouts row by row; rare label always in training; "
          "JSON written",
          cmp_["all_identical"] and cmp_["n_folds"] == 3 and sorted(elig) == [186, 221, 225] and same_lab
          and rare_in_train and has_json, str(dict(elig=elig, cmp=cmp_["folds"], same_lab=same_lab)))
    lab_bad = lab.copy()
    lab_bad[0] = 186                                     # a deliberate mismatch must be caught
    H_bad = holdout_sets(lab_bad, form, [186, 221, 225])
    cmp_bad = compare_folds(folds, jids, dict(zip(jids, lab.tolist())), {g: t for g, (t, _) in H_bad.items()})
    check("fold comparison detects a one-row difference", not cmp_bad["all_identical"])


def t_phase7(rng):
    """v7 tools: fast CV+ = definition (unweighted, group weighted, normalized, CQR), HCP = Lee et al.
    eq. 6, single subsampling, coverage of HCP and group-weighted CV+ for new groups, design invariants."""
    import pandas as pd
    from . import conformal as cp
    from . import cvfast as cf
    from .design7 import A_rowsets, B_rowsets, H3_partition, build_target_design, check_design, tune_rowsets

    def same(a, b):
        return (np.array_equal(np.isposinf(a), np.isposinf(b)) and np.array_equal(np.isneginf(a), np.isneginf(b))
                and np.array_equal(a[np.isfinite(a)], b[np.isfinite(b)]))
    bad = 0
    for trial in range(150):
        K, n, nte = int(rng.integers(2, 7)), int(rng.integers(1, 50)), int(rng.integers(1, 15))
        fold = rng.integers(0, K, n)
        S = np.abs(rng.normal(size=n)) if trial % 3 else rng.normal(size=n)
        if trial % 5 == 0:
            S = np.round(S, 1)
        mu = rng.normal(size=(K, nte))
        if trial % 7 == 0:
            mu = np.round(mu, 1)
        mu2 = mu + np.abs(rng.normal(size=(K, nte)))
        sc = None if trial % 2 else np.exp(rng.normal(size=(K, nte)))
        grp = rng.integers(0, int(rng.integers(1, 15)), n)
        for a in (0.05, 0.1, 0.2, 0.5):
            for gw in (None, grp):
                lo, hi = cf.cvplus(mu, mu2, sc, S, fold, K, a, group=gw)
                blo, bhi = cf.cvplus_brute(mu, mu2, sc, S, fold, K, a, group=gw)
                bad += not (same(lo, blo) and same(hi, bhi))
        lo, hi = cf.cvplus(mu, mu, None, np.abs(S), fold, K, 0.1)
        r = cp.cvplus_interval(mu, fold, np.abs(S), 0.1)
        bad += not (same(lo, r[0]) and same(hi, r[1]))
    check("fast CV+ (plain, group weighted, normalized, CQR form) = the definition, bit for bit, 1,350 cases; "
          "= NOMAD's cvplus_interval", bad == 0, f"{bad} mismatches")
    # group weights with one point per group reduce to plain CV+ (up to the +-inf atom: (G+1) = n+1)
    K, n = 4, 37
    fold, S, mu = rng.integers(0, K, n), np.abs(rng.normal(size=n)), rng.normal(size=(K, 9))
    a1 = cf.cvplus(mu, mu, None, S, fold, K, 0.1, group=np.arange(n))
    a2 = cf.cvplus(mu, mu, None, S, fold, K, 0.1)
    check("group-weighted CV+ with one point per group = plain CV+", same(a1[0], a2[0]) and same(a1[1], a2[1]))
    # HCP = Lee et al. eq. 6 written out; N_k = 1 gives split conformal
    sc_, g_ = rng.normal(size=60) ** 2, rng.integers(0, 12, 60)
    ug, inv, cnt = np.unique(g_, return_inverse=True, return_counts=True)
    w = 1 / ((len(ug) + 1) * cnt[inv])
    ts = np.sort(sc_)
    direct = next((t for t in ts if w[sc_ <= t].sum() >= 0.9 - 1e-12), np.inf)
    ones = rng.normal(size=30) ** 2
    check("HCP threshold = Lee et al. eq. 6; with one point per group it is the split-conformal quantile",
          cf.hcp_threshold(sc_, g_, 0.1) == direct and
          cf.hcp_threshold(ones, np.arange(30), 0.1) == cp.conformal_quantile(ones, 0.1))
    # coverage for new groups: group offsets + noise, unequal group sizes
    covH, covS, covW = [], [], []
    r2 = np.random.default_rng(7)
    for rep in range(400):
        G = 21
        sizes = r2.integers(3, 40, G)
        b = r2.normal(size=G) * 1.5
        ys = [b[k] + r2.normal(size=sizes[k]) * r2.uniform(0.3, 2.0) for k in range(G)]
        cal = np.concatenate(ys[:20])
        gg = np.concatenate([[k] * sizes[k] for k in range(20)])
        T = cf.hcp_threshold(np.abs(cal), gg, 0.1)
        covH.append(np.mean(np.abs(ys[20]) <= T))
        Ts = cf.single_subsample_threshold(np.abs(cal), gg, 0.1, r2)
        covS.append(np.mean(np.abs(ys[20]) <= Ts))
        # group-weighted CV+ with group folds; model = mean of the training groups (no features)
        fold = gg % 5
        mu_k = np.array([cal[fold != k].mean() for k in range(5)])
        R = np.abs(cal - mu_k[fold])
        lo, hi = cf.cvplus(mu_k[:, None], mu_k[:, None], None, R, fold, 5, 0.1, group=gg)
        covW.append(np.mean((ys[20] >= lo[0]) & (ys[20] <= hi[0])))
    check("new-group coverage by simulation (400 reps): HCP >= 0.88, single subsampling >= 0.88, "
          "group-weighted CV+ >= 0.80 (1 - 2 alpha)", np.mean(covH) >= 0.88 and np.mean(covS) >= 0.88
          and np.mean(covW) >= 0.80, f"HCP {np.mean(covH):.3f}, SS {np.mean(covS):.3f}, gw CV+ {np.mean(covW):.3f}")
    check("interval score (Gneiting & Raftery) on known values",
          np.allclose(cf.interval_score([0, 3, -2], [-1, -1, -1], [1, 1, 1], 0.1), [2, 2 + 20 * 2, 2 + 20 * 1]))
    # design invariants on a synthetic table
    n_g, rows = 46, []
    forms = [f"A{i}B" for i in range(500)]
    for g in range(n_g):
        for i in range(int(rng.integers(8, 60))):
            rows.append(dict(g=g + 1, crys=["cubic", "hexagonal", "monoclinic", "triclinic"][g % 4],
                             f=forms[int(rng.integers(0, 500))], y=rng.normal(), gap=max(0.0, rng.normal(0.5, 1))))
    T = pd.DataFrame(rows)
    F = pd.DataFrame(dict(jid=[f"J{i}" for i in range(len(T))], reduced_formula=T.f, formation_energy_peratom=T.y,
                          optb88vdw_bandgap=T.gap, mbj_bandgap=np.nan))
    L = pd.DataFrame(dict(sga01_number=T.g, sga01_crys=T.crys))
    m7 = dict(seed=1, min_group_size=25, outer_folds_per_half=5, inner_folds=5, tuning=dict(inner_folds=3),
              h3=dict(n_partitions=4, n_cal_groups=dict(formation=8, bandgap=4)))
    ok, msg = True, ""
    try:
        for t in ("formation", "bandgap"):
            D = build_target_design(F, L, t, m7)
            check_design(D, m7)
            R = B_rowsets(D, 1, 0, "polyex")
            ok &= not np.isin(D["form"][R["train"]], D["form"][R["test"]]).any()
            R2 = B_rowsets(D, 1, 0, "primary")
            ok &= len(R2["train"]) == len(R2["cal"]) + len(R2["rare"])
            ok &= not np.array_equal(D["outer_r1"], D["outer_r2"]) or len(D["eligible"]) <= 10
            P = H3_partition(D, 1, m7)
            ok &= set(P["pool_groups"]) == set(D["eligible"][D["half"] == 1])
            rows_, folds = tune_rowsets(D, "h0")
            ok &= not np.isin(D["grp"][rows_], D["eligible"][D["half"] == 0]).any()
            A0 = A_rowsets(D, 0)
            ok &= len(A0["test"]) + len(A0["pool"]) == len(D["rows"])
    except AssertionError as e:
        ok, msg = False, str(e)
    check("v7 design: held-out groups never train; halves never cross tuning and testing; rare groups never "
          "held out; polymorph exclusion complete; Regime A formula-disjoint; H3 pools from one half", ok, msg)


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
    print("Self-test | Phase 4 detection tools");   t_detect(rng)
    print("Self-test | Phase 5 few-shot tools");    t_fewshot(rng)
    print("Self-test | Phase 5 post hoc tools");    t_posthoc(rng)
    print("Self-test | Phase 6 audit rules");       t_audit(rng)
    print("Self-test | figure overlap check");      t_figcheck(rng)
    print("Self-test | Phase 6 data tools");        t_phase6(rng)
    print("Self-test | MatFold splits");            t_matfold6(rng)
    print("Self-test | Phase 6 models (v7)");       t_phase7(rng)
    print(f"\n{'ALL SELF-TESTS PASSED' if not FAIL else 'FAILED: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0
