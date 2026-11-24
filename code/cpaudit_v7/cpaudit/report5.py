"""Phase 5 report: summaries, pre-registered decisions (D5.0, H5, H5b), practical endpoint,
exploratory tables, and the optional training-arm summary."""
from __future__ import annotations

import glob
import json
import math
import os

import numpy as np
import pandas as pd

from .config import atomic_write_csv, atomic_write_json, atomic_write_text, get_paths
from .data import LATTICE_NAME
from .fewshot import beta_reference, exact_summary, nominal_mean
from .setup5 import prereg_paths

KEYS = ["design", "heldout", "sg", "target", "scheme", "null", "m", "method"]


def load_draws(paths) -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(paths.out, "draws", "*.csv.gz")))
    if not files:
        raise RuntimeError("No Phase 5 draws found. Run the `phase5` step first.")
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)


def summarise(D: pd.DataFrame, fams: pd.DataFrame, cfg) -> pd.DataFrame:
    a = float(cfg["alpha_primary"])
    fl = float(cfg["decisions"]["coverage_floor"])
    g = D.groupby(KEYS, sort=False)
    S = g.agg(n_draws=("cov", "size"), mean_cov=("cov", "mean"), sd_cov=("cov", "std"),
              q05=("cov", lambda x: x.quantile(0.05)), q25=("cov", lambda x: x.quantile(0.25)),
              q50=("cov", "median"), q75=("cov", lambda x: x.quantile(0.75)),
              q95=("cov", lambda x: x.quantile(0.95)), p_below=("cov", lambda x: float(np.mean(x < fl))),
              mean_cov_new=("cov_new", "mean"), mean_cov_family=("cov_family", "mean"),
              width_median=("width", "median"),
              frac_infinite=("infinite", "mean"), mean_components=("n_comp", "mean"),
              mae_raw=("mae_raw", "mean"), mae_offset=("mae_offset", "mean"), b_hat=("b_hat", "mean"),
              n_test=("n_test", "mean")).reset_index()
    S = S.merge(fams[["design", "sg", "target", "N", "G", "sd_y", "bias_full", "mae_full"]],
                on=["design", "sg", "target"], how="left")
    S["width_sd"] = S.width_median / S.sd_y
    S["nominal"] = [nominal_mean(int(m), a) for m in S.m]
    ref = [beta_reference(int(m), a, fl) for m in S.m]
    for k in ("k", "beta_sd", "beta_p_below"):
        S[k] = [r[k] for r in ref]
    ex = [exact_summary(int(N), int(m), a, fl) if (sc == "random" and meth == "split") else
          dict(exact_mean=np.nan, exact_sd=np.nan, exact_p_below=np.nan)
          for N, m, sc, meth in zip(S.N, S.m, S.scheme, S.method)]
    for k in ("exact_mean", "exact_sd", "exact_p_below"):
        S[k] = [e[k] for e in ex]
    return S


def loso_reference(paths) -> pd.DataFrame:
    rows = []
    for fn, design in (("phase2_loso_all.csv", "loso"), ("phase2_l2o_all.csv", "l2o")):
        p = os.path.join(paths.phase2, "report", fn)
        if not os.path.exists(p):
            continue
        d = pd.read_csv(p)
        d = d[(d.method == "cvplus_lattice") & (d.group == "heldout") & np.isclose(d.alpha, 0.1)]
        for _, r in d.iterrows():
            rows.append(dict(design=design, heldout=r.heldout, target=r.target, loso_coverage=r.coverage,
                             loso_width=r.width_median))
    return pd.DataFrame(rows)


def d50(S, cfg, paths=None):
    """Pre-registered code check. A cell outside the band is re-run once with 10x the draws and a
    new seed (pre-registered), so a single Monte Carlo fluke among ~170 cells is not read as a bug."""
    n = int(cfg["fewshot"]["draws"])
    k = float(cfg["decisions"]["sanity_n_se"])
    R = S[(S.scheme == "random") & (~S.null) & S.method.isin(["split", "offset_full"])].copy()
    R["se"] = np.where(R.method == "split", R.exact_sd, R.sd_cov) / math.sqrt(n)
    R["z"] = (R.mean_cov - R.nominal) / R.se
    R["ok"] = R.z.abs() <= k
    R["rechecked"] = False
    if paths is not None and (~R.ok).any():
        from .data import id_table
        from .fewshot import SET_FUNCS, in_set
        from .phase5 import _seed, load_family
        tab = id_table(cfg, paths)
        a = float(cfg["alpha_primary"])
        for i in R.index[~R.ok]:
            row = R.loc[i]
            fam = load_family(paths, tab, row.design, int(row.sg), row.target)
            rng = np.random.default_rng(_seed(cfg["seed_base"], "recheck", row.design, row.sg, row.target, row.m))
            cov = []
            for _ in range(10 * n):
                p = rng.permutation(fam.N)
                cov.append(in_set(fam.z[p[int(row.m):]], SET_FUNCS[row.method](fam.z[p[:int(row.m)]], a)).mean())
            cov = np.array(cov)
            se = cov.std(ddof=1) / math.sqrt(cov.size)
            R.loc[i, ["mean_cov", "se", "z"]] = [cov.mean(), se, (cov.mean() - row.nominal) / se]
            R.loc[i, "ok"] = abs((cov.mean() - row.nominal) / se) <= k
            R.loc[i, "rechecked"] = True
    return R[["design", "heldout", "target", "m", "method", "mean_cov", "nominal", "se", "z", "ok", "rechecked",
              "p_below", "exact_p_below"]], bool(R.ok.all())


def _boot_ratio(x, y, B, rng):
    x, y = np.asarray(x), np.asarray(y)
    r = [np.std(rng.choice(x, x.size), ddof=1) / np.std(rng.choice(y, y.size), ddof=1) for _ in range(B)]
    return np.percentile(r, [2.5, 97.5])


def dispersion(D, cfg, seed):
    """SD of coverage across draws relative to the permutation null (same scheme, same group sizes,
    residuals permuted across the family), split method, every m; bootstrap CI at the H5 m."""
    r = cfg["decisions"]
    m5 = int(r["h5_m"])
    rng = np.random.default_rng(seed)
    rows = []
    sel = D[D.scheme.isin(["clustered", "one_per_group"]) & (D.method == "split")]
    for (scheme, design, sg, target, mm), G in sel.groupby(["scheme", "design", "sg", "target", "m"]):
        x, y = G[~G.null]["cov"].to_numpy(), G[G.null]["cov"].to_numpy()
        if len(x) < 10 or len(y) < 10:
            continue
        lo, hi = _boot_ratio(x, y, 1000, rng) if mm == m5 else (np.nan, np.nan)
        rows.append(dict(scheme=scheme, design=design, heldout=LATTICE_NAME[sg], sg=sg, target=target, m=mm,
                         mean_real=x.mean(), mean_null=y.mean(), nominal=nominal_mean(int(mm), cfg["alpha_primary"]),
                         sd_real=np.std(x, ddof=1), sd_null=np.std(y, ddof=1),
                         rho=np.std(x, ddof=1) / np.std(y, ddof=1), rho_lo=lo, rho_hi=hi,
                         p_below_real=np.mean(x < r["coverage_floor"]), p_below_null=np.mean(y < r["coverage_floor"])))
    return pd.DataFrame(rows)


def h5(T, cfg):
    """H5 (one_per_group) and H5c (clustered) decisions at the pre-registered m."""
    r = cfg["decisions"]
    m5 = int(r["h5_m"])
    dec = {}
    for target, Tt in T[(T.m == m5) & (T.scheme == "one_per_group")].groupby("target"):
        sup = int((Tt.rho <= r["h5_ratio_support"]).sum())
        rej = int((Tt.rho > r["h5_ratio_reject"]).sum())
        verdict = ("supported" if sup >= r["h5_support_min"] else "rejected" if rej >= r["h5_reject_min"]
                   else "partial")
        dec[f"H5_{target}"] = dict(verdict=verdict, n_families=len(Tt), n_rho_le_support=sup,
                                   n_rho_gt_reject=rej,
                                   rho={f"{h} ({d})": round(float(v), 3) for h, d, v in zip(Tt.heldout, Tt.design, Tt.rho)})
    for target, Tt in T[(T.m == m5) & (T.scheme == "clustered")].groupby("target"):
        harm = (Tt.rho > r["h5_ratio_support"]) | (Tt.mean_real < Tt.nominal - r["h5c_mean_drop"])
        n = int(harm.sum())
        dec[f"H5c_{target}"] = dict(verdict=("supported" if n >= r["h5c_min_families"] else "not supported"),
                                    n_families=len(Tt), n_harmed=n,
                                    rho={f"{h} ({d})": round(float(v), 3) for h, d, v in zip(Tt.heldout, Tt.design, Tt.rho)},
                                    mean_cov_new={f"{h} ({d})": round(float(v), 3)
                                                  for h, d, v in zip(Tt.heldout, Tt.design, Tt.mean_real)})
    return dec


def h5b(S, cfg):
    r = cfg["decisions"]
    m = int(r["h5b_m"])
    rows, dec = [], {}
    base = S[(S.scheme == "one_per_group") & (~S.null) & (S.m == m)]
    for target in ("bandgap", "formation"):
        cells = [c.split(":") for c in r["h5b_cells"][target]]
        for design, g in cells:
            g = int(g)
            q = base[(base.design == design) & (base.sg == g) & (base.target == target)].set_index("method")
            if not {"split", "offset_full"} <= set(q.index):
                continue
            ratio = q.loc["offset_full", "width_median"] / q.loc["split", "width_median"]
            cov = q.loc["offset_full", "mean_cov"]
            rows.append(dict(target=target, design=design, heldout=LATTICE_NAME[g], m=m,
                             width_split=q.loc["split", "width_median"], width_offset=q.loc["offset_full", "width_median"],
                             width_ratio=ratio, cov_split=q.loc["split", "mean_cov"], cov_offset=cov,
                             passes=bool(ratio <= r["h5b_width_ratio"] and cov >= r["h5b_min_coverage"])))
        T = pd.DataFrame([x for x in rows if x["target"] == target])
        n = len(T)
        need = math.ceil(r["h5b_share"] * n - 1e-9)
        npass = int(T.passes.sum()) if n else 0
        nrej = int((T.width_ratio > r["h5b_width_reject"]).sum()) if n else 0
        verdict = ("supported" if n and npass >= need else "rejected" if n and nrej >= math.ceil(n / 2)
                   else "partial")
        dec[target] = dict(verdict=verdict, n_cells=n, n_pass=npass, needed=need, n_ratio_gt_reject=nrej)
    return pd.DataFrame(rows), dec


def m_star(S, cfg):
    r = cfg["decisions"]
    out = []
    sub = S[(~S.null) & S.scheme.isin(["random", "clustered", "one_per_group"])]
    for (design, sg, target, scheme, method), G in sub.groupby(["design", "sg", "target", "scheme", "method"]):
        G = G.sort_values("m")
        ok = G[(G.mean_cov >= r["practical_mean"]) & (G.p_below <= r["practical_risk"]) & (G.frac_infinite == 0)]
        out.append(dict(design=design, heldout=LATTICE_NAME[sg], target=target, scheme=scheme, method=method,
                        m_star=(int(ok.m.iloc[0]) if len(ok) else np.nan), m_max=int(G.m.max()),
                        reached=bool(len(ok))))
    return pd.DataFrame(out)


def clustering_table(S):
    a = S[(~S.null) & (S.method == "split") & S.scheme.isin(["clustered", "one_per_group"])]
    p = a.pivot_table(index=["design", "heldout", "target", "m"], columns="scheme", values="sd_cov").reset_index()
    p = p.dropna()
    p["sd_ratio_clustered_to_one_per_group"] = p["clustered"] / p["one_per_group"]
    return p


def train_summary(paths):
    files = sorted(glob.glob(os.path.join(paths.out, "train", "draw_*.csv")))
    if not files:
        return None
    T = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    return T.groupby(["design", "heldout", "target", "m", "arm", "method"]).agg(
        n_draws=("cov", "size"), mean_cov=("cov", "mean"), sd_cov=("cov", "std"),
        mean_cov_group=("cov_group", "mean"), width_median=("width", "median"), mae=("mae", "mean")).reset_index()


def _flat(df):
    df = df.copy()
    df.columns = [c if isinstance(c, str) else " ".join(x for x in c if x) for c in df.columns]
    return df


def _md(df, cols=None, digits=3):
    d = df if cols is None else df[cols]
    return d.round(digits).to_markdown(index=False)


def run_report5(cfg):
    paths = get_paths(cfg)
    js, _ = prereg_paths(paths)
    if not os.path.exists(js):
        raise RuntimeError("No Phase 5 pre-registration found.")
    rep = paths.report
    D = load_draws(paths)
    fams = pd.read_csv(os.path.join(paths.out, "phase5_families.csv"))
    S = summarise(D, fams, cfg)
    L = loso_reference(paths)
    if len(L):
        S = S.merge(L, on=["design", "heldout", "target"], how="left")
        S["loso_width_sd"] = S.loso_width / S.sd_y
    atomic_write_csv(os.path.join(rep, "phase5_summary.csv"), S)

    seed = int(cfg["seed_base"]) + 5000
    san, san_ok = d50(S, cfg, paths)
    H = dispersion(D, cfg, seed)
    hdec = h5(H, cfg)
    B, bdec = h5b(S, cfg)
    M = m_star(S, cfg)
    C = clustering_table(S)  # SD of clustered vs one_per_group at equal m
    atomic_write_csv(os.path.join(rep, "phase5_sanity.csv"), san)
    atomic_write_csv(os.path.join(rep, "phase5_dispersion.csv"), H)
    atomic_write_csv(os.path.join(rep, "phase5_offset.csv"), B)
    atomic_write_csv(os.path.join(rep, "phase5_mstar.csv"), M)
    atomic_write_csv(os.path.join(rep, "phase5_clustering.csv"), C)
    sub = S[S.scheme.isin(["sub_a2b", "sub_b2a"])][["target", "scheme", "m", "method", "mean_cov", "q05",
                                                    "p_below", "width_median", "b_hat"]]
    atomic_write_csv(os.path.join(rep, "phase5_subfamily.csv"), sub)
    TR = train_summary(paths)
    if TR is not None:
        atomic_write_csv(os.path.join(rep, "phase5_train_summary.csv"), TR)

    decisions = dict(D5_0_code_check_passed=san_ok, n_sanity_cells=int(len(san)),
                     n_sanity_failures=int((~san.ok).sum()), H5=hdec, H5b=bdec)
    atomic_write_json(os.path.join(rep, "phase5_decisions.json"), decisions)

    m5, mb = int(cfg["decisions"]["h5_m"]), int(cfg["decisions"]["h5b_m"])
    W = S[(~S.null) & (S.scheme == "one_per_group") & (S.m == mb)].pivot_table(
        index=["design", "heldout", "target"], columns="method", values=["mean_cov", "width_sd"]).reset_index()
    W.columns = ["_".join([c for c in col if c]) for col in W.columns]
    if len(L):
        W = W.merge(S.drop_duplicates(["design", "heldout", "target"])[["design", "heldout", "target",
                                                                         "loso_coverage", "loso_width_sd"]],
                    on=["design", "heldout", "target"], how="left")
    atomic_write_csv(os.path.join(rep, "phase5_widths.csv"), W)
    mae = S[(~S.null) & (S.scheme == "one_per_group") & (S.m == mb) & (S.method == "split")][
        ["design", "heldout", "target", "mae_full", "mae_raw", "mae_offset", "bias_full", "b_hat"]]

    lines = [f"# Phase 5 summary (alpha = {cfg['alpha_primary']}; {cfg['fewshot']['draws']} draws per cell)", "",
             f"## D5.0 code check (random scheme): {'PASSED' if san_ok else 'FAILED, do not interpret'} "
             f"({int((~san.ok).sum())} of {len(san)} cells outside {cfg['decisions']['sanity_n_se']} SE)", "",
             f"## H5 / H5c dispersion vs permutation null (split, m = {m5}; coverage of new compositions)", "",
             _md(H[H.m == m5], ["scheme", "heldout", "design", "target", "mean_real", "mean_null", "nominal",
                                "sd_real", "sd_null", "rho", "rho_lo", "rho_hi", "p_below_real", "p_below_null"]), "",
             "```", json.dumps(hdec, indent=2), "```", "",
             f"## H5b offset correction (one_per_group, m = {mb})", "", _md(B) if len(B) else "(no cells)", "",
             "```", json.dumps(bdec, indent=2), "```", "",
             "## Practical endpoint m* (mean >= {pm}, P(cov < {fl}) <= {pr})".format(
                 pm=cfg["decisions"]["practical_mean"], fl=cfg["decisions"]["coverage_floor"],
                 pr=cfg["decisions"]["practical_risk"]), "",
             _md(_flat(M.assign(m_star=M.m_star.where(M.reached, -1)).pivot_table(
                 index=["heldout", "design", "target"], columns=["scheme", "method"], values="m_star").reset_index())),
             "(-1 = not reached within the m grid)", "",
             f"## Widths at m = {mb} (one_per_group), in units of the family SD, next to LOSO CV+", "", _md(W), "",
             f"## MAE before / after subtracting the estimated offset (m = {mb}, one_per_group)", "", _md(mae), "",
             "## Clustered vs one-per-group dispersion (exploratory)", "", _md(C), "",
             "## Fd-3m sub-family transfer (exploratory)", "", _md(sub), ""]
    if TR is not None:
        lines += ["## Training arm (descriptive)", "", _md(TR), ""]
    txt = "\n".join(lines)
    atomic_write_text(os.path.join(rep, "PHASE5_SUMMARY.md"), txt)
    print(txt)
    print(f"\nReport written to {rep}")
