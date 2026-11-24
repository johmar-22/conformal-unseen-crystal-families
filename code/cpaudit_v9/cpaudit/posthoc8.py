"""Phase 6 Task 6.5 post hoc analyses (v8.1). NOT pre-registered: designed after the v8 results
were seen (2026-11-19). Every output is labelled as post hoc; no pre-registered decision changes.

Four questions raised by the v8 results, answered from the saved Drive files (the v7 units and
the v8 outputs are only read):

1. Is the recalibrated few-shot interval more useful than no model at all?  The label sets of
   fewshot8 are re-created with the registered seeds and the same calls (scheme, m, draw), and
   checked draw by draw against the saved fewshot8 draws (model widths and coverages must be
   identical). On each label set a no-model interval is built from the SAME labels: the
   offset-corrected full-conformal interval with mu = 0, i.e. centred on the mean of the labelled
   y values (the posthoc5 definition). It is valid for the same reason as the model interval.
   The coverage spread over draws (5 % and 95 % quantiles) is also tabulated.
2. Does the inner-model SD (MatFold's ensemble SD) also rank failing GROUPS?  Spearman of
   repeat-1 group coverage with the group median of the SD (and, for comparison, of the CV+ width
   and of the kNN distance), with a bootstrap over groups and a permutation p-value; and the
   coverage retained inside the three worst groups when the highest-SD structures are abstained.
3. Oxides: are the worst-covered oxide groups driven by ABO3 (perovskite-type) compositions?
   ABO3 = reduced composition with three elements, O3 and the two cations 1 : 1 (a composition
   rule only; the crystal structure is not checked).
4. Which group-level intervals are degenerate?  With K calibration groups the point mass at
   +infinity is 1 / (K + 1): the interval is unbounded if 1 / (K + 1) > alpha, and its upper end is
   the largest calibration term if 1 / (K + 1) = alpha (hierarchical jackknife+, group-weighted
   CV+, HCP). Single subsampling needs ceil((K + 1)(1 - alpha)) <= K (equality: the largest score).

Outputs: phase08_v1/posthoc/ (new folder). A rerun recomputes everything and stops if a table
would change; existing files are never replaced by different content.
"""
from __future__ import annotations

import io
import json
import math
import os
import time

import numpy as np
import pandas as pd

from .config import atomic_write_bytes, code_hash, sha256_file

ATOL_TIE = 1e-12


# ----------------------------------------------------------------------------- small rules
def atom_status(K: int, alpha: float) -> str:
    """Weighted quantile with a point mass 1 / (K + 1) at +inf (Lee et al. eq. 6 / eq. 20)."""
    atom = 1.0 / (K + 1)
    if atom > alpha + ATOL_TIE:
        return "unbounded"
    if abs(atom - alpha) <= ATOL_TIE:
        return "largest term"
    return "finite"


def single_subsample_status(K: int, alpha: float) -> str:
    """ceil((K + 1)(1 - alpha))-th smallest of K scores (one structure per calibration group)."""
    k = math.ceil((K + 1) * (1 - alpha) - 1e-9)
    if k > K:
        return "unbounded"
    if k == K:
        return "largest term"
    return "finite"


def is_abo3(formula: str) -> bool:
    """Reduced composition A1 B1 O3 (three elements). Composition rule only."""
    from pymatgen.core import Composition
    try:
        c = Composition(formula).reduced_composition
    except Exception:
        return False
    d = {str(e): float(v) for e, v in c.get_el_amt_dict().items()}
    if len(d) != 3 or abs(d.get("O", 0.0) - 3.0) > 1e-8:
        return False
    return sorted(round(v, 8) for k, v in d.items() if k != "O") == [1.0, 1.0]


# ----------------------------------------------------------------------------- write-once
def _bytes_csv(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False, float_format="%.10g").encode()


def write_once(path: str, data: bytes) -> str:
    """Write a new file; if it exists it must already hold exactly these bytes."""
    if os.path.exists(path):
        with open(path, "rb") as f:
            old = f.read()
        if old != data:
            raise RuntimeError(f"{path} exists with different content; nothing was overwritten. "
                               "Move the old posthoc folder aside if a recomputation is intended.")
        return "unchanged"
    atomic_write_bytes(path, data)
    return "written"


def _gz_bytes(df: pd.DataFrame) -> bytes:
    import gzip
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as g:   # mtime 0: identical bytes on rerun
        g.write(_bytes_csv(df))
    return buf.getvalue()


# ----------------------------------------------------------------------------- 1. few-shot
REAL_SCHEMES = (("random", False), ("clustered", False), ("one_per_group", False))


def nomodel_draws(fam, m_grid, draws, min_test_groups, alpha, seed, label) -> pd.DataFrame:
    """Re-create the fewshot8 label sets (run8.run_family8: same seeds, same calls, same order of
    random draws) for the real (non-permuted) schemes, and evaluate on each: the split and
    offset-corrected model intervals (z = y - mu) and the no-model interval (offset-corrected on y)."""
    from .design7 import _seed
    from .fewshot import draw_clustered, draw_one_per_group, draw_random, group_weighted
    from .phase5 import _evaluate
    allidx = np.arange(fam.N)
    out = []
    for scheme, null in REAL_SCHEMES:
        for m in m_grid:
            rng = np.random.default_rng(_seed(seed, label, scheme, m, null) % (2 ** 32 - 1))
            for d in range(int(draws)):
                new, tg = None, None
                if scheme == "random":
                    cal, rest = draw_random(fam, m, rng)
                else:
                    r = (draw_clustered if scheme == "clustered" else draw_one_per_group)(fam, m, rng, min_test_groups)
                    if r is None:
                        break
                    cal, new, tg = r
                    rest = np.setdiff1d(allidx, cal, assume_unique=True)
                new_fn = ((lambda c, _t=new, _g=tg: group_weighted(c, _t, fam, _g)) if new is not None else None)
                row = dict(scheme=scheme, null=null, m=int(m), draw=d)
                for name, v in (("model", fam.z), ("nomodel", fam.y)):
                    R = _evaluate(v[cal], (v[new] if new is not None else None), v[rest], alpha, new_fn)
                    for e in R:
                        if name == "nomodel" and e["method"] != "offset_full":
                            continue
                        key = "nomodel" if name == "nomodel" else e["method"]
                        row[f"width_{key}"] = e["width"]
                        row[f"cov_{key}"] = e["cov"]
                        row[f"covfam_{key}"] = e["cov_family"]
                out.append(row)
    return pd.DataFrame(out)


def verify_against_fewshot8(B: pd.DataFrame, saved: pd.DataFrame) -> dict:
    """Model widths and coverages of the re-created draws must equal the saved fewshot8 draws."""
    P = saved[saved.null.astype(str) == "False"]
    res = {}
    for meth in ("split", "offset_full"):
        Q = P[P.method == meth].set_index(["scheme", "m", "draw"]).sort_index()
        R = B.set_index(["scheme", "m", "draw"]).sort_index()
        if len(Q) != len(R) or not Q.index.equals(R.index):
            raise AssertionError(f"re-created draw index differs from fewshot8 ({meth}): {len(R)} vs {len(Q)}")
        dw = float(np.nanmax(np.abs(Q["width"].to_numpy() - R[f"width_{meth}"].to_numpy())))
        dc = float(np.nanmax(np.abs(Q["cov"].to_numpy() - R[f"cov_{meth}"].to_numpy())))
        df_ = float(np.nanmax(np.abs(Q["cov_family"].to_numpy() - R[f"covfam_{meth}"].to_numpy())))
        if not (dw < 1e-9 and dc < 1e-12 and df_ < 1e-12):
            raise AssertionError(f"re-created {meth} draws differ from fewshot8 (max |dwidth| {dw:.3g}, "
                                 f"max |dcov| {dc:.3g}); the label sets are not identical")
        res[meth] = dict(n_draws=int(len(Q)), max_abs_dwidth=dw, max_abs_dcov=max(dc, df_))
    return res


def fewshot_by_m(B: pd.DataFrame, sd_y: float, floor: float) -> pd.DataFrame:
    rows = []
    for (scheme, m), G in B.groupby(["scheme", "m"]):
        r = dict(scheme=scheme, m=int(m), n_draws=len(G), sd_y=sd_y)
        for k in ("split", "offset_full", "nomodel"):
            c = G[f"cov_{k}"]
            r.update({f"cov_mean_{k}": c.mean(), f"cov_q05_{k}": c.quantile(0.05), f"cov_q50_{k}": c.quantile(0.5),
                      f"cov_q95_{k}": c.quantile(0.95), f"p_below_{k}": float(np.mean(c < floor)),
                      f"width_sd_{k}": G[f"width_{k}"].median() / sd_y})
        r["ratio_offset_to_nomodel"] = G.width_offset_full.median() / G.width_nomodel.median()
        r["ratio_split_to_nomodel"] = G.width_split.median() / G.width_nomodel.median()
        r["ratio_offset_to_nomodel_paired_median"] = float(np.median(G.width_offset_full / G.width_nomodel))
        r["share_draws_model_narrower"] = float(np.mean(G.width_offset_full < G.width_nomodel))
        rows.append(r)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- 2. inner-model SD
def group_level_scores(L: pd.DataFrame, B: int, n_perm: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    from .an8 import boot_spearman, perm_spearman_mc
    G = L.groupby("group").agg(n=("y", "size"), coverage_r1=("miscovered", lambda v: 1 - v.mean()),
                               committee_median=("committee", "median"), width_median=("width", "median"),
                               knn_median=("knn", "median")).reset_index()
    G["committee_rank"] = G.committee_median.rank(ascending=False, method="min").astype(int)
    rows = []
    for k, s in enumerate(("committee_median", "width_median", "knn_median")):
        rho, lo, hi, nb = boot_spearman(G.coverage_r1, G[s], B, seed + 10 * k)
        _, p = perm_spearman_mc(G.coverage_r1.to_numpy(), G[s].to_numpy(), n_perm, seed + 10 * k + 1)
        rows.append(dict(score=s, rho=rho, ci_lo=lo, ci_hi=hi, n_boot=nb, p_perm=p, n_groups=len(G)))
    return G, pd.DataFrame(rows)


def abstention_in_worst(L: pd.DataFrame, worst: list[int], fracs) -> pd.DataFrame:
    rows = []
    s = L.committee.to_numpy(float)
    cov = ~L.miscovered.to_numpy(bool)
    grp = L.group.to_numpy(int)
    for f in fracs:
        thr = np.inf if f <= 0 else float(np.quantile(s, 1 - f))
        keep = s <= thr if f > 0 else np.ones(s.size, bool)
        rows.append(dict(abstain=f, group="ALL", n=int(s.size), share_abstained=float(1 - keep.mean()),
                         retained_coverage=float(cov[keep].mean())))
        for g in worst:
            m = grp == g
            k = keep & m
            rows.append(dict(abstain=f, group=str(g), n=int(m.sum()), share_abstained=float(1 - k.sum() / m.sum()),
                             retained_coverage=float(cov[k].mean()) if k.any() else np.nan))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- 3. oxides
def perovskite_table(ctx, target: str, min_rows: int) -> tuple[pd.DataFrame, dict]:
    """CV+ coverage (alpha primary) of oxides split into ABO3 and other oxides, mean over repeats."""
    a = float(ctx.m["alpha_primary"])
    ox = ctx.F[ctx.m8["oxide"]["rule"]].to_numpy(bool)
    forms = ctx.formula
    cache = {}
    n_out = 2 * ctx.m["outer_folds_per_half"]
    rows = []
    for design, rep in (("primary", 1), ("rep2", 2)):
        if design == "rep2" and "rep2" not in ctx.m["blocks"][target]:
            continue
        for o in range(n_out):
            U = ctx.unit(target, design, rep, o)
            y = U["y_test"].astype(float)
            lo, hi = U[f"cvplus__a{a:.2f}__lo"].astype(float), U[f"cvplus__a{a:.2f}__hi"].astype(float)
            cov = (y >= lo) & (y <= hi)
            tr = U["test_rows"]
            is_ox = ox[tr]
            abo = np.array([cache.setdefault(f, is_abo3(f)) if io_ else False for f, io_ in zip(forms[tr], is_ox)])
            for g in np.unique(U["test_grp"]):
                m = U["test_grp"] == g
                for cls, sel in (("abo3", m & is_ox & abo), ("other_oxide", m & is_ox & ~abo), ("non_oxide", m & ~is_ox)):
                    rows.append(dict(rep=rep, group=int(g), cls=cls, n=int(sel.sum()), covered=int(cov[sel].sum())))
    D = pd.DataFrame(rows)
    per_rep = D.groupby(["group", "cls", "rep"]).agg(n=("n", "sum"), covered=("covered", "sum")).reset_index()
    per_rep["coverage"] = per_rep.covered / per_rep.n.where(per_rep.n > 0)
    T = per_rep.groupby(["group", "cls"]).agg(n=("n", "first"), coverage=("coverage", "mean")).reset_index()
    W = T.pivot(index="group", columns="cls", values=["n", "coverage"])
    W.columns = [f"{a_}_{b_}" for a_, b_ in W.columns]
    W = W.reset_index()
    for c in ("n_abo3", "n_other_oxide", "n_non_oxide"):
        if c not in W:
            W[c] = 0
        W[c] = W[c].fillna(0).astype(int)
    W.insert(1, "symbol", W.group.map(lambda g: ctx.sgsym.get(int(g), str(g))))
    W["n_oxide"] = W.n_abo3 + W.n_other_oxide
    sel = W[W.n_oxide >= min_rows].copy()
    sel["coverage_oxide"] = ((sel.coverage_abo3.fillna(0) * sel.n_abo3 + sel.coverage_other_oxide.fillna(0) * sel.n_other_oxide)
                             / sel.n_oxide)
    sel = sel.sort_values("coverage_oxide").reset_index(drop=True)
    tot = per_rep.groupby(["cls", "rep"]).agg(n=("n", "sum"), covered=("covered", "sum")).reset_index()
    tot["coverage"] = tot.covered / tot.n
    pooled = tot.groupby("cls").agg(n=("n", "first"), coverage=("coverage", "mean")).reset_index()
    summ = {r.cls: dict(n=int(r.n), coverage=float(r.coverage)) for r in pooled.itertuples()}
    return sel, summ


# ----------------------------------------------------------------------------- 4. degeneracy
def degeneracy_table(ctx, alphas) -> pd.DataFrame:
    rows = []
    n_out = 2 * ctx.m["outer_folds_per_half"]
    for t in ctx.m["targets"]:
        Ks = []
        for o in range(n_out):
            U = ctx.unit(t, "primary", 1, o)
            Ks.append(int(np.unique(U["cal_grp"]).size))
        for a in alphas:
            st = sorted({atom_status(K, a) for K in Ks})
            for meth in ("hierarchical jackknife+", "group-weighted CV+"):
                rows.append(dict(target=t, method=meth, alpha=a, K_min=min(Ks), K_max=max(Ks),
                                 atom=1.0 / (min(Ks) + 1), status="; ".join(st)))
        h3 = ctx.m.get("h3", {}).get("n_cal_groups", {})
        if t in h3 and "h3" in ctx.m["blocks"][t]:
            K1 = int(h3[t])
            for a in alphas:
                rows.append(dict(target=t, method="H3 HCP", alpha=a, K_min=K1, K_max=K1, atom=1.0 / (K1 + 1),
                                 status=atom_status(K1, a)))
                rows.append(dict(target=t, method="H3 single subsampling", alpha=a, K_min=K1, K_max=K1,
                                 atom=np.nan, status=single_subsample_status(K1, a)))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- driver
def ph8(cfg) -> dict:
    from .fit7 import MOCK
    p = dict(cfg.get("posthoc8", {}))
    p.setdefault("out_subdir", "posthoc")
    p.setdefault("bootstrap", 2000)
    p.setdefault("n_perm", 100000)
    p.setdefault("abstain", [0.0, 0.1, 0.25, 0.5])
    p.setdefault("coverage_floor", 0.80)
    if MOCK:
        p.update(bootstrap=200, n_perm=2000)
    return p


def run_posthoc8(cfg):
    from .run8 import DETECT_TARGETS, Ctx8, family_from_unit
    t0 = time.time()
    ctx = Ctx8(cfg)
    P = ph8(cfg)
    rep = ctx.d8["report"]
    out = os.path.join(ctx.d8["root"], P["out_subdir"])
    os.makedirs(out, exist_ok=True)
    a = float(ctx.m["alpha_primary"])
    floor = float(P["coverage_floor"])
    f8 = ctx.m8["fewshot"]
    seed = int(ctx.m8["seed"])
    W = pd.read_csv(os.path.join(ctx.d8["fewshot"], "worst_groups.csv"))
    pre_json = os.path.join(ctx.d8["phase0"], "PREREGISTRATION_phase6_v8.json")
    inputs = {pre_json: sha256_file(pre_json), os.path.join(ctx.d8["fewshot"], "worst_groups.csv"):
              sha256_file(os.path.join(ctx.d8["fewshot"], "worst_groups.csv"))}
    with open(pre_json) as fh:
        reg_w = {(t, int(w["group"])) for t, v in json.load(fh)["worst_groups"].items() for w in v}
    if reg_w and reg_w != {(str(r.target), int(r.group)) for r in W.itertuples()}:
        raise AssertionError("worst_groups.csv differs from the worst groups fixed in the v8 addendum")
    status, checks, lines = {}, {}, []
    print("Phase 6 Task 6.5 post hoc (v8.1; NOT pre-registered)", flush=True)

    # 1. few-shot: no-model baseline and coverage spread, paired with the saved draws
    BY = []
    for _, w in W.iterrows():
        t, g = str(w.target), int(w.group)
        fam, o = family_from_unit(ctx, t, g)
        dfile = os.path.join(ctx.d8["fewshot"], f"{t}_g{g}.csv.gz")
        inputs[dfile] = sha256_file(dfile)
        t1 = time.time()
        Bd = nomodel_draws(fam, f8["m"], f8["draws"], int(f8["min_test_groups"]), a, seed, f"{t}|{g}")
        checks[f"{t}_g{g}"] = verify_against_fewshot8(Bd, pd.read_csv(dfile))
        sd = float(np.std(fam.y, ddof=1))
        M = fewshot_by_m(Bd, sd, floor)
        M.insert(0, "symbol", ctx.sgsym.get(g, str(g)))
        M.insert(0, "group", g)
        M.insert(0, "target", t)
        M["unseen_coverage"] = float(w.primary_coverage)
        BY.append(M)
        Bd.insert(0, "group", g)
        Bd.insert(0, "target", t)
        status[f"fewshot_draws_{t}_g{g}.csv.gz"] = write_once(os.path.join(out, f"fewshot_draws_{t}_g{g}.csv.gz"),
                                                              _gz_bytes(Bd))
        r20 = M[(M.scheme == "one_per_group") & (M.m == 20)].iloc[0]
        print(f"  [posthoc8] {t} {ctx.sgsym.get(g, g)}: draws re-created and equal to fewshot8 "
              f"({checks[f'{t}_g{g}']['offset_full']['n_draws']} per method); m = 20 width / no-model "
              f"{r20.ratio_offset_to_nomodel:.2f} (offset), {r20.ratio_split_to_nomodel:.2f} (split) "
              f"({time.time() - t1:.0f} s)", flush=True)
    BYM = pd.concat(BY, ignore_index=True)
    status["fewshot_by_m.csv"] = write_once(os.path.join(out, "fewshot_by_m.csv"), _bytes_csv(BYM))

    # 2. inner-model SD at group level and abstention inside the worst groups
    SDT, ABS, GL = [], [], []
    for k, t in enumerate(DETECT_TARGETS):
        f = os.path.join(rep, f"{t}_structures.csv.gz")
        inputs[f] = sha256_file(f)
        L = pd.read_csv(f)
        G, S = group_level_scores(L, int(P["bootstrap"]), int(P["n_perm"]), seed + 900 + 100 * k)
        S.insert(0, "target", t)
        G.insert(0, "target", t)
        G.insert(2, "symbol", G.group.map(lambda x: ctx.sgsym.get(int(x), str(x))))
        worst = W[W.target == t].group.astype(int).tolist()
        A = abstention_in_worst(L, worst, P["abstain"])
        A.insert(0, "target", t)
        SDT.append(S)
        GL.append(G)
        ABS.append(A)
    SDT, GL, ABS = (pd.concat(x, ignore_index=True) for x in (SDT, GL, ABS))
    status["committee_spearman.csv"] = write_once(os.path.join(out, "committee_spearman.csv"), _bytes_csv(SDT))
    status["committee_groups.csv"] = write_once(os.path.join(out, "committee_groups.csv"), _bytes_csv(GL))
    status["committee_abstention_worst.csv"] = write_once(os.path.join(out, "committee_abstention_worst.csv"),
                                                          _bytes_csv(ABS))

    # 3. oxides: ABO3 vs other oxides
    PV, PVS = [], {}
    for t in ctx.m["targets"]:
        T, s = perovskite_table(ctx, t, int(ctx.m8["oxide"]["min_rows"]))
        T.insert(0, "target", t)
        PV.append(T)
        PVS[t] = s
    PV = pd.concat(PV, ignore_index=True)
    status["oxides_abo3.csv"] = write_once(os.path.join(out, "oxides_abo3.csv"), _bytes_csv(PV))

    # 4. degeneracy flags
    DG = degeneracy_table(ctx, [float(x) for x in ctx.m8["hjplus"]["alphas"]])
    status["degeneracy.csv"] = write_once(os.path.join(out, "degeneracy.csv"), _bytes_csv(DG))

    # summary
    def md(df):
        return df.to_markdown(index=False, floatfmt=".3f")
    show = BYM[BYM.scheme == "one_per_group"][["target", "symbol", "unseen_coverage", "m", "cov_mean_split",
                                               "cov_q05_split", "cov_q95_split", "p_below_split", "width_sd_split",
                                               "width_sd_offset_full", "width_sd_nomodel", "ratio_offset_to_nomodel",
                                               "ratio_split_to_nomodel", "share_draws_model_narrower",
                                               "cov_mean_nomodel"]]
    absw = ABS[ABS.group != "ALL"].copy()
    absw["symbol"] = absw.group.map(lambda x: ctx.sgsym.get(int(x), str(x)))
    lines = ["# Phase 6 Task 6.5 post hoc analyses (v8.1; NOT pre-registered)", "",
             "Designed after the v8 results were seen (2026-11-19). Nothing here changes a pre-registered decision.", "",
             "## 1. Few-shot intervals against a no-model interval from the same labels", "",
             f"Label sets re-created with the registered fewshot8 seeds and checked draw by draw against the saved draws "
             f"(model widths and coverages identical). No-model interval = offset-corrected full conformal on the "
             f"labelled y values (centred on their mean; the posthoc5 definition). alpha = {a}; "
             f"{f8['draws']} draws per cell; one structure per composition. ratio < 1: the model interval is narrower.", "",
             md(show), "",
             "## 2. Inner-model SD (MatFold's ensemble SD) at group level", "",
             "Spearman of repeat-1 group coverage with the group median of each score; percentile bootstrap over "
             f"groups ({P['bootstrap']}) and Monte Carlo permutation p ({P['n_perm']}).", "", md(SDT), "",
             "Coverage retained inside the three worst groups when the structures with the highest SD (pooled "
             "threshold over all held-out structures of the target) are abstained:", "",
             md(absw[["target", "symbol", "abstain", "n", "share_abstained", "retained_coverage"]]), "",
             "## 3. Oxides: ABO3 compositions (composition rule only) and other oxides", "",
             "Pooled over all eligible held-out groups (mean over repeats): " +
             "; ".join(f"{t}: ABO3 {v['abo3']['coverage']:.3f} (n {v['abo3']['n']}), other oxides "
                       f"{v['other_oxide']['coverage']:.3f} (n {v['other_oxide']['n']}), non-oxides "
                       f"{v['non_oxide']['coverage']:.3f} (n {v['non_oxide']['n']})"
                       for t, v in PVS.items() if {"abo3", "other_oxide", "non_oxide"} <= set(v)) + ".", "",
             md(PV.groupby("target").head(8)[["target", "group", "symbol", "n_abo3", "coverage_abo3", "n_other_oxide",
                                  "coverage_other_oxide", "coverage_oxide"]]), "",
             "## 4. Degenerate group-level intervals", "",
             "atom = 1/(K+1) at +infinity; unbounded if atom > alpha, upper end = largest calibration term if atom = "
             "alpha. Single subsampling: rank ceil((K+1)(1-alpha)) vs K.", "", md(DG), ""]
    txt = "\n".join(lines)
    status["POSTHOC8_SUMMARY.md"] = write_once(os.path.join(out, "POSTHOC8_SUMMARY.md"), txt.encode())
    meta_path = os.path.join(out, "POSTHOC8.json")
    meta = dict(label="POST HOC, not pre-registered (v8.1)", created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                code_hash=code_hash(), inputs_sha256=inputs, verification_against_fewshot8=checks,
                abo3_pooled=PVS, files=status, minutes=(time.time() - t0) / 60)
    if not os.path.exists(meta_path):
        atomic_write_bytes(meta_path, json.dumps(meta, indent=2, default=str).encode())
    print("\n" + txt)
    print(f"[posthoc8] outputs in {out}: " + ", ".join(f"{k} ({v})" for k, v in status.items()))
    print(f"[posthoc8] done ({(time.time() - t0) / 60:.1f} min). Next: figures81")
    return dict(by_m=BYM, committee=SDT, abstention=ABS, oxides=PV, degeneracy=DG, checks=checks)
