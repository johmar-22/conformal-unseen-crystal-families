"""report7: Task 6.4 results tables and the pre-registered decisions that need only Task 6.4
outputs (G5, H1-J, G3-J, H3-J). H4-J and H5-J are evaluated in v8 (detection scores, few-shot)."""
from __future__ import annotations

import glob
import json
import os

import numpy as np
import pandas as pd

from .config import atomic_write_csv, atomic_write_json, atomic_write_text


def _load(pattern):
    fs = sorted(glob.glob(pattern))
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def _fmt(x, k=3):
    return "nan" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{k}f}"


def run_report7(cfg):
    from .fit7 import BLOCK_UNITS, Ctx
    ctx = Ctx(cfg)
    m, d = ctx.m, ctx.d
    a_p = round(float(m["alpha_primary"]), 3)
    dec = m["decisions"]
    B = _load(os.path.join(d["unitsB"], "*.csv"))
    A = _load(os.path.join(d["unitsA"], "*.csv"))
    H = _load(os.path.join(d["unitsH3"], "*.csv"))
    lines = ["# Phase 6 Task 6.4 report (JARVIS-DFT)", ""]
    decisions, tables = {}, {}
    # completeness
    lines += ["## Units", ""]
    for t in m["targets"]:
        D = ctx.design(t)
        for b in m["blocks"][t]:
            units = BLOCK_UNITS[b](D, m)
            sub = {"regimeA": "unitsA", "h3": "unitsH3"}.get(b, "unitsB")
            done = 0
            for kind, design, rep, k in units:
                stem = (f"{t}_{design}_r{rep}_o{k:02d}" if kind == "B" else
                        f"{t}_A_s{k:02d}" if kind == "A" else f"{t}_H3_p{k:02d}")
                done += os.path.exists(os.path.join(d[sub], stem + ".csv"))
            lines.append(f"* {t} {b}: {done}/{len(units)} units")
    lines.append("")
    for t in m["targets"]:
        D = ctx.design(t)
        elig = set(int(g) for g in D["eligible"])
        lines += [f"## {t}", ""]
        dt = {}
        if len(B):
            b = B[(B.target == t) & (B.group != "ALL")].copy()
            b["alpha"] = b.alpha.round(3)
            prim = b[(b.design.isin(["primary", "rep2"])) & (b.method == "cvplus") & (b.alpha == a_p)]
            if len(prim):
                pc = prim.groupby("group").agg(coverage=("coverage", "mean"), n_repeats=("rep", "nunique"),
                                               n=("n", "first"), mae=("mae", "mean"), bias=("bias", "mean"),
                                               width_sd=("width_median_sd", "mean")).reset_index()
                pc["group"] = pc.group.astype(int)
                pc = pc[pc.group.isin(elig)].sort_values("coverage")
                tables[f"{t}_primary_by_group"] = pc
                share_cov = float((pc.coverage >= dec["g5"]["coverage"]).mean())
                dt["G5"] = dict(share_covered=share_cov, holds=share_cov >= dec["g5"]["share"],
                                n_groups=int(len(pc)), share_below_0_90=float((pc.coverage < 0.90).mean()),
                                share_below_0_80=float((pc.coverage < 0.80).mean()),
                                median_coverage=float(pc.coverage.median()))
                lines += [f"Primary (CV+, space-group inner folds, alpha {a_p}; mean over repeats): {len(pc)} groups; "
                          f"median coverage {_fmt(pc.coverage.median())}; below 0.90: {dt['G5']['share_below_0_90']:.0%}; "
                          f"below 0.80: {dt['G5']['share_below_0_80']:.0%}. **G5** (>= {dec['g5']['share']:.0%} of groups "
                          f">= {dec['g5']['coverage']}): {'holds' if dt['G5']['holds'] else 'does not hold'} "
                          f"({share_cov:.0%}).", ""]
                lines += ["Worst 10 groups:", "", pc.head(10).round(3).to_markdown(index=False), ""]
            # method summary (repeat 1 designs)
            s = b[b.alpha == a_p].groupby(["design", "method"]).agg(
                median_coverage=("coverage", "median"), share_below_090=("coverage", lambda x: float((x < 0.9).mean())),
                share_below_080=("coverage", lambda x: float((x < 0.8).mean())),
                median_width_sd=("width_median_sd", "median"), groups=("group", "nunique")).reset_index()
            tables[f"{t}_method_summary"] = s
            lines += ["Per-group summary by design and method (alpha 0.10):", "", s.round(3).to_markdown(index=False), ""]
            # G3-J
            r1 = b[(b.design == "primary") & (b.rep == 1) & (b.alpha == a_p)]
            if len(r1) and t != "tbmbj":
                piv = r1.pivot_table(index="group", columns="method", values="coverage")
                fail = piv[piv["cvplus"] < dec["g3"]["failing"]]
                lifted = ((fail.get("ncvplus", -1) >= dec["g3"]["lifted"]) |
                          (fail.get("cqr_cvplus", -1) >= dec["g3"]["lifted"])) if len(fail) else pd.Series(dtype=bool)
                share = float(lifted.mean()) if len(fail) else np.nan
                dt["G3J"] = dict(n_failing=int(len(fail)), share_lifted=share,
                                 positive=bool(len(fail) and share >= dec["g3"]["share"]))
                lines += [f"**G3-J**: {len(fail)} groups below {dec['g3']['failing']} under CV+ (repeat 1); "
                          f"normalized or CQR lifts {('n/a' if not len(fail) else f'{share:.0%}')} to >= "
                          f"{dec['g3']['lifted']} -> {'positive' if dt['G3J']['positive'] else 'negative'}.", ""]
        if len(A):
            a = A[(A.target == t)].copy()
            a["alpha"] = a.alpha.round(3)
            a = a[a.alpha == a_p]
            if len(a):
                allr = a[(a.group == "ALL") & (a.method == "cvplus_random")]
                marg = float(allr.coverage.mean())
                g = a[(a.group != "ALL") & (a.group != "-1")].copy()
                g["group"] = g.group.astype(int)
                g = g[g.group.isin(elig)]
                pg = g.groupby(["group", "method"]).coverage.mean().unstack()
                low = float((pg["cvplus_random"] < dec["h1"]["low"]).mean())
                mon = float((pg["mondrian"] >= dec["h1"]["mondrian_coverage"]).mean())
                sup = (dec["h1"]["marginal_lo"] <= marg <= dec["h1"]["marginal_hi"] and low >= dec["h1"]["share_support"]
                       and mon >= dec["h1"]["mondrian_share"])
                verdict = "supported" if sup else ("heterogeneity absent" if low < dec["h1"]["share_absent"] else "not supported")
                dt["H1J"] = dict(marginal=marg, share_groups_below=low, share_mondrian_ok=mon, verdict=verdict)
                tables[f"{t}_regimeA_by_group"] = pg.reset_index()
                lines += [f"Regime A: marginal CV+ random coverage {marg:.3f}; groups below {dec['h1']['low']}: "
                          f"{low:.0%}; Mondrian >= {dec['h1']['mondrian_coverage']}: {mon:.0%}. **H1-J**: {verdict}.", ""]
                lines += [a[a.group == "ALL"].groupby("method").agg(coverage=("coverage", "mean"),
                          width=("width_median", "mean"), mae=("mae", "mean")).round(4).to_markdown(), ""]
        if len(H):
            h = H[(H.target == t)].copy()
            h["alpha"] = h.alpha.round(3)
            if len(h):
                gw = h.groupby(["method", "alpha", "partition", "reshuffle"]).coverage.mean().groupby(
                    ["method", "alpha"]).agg(["mean", "std"]).reset_index()
                cond = h[h.alpha == a_p].groupby("method").coverage.apply(lambda x: float((x < 0.8).mean()))
                p = gw[gw.alpha == a_p].set_index("method")["mean"]
                sup = all(p.get(k, 0) >= dec["h3"]["support"] for k in ("hcp", "single_subsample"))
                rej = any(p.get(k, 1) < dec["h3"]["reject"] for k in ("hcp", "single_subsample"))
                dt["H3J"] = dict(coverage=p.to_dict(), share_test_groups_below_080=cond.to_dict(),
                                 verdict="supported" if sup else ("rejected" if rej else "inconclusive"))
                tables[f"{t}_h3"] = gw
                lines += [f"H3 (valid group-level design): group-weighted coverage {json.dumps({k: round(v, 3) for k, v in p.items()})}; "
                          f"share of test groups with own coverage < 0.80: {json.dumps({k: round(v, 3) for k, v in cond.items()})}. "
                          f"**H3-J**: {dt['H3J']['verdict']}.", "", gw.round(4).to_markdown(index=False), ""]
        decisions[t] = dt
    atomic_write_json(os.path.join(d["report"], "phase6_decisions_v7.json"), decisions)
    for k, v in tables.items():
        atomic_write_csv(os.path.join(d["report"], f"{k}.csv"), v)
    txt = "\n".join(lines)
    atomic_write_text(os.path.join(d["report"], "PHASE6_TASK64_SUMMARY.md"), txt)
    print(txt)
    return decisions
