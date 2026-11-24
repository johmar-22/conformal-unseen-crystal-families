"""Phase 5 post hoc analyses (NOT pre-registered; designed after the Phase 5 results were seen).

Two questions raised by the Phase 5 results, answered from the saved Drive files only:

1. Why do clustered labels harm some families and not others?  One-way random-effects intraclass
   correlation ICC(1) of the residuals z = y - mu by composition group (ANOVA estimator for
   unequal group sizes, n0 of Donner 1986; Shrout & Fleiss 1979), and of the absolute deviations
   |z - mean z| (clustering of the error SIZE). Compared across the 14 family-target pairs with
   the pre-registered clustered dispersion ratio (Spearman, with a permutation p-value).

2. Is the recalibrated model interval more useful than no model at all?  A no-model interval
   built from the SAME labels: full conformal for y = b + e (the offset-corrected interval with
   mu = 0, i.e. centred on the mean of the labelled y values). It is valid for the same reason as
   the model interval. The label sets are re-created with the registered Phase 5 seeds and
   checked draw by draw against the saved Phase 5 draws (the model widths and coverages must be
   identical), so the comparison is paired.

Inputs (read only): phase05_v1/phase0/PREREGISTRATION_phase5.json (frozen config: seeds, grids),
phase05_v1/inputs/train.csv, phase02_v1/loso/*.npz, phase05_v1/draws/*.csv.gz,
phase05_v1/report/phase5_{dispersion,summary}.csv.
Outputs: phase05_v1/posthoc/ (new folder; nothing existing is touched).
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
import pandas as pd

from .config import (Paths, atomic_write_bytes, atomic_write_csv, atomic_write_json, atomic_write_text,
                     code_hash, sha256_file)
from .data import LATTICE_NAME, id_table
from .fewshot import Family, draw_one_per_group, group_weighted, hull_width, in_set, offset_full_set, split_set
from .phase5 import _gz, _seed, config_label, load_family, parse_config


# ----------------------------------------------------------------------------- statistics
def icc1(x: np.ndarray, groups: np.ndarray) -> float:
    """One-way random-effects ICC(1), ANOVA estimator with n0 for unequal group sizes:
    ICC = (MSB - MSW) / (MSB + (n0 - 1) MSW), n0 = (N - sum n_g^2 / N) / (G - 1)."""
    x = np.asarray(x, float)
    _, inv = np.unique(groups, return_inverse=True)
    n = np.bincount(inv).astype(float)
    G, N = n.size, x.size
    if G < 2 or N <= G:
        return float("nan")
    means = np.bincount(inv, x) / n
    ssb = float(np.sum(n * (means - x.mean()) ** 2))
    ssw = float(np.sum((x - means[inv]) ** 2))
    msb, msw = ssb / (G - 1), ssw / (N - G)
    n0 = (N - np.sum(n ** 2) / N) / (G - 1)
    den = msb + (n0 - 1) * msw
    return float((msb - msw) / den) if den > 0 else float("nan")


def spearman_perm(x, y, n_perm: int, seed: int):
    """Spearman rho with a two-sided permutation p-value (add-one estimate)."""
    from scipy.stats import rankdata, spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    rho = float(spearmanr(x, y).statistic)
    rx, ry = rankdata(x), rankdata(y)
    rng = np.random.default_rng(seed)
    obs = abs(np.corrcoef(rx, ry)[0, 1])
    hits = sum(abs(np.corrcoef(rx, rng.permutation(ry))[0, 1]) >= obs - 1e-12 for _ in range(n_perm))
    return rho, (hits + 1) / (n_perm + 1), float(spearmanr(x, y).pvalue)


# ----------------------------------------------------------------------------- paths
def phase5_paths(cfg: dict) -> tuple[Paths, dict]:
    p = cfg["paths"]
    root = p["drive_root"]
    if not os.path.isdir(root):
        raise FileNotFoundError(f"Drive folder not found: {root}\nMount Drive first.")
    p5 = os.path.join(root, p["phase5_dir"])
    js = os.path.join(p5, "phase0", "PREREGISTRATION_phase5.json")
    if not os.path.exists(js):
        raise FileNotFoundError(f"Phase 5 pre-registration not found: {js}")
    with open(js) as f:
        reg = json.load(f)
    c5 = reg["config"]
    paths = Paths(drive_root=root, out=p5, cache=p["local_cache"], phase2=os.path.join(root, p["phase2_dir"]))
    os.makedirs(p["local_cache"], exist_ok=True)
    return paths, dict(config=c5, prereg_json=js, prereg_sha256=sha256_file(js), prereg_time=reg["time_utc"])


# ----------------------------------------------------------------------------- draws
def baseline_draws(fam: Family, c5: dict, design: str, g: int, target: str) -> pd.DataFrame:
    """Re-create the Phase 5 one-per-composition label sets (same seeds, same calls, same order as
    phase5.run_family) and evaluate the model intervals and the no-model interval on each."""
    f = c5["fewshot"]
    alpha = float(c5["alpha_primary"])
    D, mtg = int(f["draws"]), int(f["min_test_groups"])
    rows = []
    for m in f["m_one_per_group"]:
        rng = np.random.default_rng(_seed(c5["seed_base"], design, g, target, "one_per_group", m, False))
        for d in range(D):
            r = draw_one_per_group(fam, m, rng, mtg)
            if r is None:
                break
            cal, new, tg = r
            out = dict(m=int(m), draw=d)
            for name, vals, fn in (("split", fam.z, split_set), ("offset", fam.z, offset_full_set),
                                   ("nomodel", fam.y, offset_full_set)):
                S = fn(vals[cal], alpha)
                out[f"width_{name}"] = hull_width(S)
                out[f"cov_{name}"] = group_weighted(in_set(vals[new], S), new, fam, tg)
            rows.append(out)
    return pd.DataFrame(rows)


def verify_against_phase5(B: pd.DataFrame, draws_file: str) -> dict:
    """The model widths and coverages of the re-created draws must equal the saved Phase 5 draws."""
    P = pd.read_csv(draws_file)
    P = P[(P.scheme == "one_per_group") & (P.null.astype(str) == "False")]
    res = {}
    for meth, col in (("split", "split"), ("offset_full", "offset")):
        Q = P[P.method == meth].set_index(["m", "draw"]).sort_index()
        R = B.set_index(["m", "draw"]).sort_index()
        if len(Q) != len(R) or not Q.index.equals(R.index):
            raise AssertionError(f"{os.path.basename(draws_file)}: draw index differs from Phase 5 ({meth})")
        dw = float(np.nanmax(np.abs(Q["width"].to_numpy() - R[f"width_{col}"].to_numpy())))
        dc = float(np.nanmax(np.abs(Q["cov"].to_numpy() - R[f"cov_{col}"].to_numpy())))
        if not (dw < 1e-9 and dc < 1e-12):
            raise AssertionError(f"{os.path.basename(draws_file)}: re-created {meth} draws differ from Phase 5 "
                                 f"(max |dwidth| {dw:.3g}, max |dcov| {dc:.3g}); label sets are not identical")
        res[meth] = dict(n_draws=int(len(Q)), max_abs_dwidth=dw, max_abs_dcov=dc)
    return res


# ----------------------------------------------------------------------------- main
def run_posthoc5(cfg: dict):
    t0 = time.time()
    paths, reg = phase5_paths(cfg)
    c5 = reg["config"]
    ph = cfg["posthoc5"]
    floor = float(c5["decisions"]["coverage_floor"])
    m_main = int(c5["decisions"]["h5_m"])
    outdir = os.path.join(paths.out, ph.get("out_subdir", "posthoc"))
    os.makedirs(outdir, exist_ok=True)
    rep = os.path.join(paths.out, "report")
    need = [os.path.join(rep, "phase5_dispersion.csv"), os.path.join(rep, "phase5_summary.csv")]
    for f in need:
        if not os.path.exists(f):
            raise FileNotFoundError(f"{f} not found (run report5 of cpaudit_v4 first)")
    inputs = {f: sha256_file(f) for f in need + [reg["prereg_json"]]}
    print(f"Phase 5 post hoc (not pre-registered) | frozen Phase 5 config from {reg['prereg_json']} "
          f"({reg['prereg_time']})", flush=True)
    tab = id_table(c5, paths)
    inputs[os.path.join(paths.inputs, "train.csv")] = sha256_file(os.path.join(paths.inputs, "train.csv"))
    Dsp = pd.read_csv(need[0])
    S5 = pd.read_csv(need[1])

    fam_rows, by_m, checks, draw_tabs = [], [], {}, []
    confs = [parse_config(c) for c in c5["configs"]]
    for design, g in confs:
        for target in ("bandgap", "formation"):
            npz = os.path.join(paths.phase2_loso, f"{target}_{design}{g}.npz")
            inputs[npz] = sha256_file(npz)
            fam = load_family(paths, tab, design, g, target)
            sd = float(np.std(fam.y, ddof=1))
            B = baseline_draws(fam, c5, design, g, target)
            dfile = os.path.join(paths.out, "draws", f"{target}_{design}{g}.csv.gz")
            if not os.path.exists(dfile):
                raise FileNotFoundError(f"Phase 5 draws not found: {dfile}")
            inputs[dfile] = sha256_file(dfile)
            checks[f"{target}_{design}{g}"] = verify_against_phase5(B, dfile)
            B.insert(0, "target", target)
            B.insert(0, "heldout", LATTICE_NAME[g])
            B.insert(0, "design", design)
            draw_tabs.append(B)
            for m, G in B.groupby("m"):
                wm, wb, ws = G.width_offset.median(), G.width_nomodel.median(), G.width_split.median()
                by_m.append(dict(design=design, heldout=LATTICE_NAME[g], target=target, m=int(m), n_draws=len(G),
                                 sd_y=sd, width_split_sd=ws / sd, width_offset_sd=wm / sd, width_nomodel_sd=wb / sd,
                                 ratio_offset_to_nomodel=wm / wb,
                                 ratio_paired_median=float(np.median(G.width_offset / G.width_nomodel)),
                                 cov_offset=G.cov_offset.mean(), cov_nomodel=G.cov_nomodel.mean(),
                                 p_below_nomodel=float(np.mean(G.cov_nomodel < floor))))
            dz = np.abs(fam.z - fam.z.mean())
            drow = Dsp[(Dsp.scheme == "clustered") & (Dsp.design == design) & (Dsp.sg == g) &
                       (Dsp.target == target) & (Dsp.m == m_main)]
            srow = S5[(S5.scheme == "clustered") & (S5.null.astype(str) == "False") & (S5.design == design) &
                      (S5.sg == g) & (S5.target == target) & (S5.m == m_main) & (S5.method == "offset_full")]
            fam_rows.append(dict(design=design, heldout=LATTICE_NAME[g], target=target, family=config_label(design, g),
                                 N=fam.N, G=fam.G, sd_y=sd,
                                 icc_resid=icc1(fam.z, fam.groups), icc_absdev=icc1(dz, fam.groups),
                                 rho_clustered=float(drow.rho.iloc[0]) if len(drow) else np.nan,
                                 cov_clustered_split=float(drow.mean_real.iloc[0]) if len(drow) else np.nan,
                                 cov_clustered_offset=float(srow.mean_cov.iloc[0]) if len(srow) else np.nan))
            print(f"  {target:9s} {config_label(design, g):32s} ICC {fam_rows[-1]['icc_resid']:.3f} | "
                  f"draws match Phase 5 ({checks[f'{target}_{design}{g}']['offset_full']['n_draws']} per method)",
                  flush=True)

    M = pd.DataFrame(by_m)
    F = pd.DataFrame(fam_rows)
    F = F.merge(M[M.m == m_main][["design", "heldout", "target", "width_offset_sd", "width_nomodel_sd",
                                  "ratio_offset_to_nomodel", "ratio_paired_median", "cov_offset", "cov_nomodel",
                                  "p_below_nomodel"]], on=["design", "heldout", "target"], how="left")
    ok = F.rho_clustered.notna()
    rho, p_perm, p_asym = spearman_perm(F.icc_resid[ok], F.rho_clustered[ok], int(ph.get("n_perm", 100000)),
                                        int(c5["seed_base"]) + 7000)
    rho2, p2_perm, _ = spearman_perm(F.icc_absdev[ok], F.rho_clustered[ok], int(ph.get("n_perm", 100000)),
                                     int(c5["seed_base"]) + 7001)
    stats = dict(spearman_icc_resid_vs_rho_clustered=dict(rho=rho, p_permutation=p_perm, p_asymptotic=p_asym,
                                                         n=int(ok.sum())),
                 spearman_icc_absdev_vs_rho_clustered=dict(rho=rho2, p_permutation=p2_perm, n=int(ok.sum())),
                 m=m_main,
                 ratio_offset_to_nomodel_range={t: [float(F[F.target == t].ratio_offset_to_nomodel.min()),
                                                    float(F[F.target == t].ratio_offset_to_nomodel.max())]
                                                for t in ("bandgap", "formation")},
                 width_nomodel_sd_range=[float(F.width_nomodel_sd.min()), float(F.width_nomodel_sd.max())],
                 cov_nomodel_range=[float(F.cov_nomodel.min()), float(F.cov_nomodel.max())])

    atomic_write_csv(os.path.join(outdir, "phase5_posthoc_family.csv"), F)
    atomic_write_csv(os.path.join(outdir, "phase5_posthoc_by_m.csv"), M)
    atomic_write_bytes(os.path.join(outdir, "phase5_posthoc_draws.csv.gz"), _gz(pd.concat(draw_tabs, ignore_index=True)))
    meta = dict(label="POST HOC, not pre-registered", created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                code_hash=code_hash(), phase5_preregistration=dict(path=reg["prereg_json"], sha256=reg["prereg_sha256"],
                                                                   time_utc=reg["prereg_time"]),
                inputs_sha256=inputs, verification_against_phase5_draws=checks, statistics=stats,
                minutes=(time.time() - t0) / 60)
    atomic_write_json(os.path.join(outdir, "phase5_posthoc.json"), meta)

    show = F[["family", "target", "icc_resid", "icc_absdev", "rho_clustered", "cov_clustered_split",
              "cov_clustered_offset", "width_offset_sd", "width_nomodel_sd", "ratio_offset_to_nomodel",
              "cov_offset", "cov_nomodel"]]
    txt = "\n".join([
        "# Phase 5 post hoc analyses (NOT pre-registered)", "",
        f"Designed after the Phase 5 results were seen. m = {m_main} labels, one per composition, alpha = "
        f"{c5['alpha_primary']}; {c5['fewshot']['draws']} draws per cell, the SAME label sets as Phase 5 "
        "(re-created with the registered seeds and checked draw by draw against the saved draws).", "",
        "* icc_resid: ICC(1) of the residuals by composition group; icc_absdev: ICC(1) of |z - mean z|.",
        "* rho_clustered, cov_clustered_split: pre-registered H5c quantities (split, clustered labels).",
        "* cov_clustered_offset: mean coverage of new compositions, offset-corrected interval, clustered labels.",
        "* width_*_sd: median interval width / family SD. nomodel = full conformal on the labelled y values "
        "(no model). ratio = median width offset / median width no-model.", "",
        show.round(3).to_markdown(index=False), "",
        f"Spearman ICC(resid) vs clustered dispersion ratio: rho = {rho:.2f}, permutation p = {p_perm:.4f} "
        f"(asymptotic {p_asym:.4f}), n = {int(ok.sum())} family-target pairs.",
        f"Spearman ICC(|dev|) vs clustered dispersion ratio: rho = {rho2:.2f}, permutation p = {p2_perm:.4f}.", "",
        "Verification: model widths and coverages of every re-created draw equal the Phase 5 draws "
        f"(max |dwidth| {max(v['offset_full']['max_abs_dwidth'] for v in checks.values()):.2g}, "
        f"max |dcov| {max(v['offset_full']['max_abs_dcov'] for v in checks.values()):.2g}).", ""])
    atomic_write_text(os.path.join(outdir, "POSTHOC_SUMMARY.md"), txt)
    print("\n" + txt)
    print(f"Post hoc outputs written to {outdir} ({(time.time() - t0) / 60:.1f} min)")
    return F, M, meta
