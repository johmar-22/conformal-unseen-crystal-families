"""Exploratory follow-up checks on the saved Phase 1 outputs (no model is refitted).

1. De-duplication sensitivity. Phase 1 removed 9 entries; the v2 rule removes more (see
   setup2/geometry). All removed entries share a composition group with the entry that is kept,
   so no duplicate ever crossed a train/test or proper/calibration boundary in Phase 1 (the splits
   are group-disjoint). The only effect is extra weight on a few repeated points. Here every
   interval is recomputed from the saved predictions and residuals with the v2-removed entries
   dropped from the test, calibration and CV+ residual sets (models unchanged), and, in a second
   variant, with the conflicting records also excluded from evaluation.
2. Cluster-aware Mondrian. Points inside a composition group are correlated, and a lattice has
   only about 10-13 calibration composition groups per split. Pooled per-point calibration is then
   only approximately valid. Dunn, Wasserman & Ramdas (JASA 118, 2491, 2023) show that drawing
   ONE calibration point per group restores exact validity for a random point of a new group
   (single subsampling), i.e. for group-weighted coverage. Its expected coverage is estimated
   here by averaging over `cluster_draws` random draws per split. If this method reaches nominal
   group-weighted coverage while pooled Mondrian falls short, the Phase 1 H1b shortfall is
   attributable to within-group correlation, not to a coding or calibration error.
3. Phase 1 results by Fd-3m sub-family (Task 2.3, in-distribution part), from the Phase 1 report.
These analyses were designed after Phase 1 results were seen; they are exploratory.
"""
from __future__ import annotations

import glob
import os
import re

import numpy as np
import pandas as pd

from . import conformal as cp
from .config import atomic_write_csv, atomic_write_text, get_paths
from .data import LATTICE_NAME, LATTICE_ORDER, build_dataset, id_table
from .setup2 import _phase1_removed

SPLIT_METHODS = ("split_global", "split_mondrian", "split_normalized_knn")


def intervals_from_npz(z, sg_cal, keep_te, keep_ca, keep_tr, alphas, a_p):
    """All Phase 1 interval methods, recomputed on the kept points (models unchanged)."""
    out = {}
    p_te, p_ca = z["pred_split_test"][keep_te], z["pred_split_cal"][keep_ca]
    y_ca = z["y_cal"][keep_ca]
    s_ca = np.abs(y_ca - p_ca)
    sg_te = z["sg_test"][keep_te]
    sc = sg_cal[keep_ca]
    sig_ca, sig_te = z["sigma_knn_cal"][keep_ca], z["sigma_knn_test"][keep_te]
    for a in alphas:
        out[("split_global", a)] = cp.split_interval(p_te, s_ca, a)[:2]
        out[("split_mondrian", a)] = cp.mondrian_interval(p_te, sg_te, s_ca, sc, a)[:2]
        out[("split_normalized_knn", a)] = cp.normalized_interval(p_te, sig_te, y_ca, p_ca, sig_ca, a)[:2]
        for fs, m in (("cvr", "cvplus_random"), ("cvl", "cvplus_lattice")):
            if f"{fs}_fold_pred_test" in z:
                out[(m, a)] = cp.cvplus_interval(z[f"{fs}_fold_pred_test"][:, keep_te],
                                                 z[f"{fs}_fold_of_train"][keep_tr], z[f"{fs}_resid"][keep_tr], a)
    if "cqr_lo_test" in z:
        out[("cqr", a_p)] = cp.cqr_interval(z["cqr_lo_test"][keep_te], z["cqr_hi_test"][keep_te], y_ca,
                                            z["cqr_lo_cal"][keep_ca], z["cqr_hi_cal"][keep_ca], a_p)[:2]
    return out


def _same(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return bool(np.array_equal(np.isinf(a), np.isinf(b)) and np.allclose(a[np.isfinite(a)], b[np.isfinite(b)],
                                                                          rtol=1e-10, atol=1e-12))


def group_weighted(covered, groups):
    s = pd.Series(covered.astype(float)).groupby(groups).mean()
    return float(s.mean())


def run_post1(cfg):
    paths = get_paths(cfg)
    tab = id_table(cfg, paths)
    ds = build_dataset(cfg, paths)
    removed = set(ds.dedup["removed_ids"])
    conflict = set(ds.dedup["conflicting_ids"])
    files = sorted(glob.glob(os.path.join(paths.phase1, "runs", "r*_*_hpo_once.npz")))
    if not files:
        raise RuntimeError(f"No Phase 1 run files found in {os.path.join(paths.phase1, 'runs')}")
    alphas = [float(a) for a in cfg["alpha"]]
    a_p = float(cfg["alpha_primary"])
    B = int(cfg["post1"]["cluster_draws"])
    rows_ded, rows_clu = [], []
    n_checked = 0
    for f in files:
        m = re.search(r"r(\d+)_(\w+?)_hpo_once\.npz$", os.path.basename(f))
        r, target = int(m.group(1)), m.group(2)
        z = dict(np.load(f))
        te_id, ca_id, tr_id = z["test_id"], z["cal_id"], z["train_id"]
        sg_cal = tab.loc[ca_id, "sg"].to_numpy(int)
        g_te, g_ca = tab.loc[te_id, "group"].to_numpy(), tab.loc[ca_id, "group"].to_numpy()
        y_te = z["y_test"]
        all_te, all_ca, all_tr = np.ones(len(te_id), bool), np.ones(len(ca_id), bool), np.ones(len(tr_id), bool)
        # (0) recomputation reproduces the stored Phase 1 intervals exactly
        base = intervals_from_npz(z, sg_cal, all_te, all_ca, all_tr, alphas, a_p)
        for (meth, a), (lo, hi) in base.items():
            key = f"{meth}__a{a:.2f}__"
            if key + "lo" in z and not (_same(lo, z[key + "lo"]) and _same(hi, z[key + "hi"])):
                raise AssertionError(f"recomputed {meth} alpha={a} differs from stored intervals in {f}")
        n_checked += 1
        # (1) de-duplication sensitivity
        kte = ~np.isin(te_id, list(removed))
        kca = ~np.isin(ca_id, list(removed))
        ktr = ~np.isin(tr_id, list(removed))
        ded = intervals_from_npz(z, sg_cal, kte, kca, ktr, alphas, a_p)
        noconf_te = kte & ~np.isin(te_id, list(conflict))
        lat_te = np.array([LATTICE_NAME[s] for s in z["sg_test"]])
        for (meth, a), (lo, hi) in base.items():
            lo2, hi2 = ded[(meth, a)]
            cov0 = (y_te >= lo) & (y_te <= hi)
            cov1 = (y_te[kte] >= lo2) & (y_te[kte] <= hi2)
            # variant 2: v2 data, conflicting records excluded from evaluation only
            sel2 = noconf_te[kte]
            for gname in ["ALL"] + LATTICE_ORDER:
                m0 = np.ones(len(te_id), bool) if gname == "ALL" else lat_te == gname
                m1 = m0[kte]
                m2 = m1 & sel2
                rows_ded.append(dict(repeat=r, target=target, method=meth, alpha=a, group=gname,
                                     cov_phase1=cov0[m0].mean(), cov_v2=cov1[m1].mean(),
                                     cov_v2_noconflict=cov1[m2].mean() if m2.any() else np.nan,
                                     width_med_phase1=float(np.median((hi - lo)[m0])),
                                     width_med_v2=float(np.median((hi2 - lo2)[m1])),
                                     n_phase1=int(m0.sum()), n_v2=int(m1.sum()), n_v2_noconflict=int(m2.sum())))
        # (2) cluster-aware Mondrian on the v2 data
        rng = np.random.default_rng(int(cfg["seed_base"]) + 97 * r + (0 if target == "bandgap" else 1))
        p_te, p_ca = z["pred_split_test"][kte], z["pred_split_cal"][kca]
        s_ca = np.abs(z["y_cal"][kca] - p_ca)
        e_te = np.abs(y_te[kte] - p_te)
        sgc, sgt = sg_cal[kca], z["sg_test"][kte]
        gca, gte = g_ca[kca], g_te[kte]
        for s in sorted(np.unique(sgt)):
            mc, mt = sgc == s, sgt == s
            sc, gc = s_ca[mc], gca[mc]
            et, gt = e_te[mt], gte[mt]
            ug, inv = np.unique(gc, return_inverse=True)
            order = np.argsort(inv, kind="mergesort")
            starts = np.r_[0, np.cumsum(np.bincount(inv))[:-1]]
            sizes = np.bincount(inv)
            for a in alphas:
                q_pool = cp.conformal_quantile(sc, a)
                cov_pool = et <= q_pool
                pw, gw, qs = [], [], []
                for _ in range(B):
                    pick = order[starts + rng.integers(0, sizes)]
                    q = cp.conformal_quantile(sc[pick], a)
                    cov = et <= q
                    pw.append(cov.mean())
                    gw.append(group_weighted(cov, gt))
                    qs.append(q)
                qs = np.asarray(qs)
                rows_clu.append(dict(repeat=r, target=target, lattice=LATTICE_NAME[s], alpha=a,
                                     n_cal=int(mc.sum()), n_cal_groups=int(ug.size), n_test=int(mt.sum()),
                                     n_test_groups=int(np.unique(gt).size),
                                     pooled_cov_point=float(cov_pool.mean()),
                                     pooled_cov_group=group_weighted(cov_pool, gt),
                                     pooled_width=2 * q_pool,
                                     single_sub_cov_point=float(np.mean(pw)),
                                     single_sub_cov_group=float(np.mean(gw)),
                                     single_sub_frac_infinite=float(np.isinf(qs).mean()),
                                     single_sub_width_median=float(np.median(2 * qs))))
        print(f"  {os.path.basename(f)} done", flush=True)
    print(f"  stored Phase 1 intervals reproduced exactly from the saved arrays in {n_checked} files")
    ded = pd.DataFrame(rows_ded)
    clu = pd.DataFrame(rows_clu)
    out = paths.post1
    atomic_write_csv(os.path.join(out, "dedup_sensitivity_by_run.csv"), ded)
    atomic_write_csv(os.path.join(out, "cluster_mondrian_by_run.csv"), clu)
    ds_sum = (ded.groupby(["target", "method", "alpha", "group"], sort=False)
              [["cov_phase1", "cov_v2", "cov_v2_noconflict", "width_med_phase1", "width_med_v2"]].mean().reset_index())
    ds_sum["delta_cov_v2"] = ds_sum.cov_v2 - ds_sum.cov_phase1
    cl_sum = (clu.groupby(["target", "lattice", "alpha"], sort=False)
              [["n_cal_groups", "n_test_groups", "pooled_cov_point", "pooled_cov_group", "single_sub_cov_point",
                "single_sub_cov_group", "single_sub_frac_infinite", "pooled_width", "single_sub_width_median"]]
              .mean().reset_index())
    atomic_write_csv(os.path.join(out, "dedup_sensitivity_summary.csv"), ds_sum)
    atomic_write_csv(os.path.join(out, "cluster_mondrian_summary.csv"), cl_sum)

    # (3) Fd-3m sub-families in Phase 1 (from the Phase 1 report)
    sub_tab = None
    f1 = os.path.join(paths.phase1, "report", "phase1_summary.csv")
    if os.path.exists(f1):
        s1 = pd.read_csv(f1)
        sub_tab = s1[(s1.tag == "hpo_once") & (s1.group.str.startswith("Fd-3m")) &
                     (s1.alpha.round(4) == round(a_p, 4))][["target", "method", "group", "cov_mean", "width_median",
                                                             "mae", "n_test_mean"]]
        atomic_write_csv(os.path.join(out, "fd3m_subfamily_phase1.csv"), sub_tab)

    # digest
    lines = [f"# Phase 1 follow-up checks (exploratory; alpha = {a_p})", "",
             f"Phase 1 files re-analysed: {len(files)}. Stored intervals reproduced exactly: yes.",
             f"v2 de-duplication removes {len(removed)} entries (Phase 1 removed {len(_phase1_removed(paths))}); "
             f"{len(conflict)} conflicting records flagged.", ""]
    for tgt in ("bandgap", "formation"):
        d = ds_sum[(ds_sum.target == tgt) & (ds_sum.alpha.round(4) == round(a_p, 4))]
        lines += [f"## {tgt}: mean coverage, Phase 1 data vs v2 de-duplicated (models unchanged)", "",
                  d.pivot(index="group", columns="method", values="delta_cov_v2").reindex(["ALL"] + LATTICE_ORDER)
                  .round(4).to_markdown(), "",
                  "Coverage on v2 data with conflicting records excluded from evaluation:", "",
                  d.pivot(index="group", columns="method", values="cov_v2_noconflict").reindex(["ALL"] + LATTICE_ORDER)
                  .round(3).to_markdown(), ""]
        c = cl_sum[(cl_sum.target == tgt) & (cl_sum.alpha.round(4) == round(a_p, 4))].set_index("lattice") \
            .reindex(LATTICE_ORDER)
        lines += [f"## {tgt}: Mondrian, pooled vs one-point-per-group calibration (Dunn et al. 2023)", "",
                  c[["n_cal_groups", "pooled_cov_point", "pooled_cov_group", "single_sub_cov_group",
                     "single_sub_frac_infinite", "pooled_width", "single_sub_width_median"]].round(3).to_markdown(), ""]
    if sub_tab is not None:
        lines += ["## Phase 1 by Fd-3m sub-family (cell setting)", "",
                  sub_tab.round(3).to_markdown(index=False), ""]
    atomic_write_text(os.path.join(out, "POST1_SUMMARY.md"), "\n".join(lines))
    print("\n".join(lines))
