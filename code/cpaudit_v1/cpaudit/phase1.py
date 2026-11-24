"""Phase 1: clean in-distribution baseline over repeated composition-grouped splits.

For every repeat r and target, one set of LightGBM fits feeds six interval methods:

  cvplus_random        CV+ with composition-grouped random folds (Regime A; valid design)
  cvplus_lattice       CV+ with one fold per lattice (the submitted manuscript's design), kept
                       as a reference on identical data so the two calibrations can be compared
  split_global         split conformal, one quantile for all lattices
  split_mondrian       split conformal, one quantile per lattice (same model as split_global)
  split_normalized_knn split conformal with |y - y_hat| / sigma_kNN(x) scores (same model)
  cqr                  conformalized quantile regression (alpha_primary only)

All outputs go to <drive_root>/<out_dir>/runs/ and are written atomically; finished
(repeat, target, tag) units are skipped on restart.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
import pandas as pd

from . import conformal as cp
from .config import (atomic_save_npz, atomic_write_csv, atomic_write_json, config_hash,
                     get_paths)
from .data import (LATTICE_NAME, build_dataset, ensure_feature_cache, load_features, make_split,
                   save_split_manifest, split_summary)
from .difficulty import KNNDifficulty
from .metrics import betabinom_prob_below, summarize
from .models import lgb_params, plan_parallel, run_tasks, tune

SUBFAMILY_NAME = {"227a": "Fd-3m (gamma~30)", "227b": "Fd-3m (gamma~120)"}


# ----------------------------------------------------------------------------- guards
def check_preregistration(cfg, paths, allow_deviation: bool):
    f = os.path.join(paths.phase0, "PREREGISTRATION.json")
    if not os.path.exists(f):
        raise RuntimeError("No pre-registration found. Run `python run.py prereg` before Phase 1.")
    with open(f) as fh:
        pre = json.load(fh)
    now = config_hash(cfg)
    if pre["config_hash"] != now:
        msg = (f"Config changed since pre-registration ({pre['config_hash'][:12]} -> {now[:12]}).")
        if not allow_deviation:
            raise RuntimeError(msg + " Re-run with --allow-deviation to proceed; the change will be logged.")
        log = os.path.join(paths.phase0, "DEVIATIONS.json")
        items = json.load(open(log)) if os.path.exists(log) else []
        items.append({"time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                      "from": pre["config_hash"], "to": now, "config": cfg["_config_path"]})
        atomic_write_json(log, items)
        print("  WARNING: " + msg + " Logged in phase0/DEVIATIONS.json")


def get_split(ds, cfg, paths, r):
    split = make_split(ds, cfg, r)
    man = os.path.join(paths.splits, f"split_{r:02d}.npz")
    if os.path.exists(man):  # same seeds must give the same split (guards against library changes)
        with np.load(man) as z:
            ids = ds.df["id"].to_numpy(int)
            for part in ("train", "test", "proper", "cal"):
                if not np.array_equal(np.sort(z[f"{part}_id"]), np.sort(ids[split[part]])):
                    raise RuntimeError(f"split {r} '{part}' differs from the saved manifest {man}; "
                                       "a library version change altered the splitter. Stop and investigate.")
    else:
        save_split_manifest(ds, split, man)
        atomic_write_csv(os.path.join(paths.splits, f"split_{r:02d}_sizes.csv"), split_summary(ds, split))
    return split


# ----------------------------------------------------------------------------- hyperparameters
def hyperparameters(ds, cfg, paths, X_path, r: int, target: str) -> list[tuple[str, dict]]:
    """Returns [(tag, hp), ...] to run for this repeat."""
    mode, checks = cfg["hpo"]["mode"], set(cfg["hpo"].get("check_splits") or [])
    out = []

    def tuned_on(rr):
        sp = get_split(ds, cfg, paths, rr)
        pr = sp["proper"]
        if not os.path.exists(os.path.join(paths.hpo, f"hp_{target}_split{rr:02d}.json")):
            print(f"  [HPO] {target}: tuning on proper-training set of split {rr} (n={len(pr)})", flush=True)
        return tune(X_path, ds.orig_row[pr], ds.y[target][pr], ds.sg[pr], ds.group[pr], cfg,
                    seed=int(cfg["splits"]["seed_base"]) + 5000 + rr,
                    out_json=os.path.join(paths.hpo, f"hp_{target}_split{rr:02d}.json"),
                    log_json=os.path.join(paths.hpo, f"trials_{target}_split{rr:02d}.json"))

    if mode == "once":
        out.append(("hpo_once", tuned_on(0)))
        if r in checks and r != 0:
            out.append(("hpo_own", tuned_on(r)))
    elif mode == "per_split":
        out.append(("hpo_own", tuned_on(r)))
    else:
        raise ValueError(f"unknown hpo.mode {mode}")
    return out


# ----------------------------------------------------------------------------- one unit
def run_unit(ds, split, target, tag, hp, cfg, paths, X_path, X):
    r = split["repeat"]
    stem = os.path.join(paths.runs, f"r{r:02d}_{target}_{tag}")
    if os.path.exists(stem + ".csv") and os.path.exists(stem + ".npz"):
        return "skipped"
    t0 = time.time()
    rows, y = ds.orig_row, ds.y[target]
    tr, te, pr, ca = split["train"], split["test"], split["proper"], split["cal"]
    workers, threads = plan_parallel(cfg)
    seed = int(cfg["splits"]["seed_base"]) + 3000 + r
    a_p = float(cfg["alpha_primary"])
    base = dict(X_path=X_path, max_bin=cfg["lgbm"]["max_bin"], n_rounds=hp["n_rounds"])

    tasks, names = [], []
    tasks.append(dict(base, rows_train=rows[pr], y_train=y[pr], params=lgb_params(hp, cfg, seed, threads),
                      predict={"cal": rows[ca], "test": rows[te]}))
    names.append("split")
    if cfg["cqr"]["enabled"]:
        for nm, q in (("cqr_lo", a_p / 2), ("cqr_hi", 1 - a_p / 2)):
            tasks.append(dict(base, rows_train=rows[pr], y_train=y[pr],
                              params=lgb_params(hp, cfg, seed, threads, quantile=q),
                              predict={"cal": rows[ca], "test": rows[te]}))
            names.append(nm)
    fold_sets = [("cvr", split["cv_random"])]
    if cfg["cvplus"]["run_spacegroup_reference"]:
        fold_sets.append(("cvl", split["cv_lattice"]))
    for fs_name, folds in fold_sets:
        for k in np.unique(folds):
            fit, oof = tr[folds != k], tr[folds == k]
            tasks.append(dict(base, rows_train=rows[fit], y_train=y[fit],
                              params=lgb_params(hp, cfg, seed + 100 + int(k), threads),
                              predict={"oof": rows[oof], "test": rows[te]}))
            names.append(f"{fs_name}_{k}")
    res = dict(zip(names, run_tasks(tasks, workers)))

    # --- difficulty for normalized scores (fitted on proper-training only) ---
    knn = KNNDifficulty(cfg["normalized"]["knn_k"], cfg["normalized"]["pca_components"],
                        cfg["normalized"]["beta"], seed).fit(np.asarray(X[rows[pr]]))
    sig_cal, sig_te = knn.sigma(np.asarray(X[rows[ca]])), knn.sigma(np.asarray(X[rows[te]]))

    y_te, y_ca = y[te], y[ca]
    p_cal, p_te = res["split"]["cal"], res["split"]["test"]
    s_cal = np.abs(y_ca - p_cal)
    intervals, preds = {}, {}

    def cvplus(fs_name, folds):
        ks = np.unique(folds)
        fold_pred = np.stack([res[f"{fs_name}_{k}"]["test"] for k in ks])
        fold_idx = np.searchsorted(ks, folds)
        R = np.empty(len(tr))
        for k in ks:
            R[folds == k] = np.abs(y[tr][folds == k] - res[f"{fs_name}_{k}"]["oof"])
        return fold_pred, fold_idx, R

    cv_parts = {fs: cvplus(fs, folds) for fs, folds in fold_sets}
    for a in cfg["alpha"]:
        a = float(a)
        intervals[("split_global", a)] = cp.split_interval(p_te, s_cal, a)[:2]
        intervals[("split_mondrian", a)] = cp.mondrian_interval(p_te, ds.sg[te], s_cal, ds.sg[ca], a)[:2]
        intervals[("split_normalized_knn", a)] = cp.normalized_interval(p_te, sig_te, y_ca, p_cal, sig_cal, a)[:2]
        for fs, method in (("cvr", "cvplus_random"), ("cvl", "cvplus_lattice")):
            if fs in cv_parts:
                fp, fi, R = cv_parts[fs]
                intervals[(method, a)] = cp.cvplus_interval(fp, fi, R, a)
    if cfg["cqr"]["enabled"]:
        intervals[("cqr", a_p)] = cp.cqr_interval(res["cqr_lo"]["test"], res["cqr_hi"]["test"], y_ca,
                                                  res["cqr_lo"]["cal"], res["cqr_hi"]["cal"], a_p)[:2]
    for m in ("split_global", "split_mondrian", "split_normalized_knn"):
        preds[m] = p_te
    for fs, method in (("cvr", "cvplus_random"), ("cvl", "cvplus_lattice")):
        if fs in cv_parts:
            preds[method] = cv_parts[fs][0].mean(axis=0)   # mean of fold models (no extra refit)
    if cfg["cqr"]["enabled"]:
        preds["cqr"] = p_te

    # --- metrics by lattice, Fd-3m cell setting, and overall ---
    sg_te, sub_te = ds.sg[te], ds.sub[te]
    masks = [("ALL", np.ones(len(te), bool))]
    masks += [(LATTICE_NAME[s], sg_te == s) for s in sorted(np.unique(ds.sg))]
    masks += [(SUBFAMILY_NAME[s], sub_te == s) for s in ("227a", "227b")]
    n_cal_by = {LATTICE_NAME[s]: int((ds.sg[ca] == s).sum()) for s in np.unique(ds.sg)}
    lat_code = {LATTICE_NAME[s]: s for s in np.unique(ds.sg)}
    grp_te, grp_ca = ds.group[te], ds.group[ca]
    out_rows = []
    for (method, a), (lo, hi) in intervals.items():
        for gname, m in masks:
            if m.sum() == 0:
                continue
            row = dict(repeat=r, target=target, tag=tag, method=method, alpha=a, group=gname)
            row.update(summarize(y_te[m], preds[method][m], lo[m], hi[m], a))
            # sampling distribution expected IF the procedure is valid for this group
            if method == "split_global" and gname == "ALL":
                row["p_below_if_valid"] = betabinom_prob_below(len(ca), int(m.sum()), a)
            elif method == "split_mondrian" and gname in n_cal_by:
                row["p_below_if_valid"] = betabinom_prob_below(n_cal_by[gname], int(m.sum()), a)
            else:
                row["p_below_if_valid"] = np.nan
            row["n_cal_group"] = n_cal_by.get(gname, len(ca)) if method.startswith("split") else np.nan
            row["n_test_comp_groups"] = int(np.unique(grp_te[m]).size)
            code = lat_code.get(gname, 227 if gname.startswith("Fd-3m (") else None)
            row["n_cal_comp_groups"] = (int(np.unique(grp_ca[ds.sg[ca] == code]).size) if code is not None
                                        else int(np.unique(grp_ca).size))
            out_rows.append(row)
    df = pd.DataFrame(out_rows)

    ids = ds.df["id"].to_numpy(int)
    arrays = dict(test_id=ids[te], cal_id=ids[ca], train_id=ids[tr], y_test=y_te, y_cal=y_ca,
                  sg_test=sg_te, sub_test=sub_te.astype("U4"), pred_split_test=p_te, pred_split_cal=p_cal,
                  sigma_knn_cal=sig_cal, sigma_knn_test=sig_te,
                  knn_pca_explained=np.array(knn.explained_))
    for fs in cv_parts:
        fp, fi, R = cv_parts[fs]
        arrays[f"{fs}_fold_pred_test"], arrays[f"{fs}_fold_of_train"], arrays[f"{fs}_resid"] = fp, fi, R
    if cfg["cqr"]["enabled"]:
        for nm in ("cqr_lo", "cqr_hi"):
            arrays[f"{nm}_cal"], arrays[f"{nm}_test"] = res[nm]["cal"], res[nm]["test"]
    for (method, a), (lo, hi) in intervals.items():
        arrays[f"{method}__a{a:.2f}__lo"], arrays[f"{method}__a{a:.2f}__hi"] = lo, hi
    atomic_save_npz(stem + ".npz", **arrays)
    atomic_write_csv(stem + ".csv", df)
    atomic_write_json(stem + ".hp.json", hp)
    return f"done in {(time.time() - t0) / 60:.1f} min"


# ----------------------------------------------------------------------------- driver
def run_phase1(cfg, repeats=None, targets=("formation", "bandgap"), allow_deviation=False):
    paths = get_paths(cfg)
    check_preregistration(cfg, paths, allow_deviation)
    X_path = ensure_feature_cache(cfg, paths)
    X = load_features(cfg, paths)
    ds = build_dataset(cfg, paths, X)
    w, t = plan_parallel(cfg)
    repeats = list(range(cfg["splits"]["n_repeats"])) if repeats is None else list(repeats)
    print(f"Phase 1 | n={ds.n} after de-duplication | repeats={repeats[0]}..{repeats[-1]} | "
          f"targets={list(targets)} | {w} parallel fits x {t} threads", flush=True)
    t_start, n_done = time.time(), 0
    units = [(r, tg) for r in repeats for tg in targets]
    for i, (r, target) in enumerate(units, 1):
        split = get_split(ds, cfg, paths, r)
        for tag, hp in hyperparameters(ds, cfg, paths, X_path, r, target):
            status = run_unit(ds, split, target, tag, hp, cfg, paths, X_path, X)
            if status != "skipped":
                n_done += 1
            el = (time.time() - t_start) / 60
            print(f"  [{i}/{len(units)}] repeat {r:02d} {target:<9} {tag:<8} {status}   (elapsed {el:.0f} min)",
                  flush=True)
    print("Phase 1 finished. Next: python run.py report")
