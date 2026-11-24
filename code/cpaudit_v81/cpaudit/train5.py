"""Phase 5, optional training arm: should new labels of a family go into the model or into the
calibration set? (Plan, Phase 5 'optional variant'; descriptive, no decision rule.)

Per family and draw: a random 40 % of the family's composition groups are the test set (new
compositions). From the remaining groups, n_cal structures are drawn for calibration and, nested,
m more for training. Arms (evaluated on the same test structures):
  train      LightGBM refit on the Phase 2 training lattices + the m structures, calibrated on n_cal
  calibrate  the fixed Phase 2 model, calibrated on the m + n_cal structures
  fixed      the fixed Phase 2 model, calibrated on n_cal only            (m = 0)
  baseline   a refit on the Phase 2 training lattices only, calibrated on n_cal (m = 0)
Methods: split and offset_full (cpaudit/fewshot.py). Hyperparameters: the Phase 2 LOSO set of the
family (the leave-two-out units reused them too). Each draw is saved as it finishes (resume).
"""
from __future__ import annotations

import json
import math
import os
import time

import numpy as np
import pandas as pd

from .config import atomic_save_npz, atomic_write_csv, get_paths
from .data import LATTICE_NAME, build_dataset, ensure_feature_cache
from .fewshot import SET_FUNCS, hull_width, in_set
from .models import lgb_params, plan_parallel, run_tasks
from .phase5 import _seed, parse_config
from .setup5 import check_preregistration

TARGETS = ("bandgap", "formation")


def _evaluate(z_cal, z_test, test_groups_idx, alpha):
    out = []
    for meth in ("split", "offset_full"):
        S = SET_FUNCS[meth](z_cal, alpha)
        c = in_set(z_test, S)
        cg = float(np.mean([c[ix].mean() for ix in test_groups_idx]))
        out.append(dict(method=meth, cov=float(c.mean()), cov_group=cg, width=hull_width(S),
                        mae=float(np.mean(np.abs(z_test)))))
    return out


def run_train5(cfg, allow_deviation=False):
    paths = get_paths(cfg)
    check_preregistration(cfg, paths, allow_deviation)
    t = cfg["train"]
    alpha = float(cfg["alpha_primary"])
    outdir = os.path.join(paths.out, "train")
    os.makedirs(outdir, exist_ok=True)
    ds = build_dataset(cfg, paths)
    X_path = ensure_feature_cache(cfg, paths)
    row_of = {int(i): r for r, i in zip(ds.orig_row, ds.ids)}
    workers, threads = plan_parallel(cfg)
    m_grid = sorted(int(m) for m in t["m"])
    t0 = time.time()
    for c in t["configs"]:
        design, g = parse_config(c)
        for target in TARGETS:
            z = np.load(os.path.join(paths.phase2_loso, f"{target}_{design}{g}.npz"), allow_pickle=False)
            hp = json.load(open(os.path.join(paths.phase2_loso, f"{target}_loso{g}.hp.json")))
            tr_ids, te_ids = z["train_id"].astype(int), z["test_id"].astype(int)
            mu_fixed = z["mu_te_l"].mean(0) if design == "loso" else z["mu_te"].mean(0)
            y_fam = z["y_test"].astype(float)
            y_tr = z["y_train"].astype(float) if "y_train" in z.files else ds.y[target][np.searchsorted(ds.ids, tr_ids)]
            if not np.allclose(y_fam, ds.y[target][np.searchsorted(ds.ids, te_ids)]):
                raise AssertionError(f"{target} {design}{g}: Phase 2 targets differ from the dataset")
            rows_tr = np.array([row_of[i] for i in tr_ids])
            rows_fam = np.array([row_of[i] for i in te_ids])
            fam_groups = ds.group[np.searchsorted(ds.ids, te_ids)]
            ug, ginv = np.unique(fam_groups, return_inverse=True)
            base_kw = dict(X_path=X_path, max_bin=cfg["lgbm"]["max_bin"], n_rounds=hp["n_rounds"])
            label = f"{target} {LATTICE_NAME[g]} ({design})"
            bfile = os.path.join(outdir, f"baseline_{target}_{design}{g}.npz")
            if os.path.exists(bfile):
                mu_base = np.load(bfile)["pred"]
            else:
                r = run_tasks([dict(base_kw, rows_train=rows_tr, y_train=y_tr,
                                    params=lgb_params(hp, cfg, _seed(cfg["seed_base"], "base", c, target) % 100000,
                                                      threads), predict={"fam": rows_fam})], 1)[0]
                mu_base = r["fam"]
                atomic_save_npz(bfile, pred=mu_base)
            for d in range(int(t["draws"])):
                fn = os.path.join(outdir, f"draw_{target}_{design}{g}_{d:02d}.csv")
                if os.path.exists(fn):
                    continue
                t1 = time.time()
                rng = np.random.default_rng(_seed(cfg["seed_base"], "train", c, target, d))
                order = rng.permutation(len(ug))
                n_test_g = int(math.ceil(t["test_group_fraction"] * len(ug)))
                test_g = order[:n_test_g]
                test = np.where(np.isin(ginv, test_g))[0]
                pool = rng.permutation(np.where(~np.isin(ginv, test_g))[0])
                ncal = int(t["n_cal"])
                if len(pool) < ncal + m_grid[-1]:
                    raise RuntimeError(f"{label}: label pool {len(pool)} < n_cal + max m")
                cal, add = pool[:ncal], pool[ncal:]
                groups_idx = [np.where(ginv[test] == gg)[0] for gg in test_g]
                tasks = []
                for m in m_grid:
                    a = add[:m]
                    tasks.append(dict(base_kw, rows_train=np.concatenate([rows_tr, rows_fam[a]]),
                                      y_train=np.concatenate([y_tr, y_fam[a]]),
                                      params=lgb_params(hp, cfg, _seed(cfg["seed_base"], "fit", c, target, d, m) % 100000,
                                                        threads), predict={"fam": rows_fam}))
                res = run_tasks(tasks, workers)
                rows = []

                def add_rows(arm, m, mu, cal_idx):
                    zz = y_fam - mu
                    for rr in _evaluate(zz[cal_idx], zz[test], groups_idx, alpha):
                        rows.append(dict(design=design, heldout=LATTICE_NAME[g], target=target, draw=d, m=m,
                                         arm=arm, n_test=len(test), **rr))
                add_rows("fixed", 0, mu_fixed, cal)
                add_rows("baseline", 0, mu_base, cal)
                for m, r in zip(m_grid, res):
                    add_rows("train", m, r["fam"], cal)
                    add_rows("calibrate", m, mu_fixed, np.concatenate([cal, add[:m]]))
                atomic_write_csv(fn, pd.DataFrame(rows))
                print(f"  {label} draw {d + 1}/{t['draws']} ({time.time() - t1:.0f} s; "
                      f"total {(time.time() - t0) / 60:.1f} min)", flush=True)
    print("train5 finished. Run report5 to include the training arm in the summary.")
