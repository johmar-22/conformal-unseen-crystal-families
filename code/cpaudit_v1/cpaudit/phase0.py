"""Phase 0: freeze the current state, verify the inputs, and pre-register the Phase 1 analysis.

Task 0.2  freeze_manifest()  SHA-256 of every existing file in the Drive folder (except the
                             2,400 geometry folders, which are counted) + software environment.
          integrity()        reproduces the v5 split exactly, checks SOAP row alignment against
                             data_prep_v5.npz, re-derives the duplicate groups and compares them
                             with duplicate_entries_v7.csv, and reports group / sub-family counts.
Task 0.3  preregister()      writes an immutable PREREGISTRATION.md/.json holding the hypotheses,
                             endpoints, decision rules and the hashes of config and code.
Task 0.1  (outside expert read) is not code; see README.
"""
from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd

from .config import (atomic_write_csv, atomic_write_json, atomic_write_text, code_hash, config_hash,
                     get_paths, sha256_file)
from .data import (LATTICE_NAME, build_dataset, composition_groups, ensure_feature_cache,
                   find_duplicates, load_features, load_train_csv)
from .models import n_cpu

_SKIP_DIRS = {"train", "__pycache__"}


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


# ----------------------------------------------------------------------------- environment
def environment() -> dict:
    import importlib
    pk = {}
    for m in ("numpy", "scipy", "pandas", "sklearn", "lightgbm", "optuna", "crepes", "mapie", "yaml", "joblib"):
        try:
            pk[m] = importlib.import_module(m).__version__
        except Exception:
            pk[m] = None
    gpu = None
    try:
        gpu = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=10).stdout.strip() or None
    except Exception:
        pass
    mem = None
    try:
        with open("/proc/meminfo") as f:
            mem = round(int(f.readline().split()[1]) / 1e6, 1)
    except Exception:
        pass
    return dict(time_utc=_now(), python=sys.version.split()[0], platform=platform.platform(),
                cpus=n_cpu(), ram_gb=mem, gpu=gpu, packages=pk, lightgbm_gpu_ok=_lgbm_gpu_probe())


def _lgbm_gpu_probe() -> bool:
    try:
        import lightgbm as lgb
        X = np.random.default_rng(0).standard_normal((200, 5))
        lgb.train({"device_type": "gpu", "verbosity": -1, "objective": "regression"},
                  lgb.Dataset(X, X[:, 0]), num_boost_round=2)
        return True
    except Exception:
        return False


# ----------------------------------------------------------------------------- freeze
def freeze_manifest(cfg, hash_large: bool = True) -> str:
    paths = get_paths(cfg)
    out = os.path.join(paths.phase0, "freeze_manifest_v7.json")
    if os.path.exists(out):
        print(f"  freeze manifest already exists (kept unchanged): {out}")
        return out
    root, own = paths.drive_root, os.path.basename(paths.out)
    files, n_geom = [], 0
    for dirpath, dirnames, filenames in os.walk(root):
        rel = os.path.relpath(dirpath, root)
        top = rel.split(os.sep)[0]
        if top == own:
            dirnames[:] = []
            continue
        if top in _SKIP_DIRS:
            if top == "train":
                n_geom += len(filenames)
            continue
        for fn in sorted(filenames):
            p = os.path.join(dirpath, fn)
            size = os.path.getsize(p)
            h = sha256_file(p) if (hash_large or size < 50e6) else None
            files.append(dict(path=os.path.relpath(p, root), bytes=size, sha256=h))
            print(f"    hashed {files[-1]['path']}", flush=True)
    atomic_write_json(out, dict(created_utc=_now(), drive_root=root, n_files=len(files),
                                geometry_files_counted=n_geom, files=files, environment=environment()))
    print(f"  wrote {out}  ({len(files)} files hashed, {n_geom} geometry files counted)")
    return out


# ----------------------------------------------------------------------------- integrity
def integrity(cfg) -> dict:
    paths = get_paths(cfg)
    rep = {"time_utc": _now()}
    df, csv_hash = load_train_csv(cfg, paths)
    rep["train_csv_sha256"] = csv_hash
    X = load_features(cfg, paths)
    rep["soap_shape"] = list(X.shape)
    ensure_feature_cache(cfg, paths)

    # (1) the v5 split is reproduced exactly from train.csv order
    from sklearn.model_selection import StratifiedShuffleSplit
    sg = df["spacegroup"].to_numpy(int)
    tr, te = next(StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42).split(np.zeros(len(df)), sg))
    prep = paths.drive_file(cfg["paths"]["data_prep_v5"])
    if os.path.exists(prep):
        with np.load(prep, allow_pickle=False) as z:
            same_idx = bool(np.array_equal(z["train_idx"], tr) and np.array_equal(z["test_idx"], te))
            y_ok = bool(np.allclose(z["y_test_band"], df["bandgap_energy_ev"].to_numpy()[te], atol=1e-4)
                        and np.allclose(z["y_test_form"], df["formation_energy_ev_natom"].to_numpy()[te], atol=1e-4))
            cols = np.arange(0, X.shape[1], max(1, X.shape[1] // 300))
            Xtr = np.asarray(X[tr][:, cols], dtype=np.float64)
            mu, sd = Xtr.mean(0), Xtr.std(0)
            sd[sd == 0] = 1.0
            expect = (np.asarray(X[te[:50]][:, cols], np.float64) - mu) / sd
            got = z["X_test_scaled"][:50, cols].astype(np.float64)
            x_ok = bool(np.allclose(expect, got, atol=2e-3))
        rep["v5_split_reproduced"] = same_idx
        rep["v5_targets_aligned"] = y_ok
        rep["soap_rows_aligned_with_v5"] = x_ok
    else:
        rep["v5_split_reproduced"] = "data_prep_v5.npz not found (skipped)"

    # (2) duplicates re-derived and compared with the v7 table
    dups = find_duplicates(df, X)
    ids = [[int(df.loc[r, "id"]) for r in g] for g in dups]
    rep["duplicate_groups"] = len(ids)
    rep["duplicate_entries"] = int(sum(len(g) for g in ids))
    rep["entries_removed"] = int(sum(len(g) - 1 for g in ids))
    from .data import MAX_DUP_SOAP_RELDIFF
    rep["max_soap_rel_diff_within_duplicates"] = MAX_DUP_SOAP_RELDIFF["value"]
    v7 = paths.drive_file(cfg["paths"]["duplicates_v7"])
    if os.path.exists(v7):
        d7 = pd.read_csv(v7)
        set7 = {tuple(sorted(g["id"].tolist())) for _, g in d7.groupby("group")}
        setn = {tuple(sorted(g)) for g in ids}
        rep["duplicates_match_v7"] = bool(set7 == setn)
        rep["duplicates_only_in_v7"] = sorted(map(list, set7 - setn))
        rep["duplicates_only_new"] = sorted(map(list, setn - set7))
    atomic_write_csv(os.path.join(paths.phase0, "duplicate_groups.csv"),
                     pd.DataFrame([dict(group=i, ids=" ".join(map(str, g)), kept=g[0]) for i, g in enumerate(ids)]))

    # (3) composition groups (identical labels to the v7 leakage table) and sub-families
    ds = build_dataset(cfg, paths, X)
    rep["n_after_dedup"] = ds.n
    vc = pd.Series(ds.group).value_counts()
    rep["composition_groups"] = int(vc.size)
    rep["composition_group_size_max"] = int(vc.max())
    rep["composition_singletons"] = int((vc == 1).sum())
    g_v5_test = set(composition_groups(df.iloc[te]))
    rep["v5_test_composition_groups"] = len(g_v5_test)   # manuscript reports 231
    per = []
    for s in sorted(np.unique(ds.sg)):
        m = ds.sg == s
        row = dict(lattice=LATTICE_NAME[s], sg=int(s), n=int(m.sum()), composition_groups=int(np.unique(ds.group[m]).size))
        if s == 227:
            row["n_227a_gamma_lt_split"] = int((ds.sub == "227a").sum())
            row["n_227b_gamma_ge_split"] = int((ds.sub == "227b").sum())
        per.append(row)
    rep["per_lattice"] = per
    atomic_write_json(os.path.join(paths.phase0, "integrity_report.json"), rep)
    return rep


# ----------------------------------------------------------------------------- pre-registration
PREREG_TEXT = """# Pre-registration: Phase 1 (in-distribution baseline)

Registered (UTC): {time}
Config hash: {chash}
Code hash:   {khash}
Config file: {cpath}

This file is written once and never overwritten. Any later change to the scientific
settings is refused by the runner unless explicitly allowed, and is then logged in
DEVIATIONS.json next to this file.

## Design
* Data: NOMAD 2018 public training set (2,400 entries), exact duplicates removed
  (identical tabular record and identical SOAP vector; lowest id kept).
* Unit of exchangeability: composition group (lattice, Al and Ga fraction rounded to 0.1).
  Outer test set, calibration set and all CV folds are group-disjoint and stratified by lattice.
* {R} independent outer splits (test = 20 %), calibration = 25 % of the training pool.
* Model: LightGBM (L2 objective), hyperparameters tuned with Optuna TPE ({ntrials} trials, MAE,
  {inner}-fold composition-grouped CV) on the proper-training set of split 0 (mode = {hmode});
  splits {checks} are also re-tuned on their own proper-training set to measure sensitivity.
* Interval methods: cvplus_random (K={K}), cvplus_lattice (one fold per lattice, the submitted
  design, reference only), split_global, split_mondrian (per lattice), split_normalized_knn
  (k={knnk}, {pca} PCA components), cqr (alpha = {ap} only).
* alpha in {alphas}; primary alpha = {ap}.

## Hypotheses (Phase 1 part of H1)
H1a  Under a valid in-distribution calibration (cvplus_random and split_global), per-lattice
     coverage is heterogeneous and Fd-3m is under-covered.
H1b  Mondrian (per-lattice) calibration restores per-lattice coverage to nominal, at a width cost.
H1c  (descriptive) The lattice-fold CV+ design of the submitted manuscript over-covers every
     lattice relative to cvplus_random, on identical data.

## Primary endpoint
Mean over the {R} splits of the band-gap coverage of Fd-3m test points under cvplus_random at
alpha = {ap}. Same endpoint for formation energy is co-primary.

## Decision rules (fixed before any Phase 1 result is seen)
* H1a supported   if the primary endpoint < 0.85.
* H1a rejected    if the primary endpoint >= 0.88  -> decision gate G1: the in-distribution
                  heterogeneity claim is dropped and reported as a negative result.
* 0.85-0.88       inconclusive; reported as such, no claim made.
* H1b supported   if, under split_mondrian at alpha = {ap}, the mean coverage of every lattice
                  is >= 0.88 and the observed fraction of splits below 0.90 is not larger than
                  the Beta-binomial expectation for a valid procedure + 0.15.
* Every endpoint, method, alpha and lattice is reported whichever way it falls.

## Sampling-noise reference
For split-conformal methods the expected fraction of splits with coverage below 1 - alpha under
exact validity is computed per split from Beta-binomial(n_test, k, n_cal + 1 - k) with
k = ceil((1 - alpha)(n_cal + 1)) and stored as p_below_if_valid. Because points within a
composition group are correlated, the true split-to-split spread is larger than this
independent-point reference; it is therefore a lower bound on the noise, and the number of
test and calibration composition groups is reported alongside every coverage value.
"""


def preregister(cfg) -> str:
    paths = get_paths(cfg)
    js = os.path.join(paths.phase0, "PREREGISTRATION.json")
    md = os.path.join(paths.phase0, "PREREGISTRATION.md")
    if os.path.exists(js):
        with open(js) as f:
            pre = json.load(f)
        same = pre["config_hash"] == config_hash(cfg)
        print(f"  pre-registration already exists ({pre['time_utc']}); config unchanged: {same}")
        return md
    if any(f.endswith(".csv") for f in os.listdir(paths.runs)):
        raise RuntimeError("Phase 1 results already exist; a pre-registration written now would not be one.")
    info = dict(time_utc=_now(), config_hash=config_hash(cfg), code_hash=code_hash(),
                config_path=cfg["_config_path"], config={k: v for k, v in cfg.items() if not k.startswith("_")})
    h = cfg["hpo"]
    text = PREREG_TEXT.format(time=info["time_utc"], chash=info["config_hash"], khash=info["code_hash"],
                              cpath=info["config_path"], R=cfg["splits"]["n_repeats"], ntrials=h["n_trials"],
                              inner=h["inner_folds"], hmode=h["mode"], checks=h.get("check_splits"),
                              K=cfg["cvplus"]["k_random"], knnk=cfg["normalized"]["knn_k"],
                              pca=cfg["normalized"]["pca_components"], ap=cfg["alpha_primary"],
                              alphas=cfg["alpha"])
    atomic_write_text(md, text)
    atomic_write_json(js, info)
    print(f"  wrote {md}")
    return md


def run_phase0(cfg, hash_large=True):
    print("Phase 0 | Task 0.2 freeze manifest")
    freeze_manifest(cfg, hash_large=hash_large)
    print("Phase 0 | integrity checks")
    rep = integrity(cfg)
    keys = ["soap_shape", "v5_split_reproduced", "v5_targets_aligned", "soap_rows_aligned_with_v5",
            "duplicate_groups", "duplicate_entries", "entries_removed", "duplicates_match_v7",
            "n_after_dedup", "composition_groups", "v5_test_composition_groups"]
    for k in keys:
        print(f"  {k:<28} {rep.get(k)}")
    for row in rep["per_lattice"]:
        print("  ", row)
    bad = [k for k in ("v5_split_reproduced", "v5_targets_aligned", "soap_rows_aligned_with_v5",
                       "duplicates_match_v7") if rep.get(k) is False]
    if bad:
        print(f"  ATTENTION: failed checks {bad}. Do not start Phase 1 until these are understood.")
    else:
        print("  all integrity checks passed")
    env = environment()
    atomic_write_json(os.path.join(get_paths(cfg).phase0, "environment.json"), env)
    print(f"  environment: {env['cpus']} CPUs, {env['ram_gb']} GB RAM, GPU={env['gpu']}, "
          f"LightGBM GPU build usable={env['lightgbm_gpu_ok']}")
    return rep
