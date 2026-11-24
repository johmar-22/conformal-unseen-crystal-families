"""Setup for Phase 2: environment, data integrity under the v2 de-duplication rule, and the
immutable Phase 2/3 pre-registration (written before any Phase 2 or Phase 3 result exists)."""
from __future__ import annotations

import json
import os
import platform
import sys
import time

import numpy as np
import pandas as pd

from .config import (atomic_write_csv, atomic_write_json, atomic_write_text, code_hash, config_hash,
                     get_paths)
from .data import LATTICE_NAME, build_dataset, load_features, load_train_csv
from .models import n_cpu


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def environment() -> dict:
    import importlib
    pk = {}
    for m in ("numpy", "scipy", "pandas", "sklearn", "lightgbm", "optuna", "pymatgen", "yaml", "joblib"):
        try:
            mod = importlib.import_module(m if m != "pymatgen" else "pymatgen.core")
            pk[m] = getattr(mod, "__version__", None)
        except Exception:
            pk[m] = None
    if pk.get("pymatgen") is None:
        try:
            from importlib.metadata import version
            pk["pymatgen"] = version("pymatgen")
        except Exception:
            pass
    mem = None
    try:
        with open("/proc/meminfo") as f:
            mem = round(int(f.readline().split()[1]) / 1e6, 1)
    except Exception:
        pass
    return dict(time_utc=_now(), python=sys.version.split()[0], platform=platform.platform(),
                cpus=n_cpu(), ram_gb=mem, packages=pk)


# ----------------------------------------------------------------------------- integrity
def _phase1_removed(paths) -> list[int]:
    f = os.path.join(paths.phase1, "phase0", "duplicate_groups.csv")
    if not os.path.exists(f):
        return []
    d = pd.read_csv(f)
    out = []
    for _, r in d.iterrows():
        ids = [int(v) for v in str(r["ids"]).split()]
        out += [i for i in ids if i != int(r["kept"])]
    return sorted(out)


def integrity(cfg) -> dict:
    paths = get_paths(cfg)
    df, h = load_train_csv(cfg, paths)
    X = load_features(cfg, paths)
    ds = build_dataset(cfg, paths)              # computes and caches phase0/dedup_v2.json
    dd = ds.dedup
    by_id = df.set_index("id", drop=False)
    rep = dict(time_utc=_now(), train_csv_sha256=h, soap_shape=list(X.shape), rule=dd["rule"],
               matcher=dd["matcher"], target_tol=dd["target_tol"])

    # feature scale for reporting SOAP differences inside clusters
    sd = np.asarray(X[::4], np.float64).std(axis=0)
    sd[sd == 0] = 1.0
    rows = []
    status_of = {}
    for g in dd["duplicate_groups"]:
        for i in g:
            status_of[i] = "duplicate"
    for c in dd["conflicting_clusters"]:
        for s in c:
            for i in s:
                status_of[i] = "conflicting"
    for c in dd["identical_input_clusters"]:
        sub = by_id.loc[c]
        Xc = np.asarray(X[np.asarray(c) - 1], np.float64)
        dz = np.abs(Xc - Xc[0]) / sd
        rows.append(dict(ids=" ".join(map(str, c)), sg=int(sub.spacegroup.iloc[0]),
                         lattice=LATTICE_NAME[int(sub.spacegroup.iloc[0])], n=len(c),
                         d_formation=float(np.ptp(sub.formation_energy_ev_natom)),
                         d_bandgap=float(np.ptp(sub.bandgap_energy_ev)),
                         status="conflicting" if any(status_of.get(i) == "conflicting" for i in c) else "duplicate",
                         soap_max_abs_diff_in_sd=float(dz.max()),
                         soap_rel_l2_diff=float(np.max(np.linalg.norm(Xc - Xc[0], axis=1)) /
                                                max(np.linalg.norm(Xc[0]), 1e-12))))
    table = pd.DataFrame(rows).sort_values(["status", "d_bandgap"], ascending=[True, False])
    atomic_write_csv(os.path.join(paths.phase0, "dedup_table_v2.csv"), table)

    rep["identical_input_clusters"] = len(dd["identical_input_clusters"])
    rep["identical_input_entries"] = int(sum(len(c) for c in dd["identical_input_clusters"]))
    rep["duplicate_groups"] = len(dd["duplicate_groups"])
    rep["entries_removed"] = len(dd["removed_ids"])
    rep["conflicting_clusters"] = len(dd["conflicting_clusters"])
    rep["conflicting_entries"] = len(dd["conflicting_ids"])
    rep["conflicts"] = table[table.status == "conflicting"][["ids", "lattice", "d_formation", "d_bandgap",
                                                             "soap_max_abs_diff_in_sd"]].to_dict("records")
    rep["max_soap_diff_in_sd_within_clusters"] = float(table.soap_max_abs_diff_in_sd.max()) if len(table) else 0.0
    chk = dd["independent_check"]
    rep["independent_check_pairs"] = len(chk)
    rep["independent_check_same_setting_matches"] = int(sum(c["same_setting_match"] for c in chk))
    rep["n_after_dedup"] = ds.n
    removed = set(dd["removed_ids"])
    p1 = set(_phase1_removed(paths))
    rep["phase1_removed"] = len(p1)
    rep["phase1_removed_subset_of_v2"] = p1 <= removed
    extra = sorted(removed - p1)
    rep["extra_removed_vs_phase1"] = len(extra)
    rep["extra_removed_by_lattice"] = {LATTICE_NAME[s]: int(c) for s, c in
                                       pd.Series(by_id.loc[extra, "spacegroup"]).value_counts().items()} if extra else {}
    v7 = paths.drive_file(cfg["paths"]["duplicates_v7"])
    if os.path.exists(v7):
        d7 = pd.read_csv(v7)
        rem7 = set()
        for _, g in d7.groupby("group"):
            ids = sorted(g["id"].astype(int))
            rem7 |= set(ids[1:])
        rep["v7_removed"] = len(rem7)
        rep["v7_removed_subset_of_v2"] = rem7 <= removed
    vc = pd.Series(ds.group).value_counts()
    rep["composition_groups"] = int(vc.size)
    per = []
    for s in sorted(np.unique(ds.sg)):
        m = ds.sg == s
        row = dict(lattice=LATTICE_NAME[s], sg=int(s), n=int(m.sum()), composition_groups=int(np.unique(ds.group[m]).size),
                   conflicting=int((ds.conflict & m).sum()))
        if s == 227:
            row["n_227a"] = int((ds.sub == "227a").sum())
            row["n_227b"] = int((ds.sub == "227b").sum())
        per.append(row)
    rep["per_lattice"] = per
    atomic_write_json(os.path.join(paths.phase0, "integrity_report_v2.json"), rep)
    return rep


def run_integrity(cfg):
    paths = get_paths(cfg)
    env = environment()
    atomic_write_json(os.path.join(paths.phase0, "environment.json"), env)
    print(f"  environment: {env['cpus']} CPUs, {env['ram_gb']} GB RAM, packages {env['packages']}")
    rep = integrity(cfg)
    for k in ("identical_input_clusters", "identical_input_entries", "duplicate_groups", "entries_removed",
              "conflicting_clusters", "conflicting_entries", "n_after_dedup", "phase1_removed",
              "phase1_removed_subset_of_v2", "extra_removed_vs_phase1", "extra_removed_by_lattice",
              "v7_removed", "v7_removed_subset_of_v2", "independent_check_pairs",
              "independent_check_same_setting_matches", "max_soap_diff_in_sd_within_clusters",
              "composition_groups"):
        print(f"  {k:<40} {rep.get(k)}")
    for c in rep["conflicts"]:
        print(f"  conflicting record: ids {c['ids']:<12} {c['lattice']:<8} dEf={c['d_formation']:.4f} "
              f"dEg={c['d_bandgap']:.4f}  SOAP diff (sd units)={c['soap_max_abs_diff_in_sd']:.2e}")
    for row in rep["per_lattice"]:
        print("  ", row)
    return rep


# ----------------------------------------------------------------------------- pre-registration
PREREG_TEXT = """# Pre-registration: Phase 2 (LOSO mechanism) and Phase 3 (overlap diagnostics)

Registered (UTC): {time}
Config hash: {chash}
Code hash:   {khash}

Written once, before any Phase 2 or Phase 3 result exists; never overwritten. A later change
to the scientific settings is refused unless explicitly allowed and is then logged in
DEVIATIONS.json next to this file. Phase 1 results were known when this was written; the
Phase 1 follow-up checks (step post1) are therefore labelled exploratory, not confirmatory.

## Data
* NOMAD 2018 training set. v2 de-duplication: same input structure (pymatgen StructureMatcher,
  ltol = stol = {ltol}, angle_tol = {atol} deg, on the Vegard geometry files) AND both targets within
  {ttol} eV(/atom) -> one record (lowest id kept). n = {n} entries after de-duplication.
  The tolerance was chosen after inspecting the distribution of target differences inside
  identical-input clusters (largest within-tolerance difference 0.0015 eV; every conflicting cluster
  differs by >= 0.0125 eV(/atom) in at least one target),
  before any Phase 2 model was fitted.
* Same-input entries whose targets differ by more ({nconf} entries in {nconfc} clusters) are
  kept and flagged; every Phase 2 coverage is also reported with them excluded from evaluation.

## Design (nested leave-one-lattice-out, LOSO)
* For each held-out lattice g (order {order}) and target: hyperparameters tuned with Optuna TPE
  ({ntrials} trials, {nstart} random start-up, MAE, early stopping {es}, max {mr} rounds) using
  leave-one-lattice-out inner folds over the five training lattices only (g is never seen).
* CV+ with the same five lattice folds (cvplus_lattice). Nested scores with the same folds
  (Gupta et al. 2022): normalized by a kNN distance sigma fitted on each fold's training lattices
  (ncvplus_knn; k = {knnk}, {pca} PCA components) and CQR (cqr_cvplus, alpha = {ap}).
  Contrast: CV+ with {krand} composition-grouped random folds on the same training lattices
  (cvplus_random_in): in-distribution residuals only.
* Calibration-only ablation (no refit): the same fold models and residuals with the Fd-3m fold
  removed from the calibration terms (suffix _noFd3m_cal), for every g other than Fd-3m.
* Leave-two-out refit (Task 2.2): for every X other than Fd-3m, train and calibrate on the four
  lattices other than X and Fd-3m (4 lattice folds, hyperparameters of LOSO-X), test on X.
* Group-level conformal (Task 2.4, illustrative): family score = {tau}-quantile of a lattice's
  out-of-fold residuals; threshold = ceil((K+1)(1-alpha))-th smallest of K = 5 family scores
  (Dunn et al. 2023). alpha = 0.1 gives an unbounded interval by construction; alpha = 0.2 the
  largest of the five.
* alpha in {alphas}; primary alpha = {ap}.

## Phase 3 overlap gate (fixed now; diagnostics are computed before LOSO results)
* Classifiers: held-out lattice vs the other five, logistic regression (C = {C}) and kNN (k = {kk})
  on {opca} PCA components of the standardized SOAP vector, {ocv}-fold composition-grouped
  out-of-fold probabilities. AUROC of each; density-ratio weights w = p/(1-p) * n_rest/n_g
  for the five calibration lattices; Kish ESS = (sum w)^2 / sum w^2; test-point weight share
  p_test = w(x) / (sum_i w_i + w(x)) (weighted CP is unbounded when p_test > alpha).
* Weighted CV+ is run for holdout g only if max(AUROC) <= {amax} AND ESS >= {essmin} AND the
  median p_test < {ap}. Otherwise weighted CP is reported as inapplicable for that unseen
  lattice, with these numbers.

## Decision rules (alpha = {ap}, each target separately)
* P2.1 (Task 2.1 prediction): for every holdout g other than Fd-3m, Fd-3m's share of the
  calibration residuals at or above the pooled conformal quantile, divided by its share of
  calibration points, is >= 2. Supported if this holds in >= 4 of the 5 holdouts.
* H2 (primary contrast = calibration-only ablation, cvplus_lattice): for X other than Fd-3m let
  dcov = coverage(ablation) - coverage(LOSO) and r = median width(ablation) / median width(LOSO).
  H2 supported if r <= 0.80 and dcov <= -0.05 in >= 4 of 5 X.
  H2 rejected (gate G2) if r >= 0.90 and |dcov| <= 0.03 in >= 4 of 5 X. Otherwise partial.
  The same quantities for the leave-two-out refit are reported as the plan's Task 2.2 result.
* G3: positive if, under LOSO-Fd-3m, ncvplus_knn or cqr_cvplus covers >= 0.85 of Fd-3m.
* Task 2.3 (descriptive): a sub-family 'carries the failure' if its LOSO-Fd-3m cvplus_lattice
  coverage is < 0.75 while the other's is >= 0.85.
* Every lattice, method and alpha is reported whichever way it falls.
"""


def prereg_path(paths):
    return (os.path.join(paths.phase0, "PREREGISTRATION_phase2.json"),
            os.path.join(paths.phase0, "PREREGISTRATION_phase2.md"))


def preregister(cfg) -> str:
    paths = get_paths(cfg)
    js, md = prereg_path(paths)
    if os.path.exists(js):
        with open(js) as f:
            pre = json.load(f)
        print(f"  pre-registration already exists ({pre['time_utc']}); config unchanged: "
              f"{pre['config_hash'] == config_hash(cfg)}")
        return md
    if os.listdir(paths.loso) or os.path.exists(os.path.join(paths.overlap, "overlap_diagnostics.csv")):
        raise RuntimeError("Phase 2/3 results already exist; a pre-registration written now would not be one.")
    rep_f = os.path.join(paths.phase0, "integrity_report_v2.json")
    if not os.path.exists(rep_f):
        raise RuntimeError("Run the integrity step first (setup does this).")
    rep = json.load(open(rep_f))
    info = dict(time_utc=_now(), config_hash=config_hash(cfg), code_hash=code_hash(),
                config_path=cfg["_config_path"], config={k: v for k, v in cfg.items() if not k.startswith("_")},
                n_after_dedup=rep["n_after_dedup"])
    L, d, o = cfg["loso"], cfg["data"]["dedup"], cfg["overlap"]
    text = PREREG_TEXT.format(
        time=info["time_utc"], chash=info["config_hash"], khash=info["code_hash"],
        ltol=d["matcher"]["ltol"], atol=d["matcher"]["angle_tol"], ttol=d["target_tol"], n=rep["n_after_dedup"],
        nconf=rep["conflicting_entries"], nconfc=rep["conflicting_clusters"],
        order=[LATTICE_NAME[s] for s in L["order"]], ntrials=L["hpo"]["n_trials"],
        nstart=L["hpo"]["n_startup_trials"], es=L["hpo"]["early_stopping_rounds"], mr=L["hpo"]["max_rounds"],
        knnk=L["normalized"]["knn_k"], pca=L["normalized"]["pca_components"], ap=cfg["alpha_primary"],
        krand=L["random_cv_k"], tau=L["group_conformal_tau"], alphas=cfg["alpha"], C=o["logistic_C"],
        kk=o["knn_k"], opca=o["pca_components"], ocv=o["cv_folds"], amax=o["auroc_max"], essmin=o["ess_min"])
    atomic_write_text(md, text)
    atomic_write_json(js, info)
    print(f"  wrote {md}")
    return md


def check_preregistration(cfg, paths, allow_deviation: bool):
    js, _ = prereg_path(paths)
    if not os.path.exists(js):
        raise RuntimeError("No Phase 2 pre-registration found. Run the `setup` step first.")
    pre = json.load(open(js))
    now = config_hash(cfg)
    if pre["config_hash"] != now:
        msg = f"Config changed since pre-registration ({pre['config_hash'][:12]} -> {now[:12]})."
        if not allow_deviation:
            raise RuntimeError(msg + " Re-run with --allow-deviation to proceed; the change will be logged.")
        log = os.path.join(paths.phase0, "DEVIATIONS.json")
        items = json.load(open(log)) if os.path.exists(log) else []
        items.append({"time_utc": _now(), "from": pre["config_hash"], "to": now, "config": cfg["_config_path"]})
        atomic_write_json(log, items)
        print("  WARNING: " + msg + " Logged in phase0/DEVIATIONS.json")
