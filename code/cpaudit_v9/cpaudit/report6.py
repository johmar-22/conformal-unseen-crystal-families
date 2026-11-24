"""Phase 6 data report: one summary of the freeze, the features and the splits (no targets).
v6.1 writes report/PHASE6_DATA_SUMMARY_v2.md; the v6 summary, if present, is left untouched."""
from __future__ import annotations

import json
import os

import pandas as pd

from .config import atomic_write_text
from .data6 import phase6_dirs


def _load(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found: run the earlier v6 steps first.")
    with open(path) as f:
        return json.load(f)


def run_report6(cfg: dict):
    d = phase6_dirs(cfg)
    fr = _load(os.path.join(d["data"], "freeze.json"))
    fe = _load(os.path.join(d["features"], "features_v1.json"))
    sp = _load(os.path.join(d["splits"], "splits_v2.json"))
    mf = _load(os.path.join(d["root"], "splits_matfold", "matfold6.json"))
    M = pd.read_csv(os.path.join(d["splits"], "holdouts_manifest_v2.csv"))
    fam = sp["families"]
    lo = M[M.kind == fam["primary_losgo"]]
    q = fe["qa"]
    lines = [
        "# Phase 6 data report (Tasks 6.2-6.3; no target used)", "",
        "## Freeze", "",
        f"* {fr['n_rows']} entries of scope `{fr['scope']}` ({fr['n_rows_in_eligible']} in {fr['n_eligible_groups']} "
        f"eligible space groups); TBmBJ gap present for {fr['n_with_mbj']}. Release {fr['release']['file']} "
        f"(sha256 {fr['release']['sha256'][:12]}). Frozen table sha256 {fr['frozen_table']['sha256'][:12]}.",
        f"* Check: {fr['checks']}.", "",
        "## Features", "",
        f"* Matrix {fe['shape'][0]} x {fe['shape'][1]}: composition {fe['blocks']['composition']} (Magpie, Ward 2016), "
        f"density {fe['blocks']['density']}, CrystalNN {fe['blocks']['crystalnn']} (matminer {fe['matminer_version']}, "
        f"pymatgen {fe['pymatgen_version']}); {fe['minutes']:.0f} min.",
        f"* Missing values: {q['columns_with_nan']} columns contain NaN (max share {q['max_nan_fraction']:.4f}); "
        f"{q['n_columns_all_nan']} all-NaN columns; rows with any NaN by block {q['rows_with_nan_by_block']}.",
        f"* Constant columns: {q['n_constant_columns']} {q['constant_columns']}.",
        f"* Structure warnings: {q['n_structure_warnings']} (CrystalNN failed for {q['n_crystalnn_failed']}).",
        f"* Feature-identical entries: {q['n_rows_in_feature_identical_sets']} rows in {q['n_feature_identical_sets']} "
        f"sets ({q['n_feature_identical_sets_spanning_space_groups']} sets span more than one space group).", "",
        "## Splits (MatFold protocol, Witman & Schindler 2025)", "",
        f"* Labels: SpacegroupAnalyzer(symprec = 0.1), as MatFold computes them; primary = MatFold labels "
        f"(`{fam['primary_losgo']}`, `{fam['primary_locso']}`), JARVIS database labels as sensitivity "
        f"(`{fam['sensitivity_losgo']}`, `{fam['sensitivity_locso']}`).",
        f"* LOSGO over {sp['n_losgo']['mf']} space groups with at least {sp['min_group_size']} entries "
        f"(MatFold fraction_lower_limit {sp['matfold_fraction_lower_limit']:.6f}; rarer groups always in training); "
        f"same groups under both label systems and the audit: {sp['eligible_same_in_both_systems']}. "
        f"LOCSO over {sp['n_locso']} crystal systems.",
        f"* MatFold itself, run on all {mf['label_check_all_rows']['n']} structures: label mismatches "
        f"{mf['label_check_all_rows']}; its leave-one-out splits equal ours row by row: "
        f"space groups {mf['splits']['sgnum']['all_identical']} ({mf['splits']['sgnum']['n_folds']} folds), "
        f"crystal systems {mf['splits']['crystalsys']['all_identical']} ({mf['splits']['crystalsys']['n_folds']} folds).",
        f"* MatFold label = JARVIS label for {sp['label_agreement_overall']:.4f} of entries (crystal system "
        f"{sp['crystal_system_agreement_overall']:.4f}); per-group test-set Jaccard, MatFold vs JARVIS labels: median "
        f"{sp['jaccard_test_mf_vs_db']['median']:.3f}, minimum {sp['jaccard_test_mf_vs_db']['min']:.3f} "
        f"(space group {sp['jaccard_test_mf_vs_db']['argmin']}).",
        f"* Near-duplicate twins (same formula and MatFold space group, different JARVIS space group): "
        f"{sp['near_duplicates']['n_rows']} rows; with JARVIS labels "
        f"{sp['near_duplicates']['n_held_out_with_twin_in_training_db_labels']} held-out rows have their twin in "
        f"training; with MatFold labels each such pair stays on one side of every split.",
        f"* Polymorph exclusion (our extension, not part of MatFold) drops a median of "
        f"{lo.n_dropped_polymorph.median():.0f} training entries per held-out group "
        f"(range {lo.n_dropped_polymorph.min()}-{lo.n_dropped_polymorph.max()}).", "",
        lo.drop(columns="kind").round(3).to_markdown(index=False), "",
        M[M.kind == fam["primary_locso"]][["label", "n_test", "n_train", "n_dropped_polymorph",
                                           "n_train_polymorph_excluded"]].to_markdown(index=False), ""]
    txt = "\n".join(lines)
    atomic_write_text(os.path.join(d["report"], "PHASE6_DATA_SUMMARY_v2.md"), txt)
    print(txt)
    return txt
