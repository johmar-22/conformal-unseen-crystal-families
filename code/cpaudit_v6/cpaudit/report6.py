"""Phase 6 data report: one summary of the freeze, the features and the splits (no targets)."""
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
    sp = _load(os.path.join(d["splits"], "splits_v1.json"))
    M = pd.read_csv(os.path.join(d["splits"], "holdouts_manifest.csv"))
    lo = M[M.kind == "losgo"]
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
        "## Splits", "",
        f"* LOSGO over {sp['n_losgo']} eligible space groups (JARVIS labels); LOCSO over {sp['n_locso']} crystal systems.",
        f"* MatFold check on {sp['matfold_check']['n']} random structures: our re-detection = MatFold's labels for "
        f"{sp['matfold_check']['identical_space_group_symbol']:.4f} (space group) and "
        f"{sp['matfold_check']['identical_crystal_system']:.4f} (crystal system).",
        f"* MatFold's label (symprec 0.1) equals the JARVIS label for {sp['label_agreement_overall']:.4f} of entries "
        f"(crystal system {sp['crystal_system_agreement_overall']:.4f}); eligible groups under MatFold labels: "
        f"{sp['K_eligible_under_matfold_labels']}.",
        f"* Polymorph exclusion drops a median of {lo.n_dropped_polymorph.median():.0f} training entries per held-out "
        f"group (range {lo.n_dropped_polymorph.min()}-{lo.n_dropped_polymorph.max()}).", "",
        lo.drop(columns="kind").round(3).to_markdown(index=False), "",
        M[M.kind == "locso"].drop(columns=[c for c in M.columns if c in ("kind", "spg_symbol", "crys",
                                                                       "share_matfold_label_differs",
                                                                       "n_matfold_moves_in")])
        .to_markdown(index=False), ""]
    txt = "\n".join(lines)
    atomic_write_text(os.path.join(d["report"], "PHASE6_DATA_SUMMARY.md"), txt)
    print(txt)
    return txt
