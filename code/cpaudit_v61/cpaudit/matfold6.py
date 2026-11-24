"""Phase 6, Task 6.3 check: run MatFold itself (Witman & Schindler 2025; MatFold 1.3.3) on all
structures of the frozen table and verify that our holdouts are exactly MatFold's splits.

What MatFold does (read from the paper and the 1.3.3 source):
* labels: SpacegroupAnalyzer(structure, symprec = 0.1) -> space-group symbol ("sgnum") and crystal
  system ("crystalsys"); composition = pymatgen reduced formula;
* create_nested_splits(split_type, n_outer_splits = 0 -> leave-one-label-out, n_inner_splits = 1 ->
  no inner folds), fraction_lower_limit: labels rarer than the limit are never held out and always
  sit in the training set; every non-test entry is training data;
* it writes one train and one test CSV per outer fold, a summary CSV and a JSON file from which the
  splits can be recreated (MatFold.from_json), including MD5 checksums of the input data.

Checks (any failure stops the step): MatFold's labels and compositions equal ours for every row;
MatFold holds out exactly our eligible space groups and the 7 crystal systems; for every fold, its
test set equals our test set row for row and its training set equals the complement.
Output (new folder, never overwritten): phase06_v1/splits_matfold/{sgnum,crystalsys}/ and matfold6.json.
"""
from __future__ import annotations

import glob
import json
import os
import re
import time
import warnings

import numpy as np
import pandas as pd

from .config import atomic_write_json, code_hash, sha256_file
from .data6 import audit_decision, jarvis_structure, load_frozen, load_release, phase6_dirs, release_info

SPLIT_TYPES = (("sgnum", "losgo_mf"), ("crystalsys", "locso_mf"))


def matfold_folds(outdir: str, split_type: str) -> list[tuple[int, np.ndarray, np.ndarray]]:
    """(fold index, train jids, test jids) for every outer fold MatFold wrote in outdir."""
    out = []
    for f in glob.glob(os.path.join(outdir, f"mf.{split_type}.k*_outer.test.csv")):
        k = int(re.search(r"\.k(\d+)_outer\.test\.csv$", f).group(1))
        te = pd.read_csv(f, keep_default_na=False, dtype=str).jid.to_numpy()
        tr = pd.read_csv(f.replace(".test.csv", ".train.csv"), keep_default_na=False, dtype=str).jid.to_numpy()
        out.append((k, tr, te))
    return sorted(out, key=lambda t: t[0])


def compare_folds(folds, jids: np.ndarray, label_of_jid: dict, ours: dict) -> dict:
    """Match each MatFold fold to our holdout by its (single) label and compare row sets exactly."""
    all_j = set(jids)
    res, seen = [], set()
    for k, tr, te in folds:
        labs = {label_of_jid[j] for j in te}
        if len(labs) != 1:
            raise AssertionError(f"MatFold fold k{k}: test set holds {len(labs)} labels, expected 1")
        g = labs.pop()
        if g not in ours:
            raise AssertionError(f"MatFold fold k{k} holds out label {g!r}, which is not one of our held-out groups")
        our_te = set(jids[ours[g]])
        same_test = set(te) == our_te and len(te) == len(our_te)
        same_train = set(tr) == all_j - our_te and len(tr) == len(all_j) - len(our_te)
        res.append(dict(fold=k, label=g, n_test=len(te), n_train=len(tr), same_test=same_test, same_train=same_train))
        seen.add(g)
    missing = sorted(set(ours) - seen, key=str)
    return dict(folds=res, n_folds=len(res), missing_labels=missing,
                all_identical=bool(res) and not missing and all(r["same_test"] and r["same_train"] for r in res))


def run_matfold6(cfg: dict):
    from MatFold import MatFold
    import importlib.metadata as im
    t0 = time.time()
    dirs = phase6_dirs(cfg)
    F = load_frozen(dirs)
    meta = audit_decision(dirs)
    info = release_info(dirs, meta)
    s = cfg["splits"]
    sp_file = os.path.join(dirs["splits"], "splits_v2.json")
    H_file = os.path.join(dirs["splits"], "holdouts_v2.npz")
    for f in (sp_file, H_file):
        if not os.path.exists(f):
            raise FileNotFoundError(f"{f} not found: run splits6 (v6.1) first")
    sp = json.load(open(sp_file))
    if sha256_file(H_file) != sp["holdouts_file"]["sha256"]:
        raise RuntimeError("holdouts_v2.npz differs from the SHA-256 recorded in splits_v2.json")
    L = pd.read_csv(os.path.join(dirs["splits"], "labels.csv.gz"), keep_default_na=False, na_values=[""])
    jids = F.jid.astype(str).to_numpy()
    if not np.array_equal(L.jid.astype(str).to_numpy(), jids):
        raise AssertionError("labels.csv.gz rows differ from the frozen table")
    n = len(F)
    root = os.path.join(dirs["root"], "splits_matfold")
    os.makedirs(root, exist_ok=True)

    print(f"[matfold6] building {n} pymatgen structures and running MatFold {im.version('MatFold')} "
          f"(it re-detects every space group itself; single process, about 5-10 min) ...", flush=True)
    by_jid = {e["jid"]: e["atoms"] for e in load_release(info)}
    bulk = {j: jarvis_structure(by_jid[j]).as_dict() for j in jids}
    del by_jid
    t1 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mf = MatFold(pd.DataFrame({"jid": jids}), bulk, cols_to_keep=["jid"], write_data_checksums=True,
                     seed=int(s["seed"]))
    del bulk
    t_init = time.time() - t1
    m = mf.df.set_index("structureid").loc[jids]
    lab_check = dict(
        n=n,
        sgnum_mismatch=int((m.sgnum.astype(str).to_numpy() != L.sga01_symbol.astype(str).to_numpy()).sum()),
        crystalsys_mismatch=int((m.crystalsys.astype(str).str.lower().to_numpy()
                                 != L.sga01_crys.astype(str).str.lower().to_numpy()).sum()),
        composition_mismatch=int((m.composition.astype(str).to_numpy() != F.reduced_formula.astype(str).to_numpy()).sum()))
    print(f"[matfold6] MatFold init {t_init / 60:.1f} min; label check on all rows: {lab_check}", flush=True)
    if lab_check["sgnum_mismatch"] or lab_check["crystalsys_mismatch"] or lab_check["composition_mismatch"]:
        raise AssertionError(f"MatFold's labels differ from ours: {lab_check}")

    min_size = int(sp["min_group_size"])
    lower = float(sp["matfold_fraction_lower_limit"])
    with np.load(H_file) as Z:
        H = {k: Z[k] for k in Z.files}
    sym2num = dict(zip(L.sga01_symbol.astype(str), L.sga01_number.astype(int)))
    if len(set(sym2num.values())) != len(sym2num):
        raise AssertionError("space-group symbol -> number is not one to one")
    results = {}
    for split_type, fam in SPLIT_TYPES:
        outdir = os.path.join(root, split_type)
        if split_type == "sgnum":
            label_of = dict(zip(jids, L.sga01_number.astype(int)))
            ours = {int(k.split("|")[1]): H[k] for k in H if k.startswith(f"{fam}|") and k.endswith("|test")}
            kw = dict(fraction_lower_limit=lower)
        else:
            label_of = dict(zip(jids, L.sga01_crys.astype(str).str.lower()))
            ours = {k.split("|")[1]: H[k] for k in H if k.startswith(f"{fam}|") and k.endswith("|test")}
            kw = {}
        existing = glob.glob(os.path.join(outdir, f"mf.{split_type}.k*_outer.test.csv"))
        if existing:
            print(f"[matfold6] {outdir} already holds {len(existing)} MatFold folds; verifying them, not rewriting",
                  flush=True)
        else:
            os.makedirs(outdir, exist_ok=True)
            t2 = time.time()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mf.create_nested_splits(split_type, n_inner_splits=1, n_outer_splits=0, write_base_str="mf",
                                        output_dir=outdir, verbose=False, **kw)
            print(f"[matfold6] MatFold wrote {split_type} splits in {(time.time() - t2) / 60:.1f} min", flush=True)
        folds = matfold_folds(outdir, split_type)
        if split_type == "sgnum":   # MatFold's test label is the symbol; ours the number
            label_of = {j: int(label_of[j]) for j in jids}
        cmp_ = compare_folds(folds, jids, label_of, ours)
        js = os.path.join(outdir, f"mf.{split_type}.json")
        cmp_.update(matfold_json=dict(path=js, sha256=sha256_file(js) if os.path.exists(js) else None),
                    options=dict(n_outer_splits=0, n_inner_splits=1, seed=int(s["seed"]), **kw),
                    n_csv_files=len(glob.glob(os.path.join(outdir, "*.csv"))))
        results[split_type] = cmp_
        print(f"[matfold6] {split_type}: {cmp_['n_folds']} MatFold folds; identical to our {fam} holdouts: "
              f"{cmp_['all_identical']}; missing labels: {cmp_['missing_labels']}", flush=True)
        if not cmp_["all_identical"]:
            bad = [r for r in cmp_["folds"] if not (r["same_test"] and r["same_train"])][:10]
            raise AssertionError(f"MatFold {split_type} splits differ from ours: {bad} missing {cmp_['missing_labels']}")
    rec = dict(created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), code_hash=code_hash(),
               matfold_version=im.version("MatFold"), pymatgen_version=im.version("pymatgen"),
               holdouts_v2_sha256=sp["holdouts_file"]["sha256"], min_group_size=min_size,
               fraction_lower_limit=lower, label_check_all_rows=lab_check, splits=results,
               matfold_init_minutes=t_init / 60, minutes=(time.time() - t0) / 60)
    atomic_write_json(os.path.join(root, "matfold6.json"), rec)
    print(f"[matfold6] all checks passed; written {os.path.join(root, 'matfold6.json')} "
          f"in {(time.time() - t0) / 60:.1f} min", flush=True)
    return rec
