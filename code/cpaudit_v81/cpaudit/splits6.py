"""Phase 6, Task 6.3: group-holdout splits for JARVIS-DFT (no target is used).

MatFold protocol (Witman & Schindler 2025, Digital Discovery 4, 625; MatFold 1.3.3 source):
* Split labels are re-detected from each structure with pymatgen SpacegroupAnalyzer(symprec = 0.1):
  space-group symbol ("sgnum") and crystal system ("crystalsys"). These MatFold labels ("mf") are the
  PRIMARY labels here. The JARVIS database labels (spg_number, crys; "db") are kept as a sensitivity
  family, because 2.2 % of entries change space group between the two (mostly P1, Pm, P2/m, Cm).
* Leave-one-out outer folds (MatFold n_outer_splits = 0): one held-out label per fold. Labels with
  fewer than the audit's minimum group size (200 entries) are never held out and are always in the
  training set, which is MatFold's fraction_lower_limit rule (limit = (200 - 0.5) / N). Crystal
  systems: all 7, no limit. Every other entry is training data (MatFold's default-train rule).
* Our extension, not part of MatFold: a polymorph-exclusion variant in which training additionally
  drops every entry whose reduced formula occurs in the held-out group (stored per group).
* matfold6.py runs MatFold itself on all 81,884 structures and checks these holdouts row by row.
"""
from __future__ import annotations

import json
import os
import tempfile
import time
import warnings

import numpy as np
import pandas as pd

from .config import atomic_save_npz, atomic_write_csv, atomic_write_json, code_hash, sha256_file
from .data6 import audit_decision, jarvis_structure, load_frozen, load_release, phase6_dirs, release_info


# ----------------------------------------------------------------------------- labels
def _sga(atoms: dict, symprec: float = 0.1):
    """(number, symbol, crystal system) from pymatgen SpacegroupAnalyzer, as MatFold computes them."""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    warnings.simplefilter("ignore")
    try:
        a = SpacegroupAnalyzer(jarvis_structure(atoms), symprec)
        return int(a.get_space_group_number()), str(a.get_space_group_symbol()), str(a.get_crystal_system())
    except Exception as e:  # pragma: no cover
        return -1, f"error: {type(e).__name__}", "unknown"


def _sga01(atoms):
    return _sga(atoms, 0.1)


def matfold_check(jids, by_jid, ours: pd.DataFrame, n: int, seed: int) -> dict:
    """Run MatFold itself on a random subsample and compare its labels with our re-detection."""
    from MatFold import MatFold
    rng = np.random.default_rng(seed)
    pick = rng.choice(len(jids), size=min(n, len(jids)), replace=False)
    sub = [str(jids[i]) for i in sorted(pick)]
    bulk = {j: jarvis_structure(by_jid[j]).as_dict() for j in sub}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mf = MatFold(pd.DataFrame({"jid": sub}), bulk, write_data_checksums=False, cols_to_keep=["jid"])
    m = mf.df.set_index("structureid")
    o = ours.set_index("jid").loc[sub]
    same_sg = (m.loc[sub, "sgnum"].astype(str).to_numpy() == o.sga01_symbol.astype(str).to_numpy())
    same_cs = (m.loc[sub, "crystalsys"].astype(str).to_numpy() == o.sga01_crys.astype(str).to_numpy())
    return dict(n=len(sub), seed=seed, matfold_version=_version("MatFold"),
                identical_space_group_symbol=float(same_sg.mean()), identical_crystal_system=float(same_cs.mean()))


def _version(mod):
    try:
        import importlib.metadata as im
        return im.version(mod)
    except Exception:
        return "unknown"


# ----------------------------------------------------------------------------- splits
def holdout_sets(labels: np.ndarray, formulas: np.ndarray, groups) -> dict:
    """For each held-out label: test rows, and rows dropped from training in the polymorph-exclusion
    variant (other groups' rows whose reduced formula occurs in the held-out group)."""
    out = {}
    for g in groups:
        test = np.where(labels == g)[0]
        f = set(formulas[test])
        other = labels != g
        polyex = np.where(other & np.isin(formulas, list(f)))[0]
        out[g] = (test.astype(np.int32), polyex.astype(np.int32))
    return out


def manifest(sets: dict, n: int, kind: str, extra: dict | None = None) -> pd.DataFrame:
    rows = []
    for g, (test, polyex) in sets.items():
        r = dict(kind=kind, label=g, n_test=len(test), n_train=n - len(test), n_dropped_polymorph=len(polyex),
                 n_train_polymorph_excluded=n - len(test) - len(polyex))
        if extra and g in extra:
            r.update(extra[g])
        rows.append(r)
    return pd.DataFrame(rows)


def label_table(F: pd.DataFrame, L: pd.DataFrame) -> dict:
    """Group labels per label system: (space-group number, crystal system), lower-case crystal system."""
    return {"mf": (L.sga01_number.to_numpy().astype(int), L.sga01_crys.astype(str).str.lower().to_numpy()),
            "db": (F.spg_number.to_numpy().astype(int), F.crys.astype(str).str.lower().to_numpy())}


def eligible_groups(spg: np.ndarray, min_size: int) -> list:
    vc = pd.Series(spg).value_counts()
    return sorted(int(g) for g in vc[vc >= min_size].index)


def run_splits6(cfg: dict):
    t0 = time.time()
    dirs = phase6_dirs(cfg)
    F = load_frozen(dirs)
    meta = audit_decision(dirs)
    info = release_info(dirs, meta)
    s = cfg["splits"]
    n = len(F)
    min_size = int(meta["audit_settings"]["min_group_size"])
    primary = str(s.get("primary_labels", "mf"))
    if primary not in ("mf", "db"):
        raise ValueError("splits.primary_labels must be 'mf' (MatFold) or 'db' (JARVIS database)")
    by_jid = {e["jid"]: e["atoms"] for e in load_release(info)}
    lab_file = os.path.join(dirs["splits"], "labels.csv.gz")
    if os.path.exists(lab_file):
        L = pd.read_csv(lab_file, keep_default_na=False, na_values=[""])
        if not np.array_equal(L.jid.to_numpy().astype(str), F.jid.to_numpy().astype(str)):
            raise AssertionError("labels.csv.gz does not match the frozen table; delete it and rerun")
        print(f"[splits6] reusing {lab_file}", flush=True)
    else:
        import multiprocessing as mp
        from .features6 import _workers
        ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
        W = _workers(cfg)
        print(f"[splits6] re-detecting space groups as MatFold does (SpacegroupAnalyzer, symprec "
              f"{s['symprec_matfold']}) for {n} structures, {W} processes ...", flush=True)
        if float(s["symprec_matfold"]) != 0.1:
            raise ValueError("MatFold 1.3 uses symprec = 0.1; keep it for comparability")
        with ctx.Pool(W) as pool:
            res = pool.map(_sga01, [by_jid[j] for j in F.jid], chunksize=200)
        L = pd.DataFrame(dict(row=F.row, jid=F.jid, spg_jarvis=F.spg_number, crys_jarvis=F.crys,
                              sga01_number=[r[0] for r in res], sga01_symbol=[r[1] for r in res],
                              sga01_crys=[r[2] for r in res]))
        from .phase5 import _gz
        from .config import atomic_write_bytes
        atomic_write_bytes(lab_file, _gz(L))
    if (L.sga01_number < 0).any():
        raise RuntimeError(f"SpacegroupAnalyzer failed for {(L.sga01_number < 0).sum()} structures")
    chk = matfold_check(F.jid.to_numpy(), by_jid, L, int(s["matfold_check_n"]), int(s["seed"]))
    del by_jid
    print(f"[splits6] MatFold itself on {chk['n']} random structures: same space-group symbol "
          f"{chk['identical_space_group_symbol']:.4f}, same crystal system {chk['identical_crystal_system']:.4f}"
          f" (all {n} structures are checked by matfold6)", flush=True)

    form = F.reduced_formula.astype(str).to_numpy()
    labs = label_table(F, L)
    elig_audit = sorted(int(g) for g in F.loc[F.eligible, "spg_number"].unique())
    arrays, manifests, groups = {}, [], {}
    for sysname, (spg, crys) in labs.items():
        elig = eligible_groups(spg, min_size)
        groups[sysname] = elig
        other = labs["db" if sysname == "mf" else "mf"][0]
        extra = {g: dict(share_label_differs_other_system=float((other[spg == g] != g).mean()),
                         n_moves_in_vs_other_system=int(((spg == g) & (other != g)).sum()))
                 for g in elig}
        sym = L.sga01_symbol.to_numpy() if sysname == "mf" else F.spg_symbol.astype(str).to_numpy()
        for g in elig:
            extra[g].update(spg_symbol=str(sym[spg == g][0]), crys=str(crys[spg == g][0]))
        losgo = holdout_sets(spg, form, elig)
        locso = holdout_sets(crys, form, sorted(set(crys)))
        for kind, sets in ((f"losgo_{sysname}", losgo), (f"locso_{sysname}", locso)):
            manifests.append(manifest(sets, n, kind, extra if kind.startswith("losgo") else None))
            for g, (test, polyex) in sets.items():
                arrays[f"{kind}|{g}|test"] = test
                arrays[f"{kind}|{g}|polymorph_drop"] = polyex
    M = pd.concat(manifests, ignore_index=True)
    # test-set overlap between the two label systems, per space group eligible in both
    both = sorted(set(groups["mf"]) & set(groups["db"]))
    jac = {g: len(np.intersect1d(arrays[f"losgo_mf|{g}|test"], arrays[f"losgo_db|{g}|test"])) /
              len(np.union1d(arrays[f"losgo_mf|{g}|test"], arrays[f"losgo_db|{g}|test"])) for g in both}
    M["jaccard_test_mf_vs_db"] = [jac.get(g, np.nan) if k.startswith("losgo") else np.nan
                                  for k, g in zip(M.kind, M.label)]
    # near-duplicate pairs: same formula and same MatFold space group but different JARVIS space group
    key = pd.Series(form) + "|" + pd.Series(labs["mf"][0]).astype(str)
    dup = key.duplicated(keep=False).to_numpy()
    out = os.path.join(dirs["splits"], "holdouts_v2.npz")
    if os.path.exists(out):     # never overwrite: a rerun must reproduce the stored holdouts exactly
        with np.load(out) as Z:
            same = set(Z.files) == set(arrays) and all(np.array_equal(Z[k], arrays[k]) for k in arrays)
        if not same:
            raise RuntimeError(f"{out} exists and differs from the holdouts rebuilt now; move it away to rebuild")
        print(f"[splits6] {out} exists and is identical to the rebuilt holdouts; kept", flush=True)
    else:
        atomic_save_npz(out, **arrays)
    atomic_write_csv(os.path.join(dirs["splits"], "holdouts_manifest_v2.csv"), M)
    rec = dict(created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), code_hash=code_hash(),
               protocol="MatFold (Witman & Schindler 2025) leave-one-label-out, labels from "
                        "SpacegroupAnalyzer(symprec=0.1); fraction_lower_limit = (min_group_size - 0.5) / N",
               primary_labels=primary, families=dict(primary_losgo=f"losgo_{primary}", primary_locso=f"locso_{primary}",
                                                     sensitivity_losgo=f"losgo_{'db' if primary == 'mf' else 'mf'}",
                                                     sensitivity_locso=f"locso_{'db' if primary == 'mf' else 'mf'}"),
               frozen_table_sha256=sha256_file(os.path.join(dirs["data"], "jarvis_frozen.csv.gz")),
               labels_sha256=sha256_file(lab_file),
               holdouts_file=dict(path=out, sha256=sha256_file(out)), n_rows=n, min_group_size=min_size,
               matfold_fraction_lower_limit=(min_size - 0.5) / n,
               eligible_groups=groups, eligible_audit_db_labels=elig_audit,
               eligible_same_in_both_systems=groups["mf"] == groups["db"] == elig_audit,
               n_losgo=dict(mf=len(groups["mf"]), db=len(groups["db"])), n_locso=7, matfold_check=chk,
               label_agreement_overall=float((labs["mf"][0] == labs["db"][0]).mean()),
               crystal_system_agreement_overall=float((labs["mf"][1] == labs["db"][1]).mean()),
               jaccard_test_mf_vs_db=dict(median=float(np.median(list(jac.values()))),
                                          min=float(min(jac.values())),
                                          argmin=int(min(jac, key=jac.get))),
               near_duplicates=dict(
                   n_rows=int(dup.sum()),
                   n_held_out_with_twin_in_training_db_labels=int((dup & np.isin(labs["db"][0], groups["db"])).sum()),
                   n_held_out_with_twin_in_training_mf_labels=0,
                   note="rows that share reduced formula and MatFold space group with another row (the JARVIS "
                        "dedup keeps both because their JARVIS space groups differ). Under the JARVIS labels a "
                        "held-out row of such a pair has its twin in training; under MatFold labels both sit in "
                        "the same group by construction. The polymorph-exclusion variant removes the twins."),
               polymorph_share_test_formulas=dict(
                   median=float(M.loc[M.kind == f"losgo_{primary}", "n_dropped_polymorph"].median() / n),
                   note="median fraction of all rows dropped from training in the polymorph-exclusion variant"),
               minutes=(time.time() - t0) / 60)
    if not rec["eligible_same_in_both_systems"]:
        print("[splits6] NOTE: the eligible space groups differ between label systems or from the audit; "
              "see eligible_groups in splits_v2.json", flush=True)
    atomic_write_json(os.path.join(dirs["splits"], "splits_v2.json"), rec)
    cols = ["label", "spg_symbol", "crys", "n_test", "n_dropped_polymorph", "share_label_differs_other_system",
            "jaccard_test_mf_vs_db"]
    print(M[M.kind == f"losgo_{primary}"][cols].round(3).to_string(index=False))
    print(M[M.kind.str.startswith("locso")][["kind", "label", "n_test", "n_dropped_polymorph"]].to_string(index=False))
    print(json.dumps({k: v for k, v in rec.items() if k not in ("holdouts_file", "eligible_groups",
                                                                 "eligible_audit_db_labels")}, indent=2, default=str))
    print(f"[splits6] written {out} in {(time.time() - t0) / 60:.1f} min", flush=True)
    return rec
