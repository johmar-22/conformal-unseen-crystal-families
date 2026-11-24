"""Phase 6, Task 6.3: group-holdout splits for JARVIS-DFT (no target is used).

* Leave-one-space-group-out (LOSGO) over the eligible space groups chosen by the Task 6.1 audit,
  and leave-one-crystal-system-out (LOCSO, secondary; 7 groups), as defined in MatFold (Witman &
  Schindler 2025, Digital Discovery 4, 625). Training = every other entry of the frozen table,
  including entries of non-eligible space groups (they are never held out).
* Polymorph-exclusion variant (plan Task 6.3): training additionally drops every entry whose reduced
  formula occurs in the held-out group; the dropped rows are stored per group.
* Group labels are the JARVIS database labels (spg_number, crys), the ones the audit counted. MatFold
  instead re-detects the space group with pymatgen SpacegroupAnalyzer(symprec = 0.1). Both labels are
  stored; we verify on a random subsample, by running MatFold itself, that our re-detection equals
  MatFold's, and report how often MatFold's label differs from the database label, per group.
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


def run_splits6(cfg: dict):
    t0 = time.time()
    dirs = phase6_dirs(cfg)
    F = load_frozen(dirs)
    meta = audit_decision(dirs)
    info = release_info(dirs, meta)
    s = cfg["splits"]
    n = len(F)
    by_jid = {e["jid"]: e["atoms"] for e in load_release(info)}
    lab_file = os.path.join(dirs["splits"], "labels.csv.gz")
    if os.path.exists(lab_file):
        L = pd.read_csv(lab_file)
        if not np.array_equal(L.jid.to_numpy().astype(str), F.jid.to_numpy().astype(str)):
            raise AssertionError("labels.csv.gz does not match the frozen table; delete it and rerun")
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
    chk = matfold_check(F.jid.to_numpy(), by_jid, L, int(s["matfold_check_n"]), int(s["seed"]))
    del by_jid
    print(f"[splits6] MatFold itself on {chk['n']} random structures: same space-group symbol "
          f"{chk['identical_space_group_symbol']:.4f}, same crystal system {chk['identical_crystal_system']:.4f}",
          flush=True)

    spg, crys, form = F.spg_number.to_numpy(), F.crys.astype(str).to_numpy(), F.reduced_formula.astype(str).to_numpy()
    elig = sorted(F.loc[F.eligible, "spg_number"].unique().tolist())
    agree = (L.sga01_number.to_numpy() == spg)
    extra = {g: dict(spg_symbol=str(F.spg_symbol[spg == g].iloc[0]), crys=str(F.crys[spg == g].iloc[0]),
                     share_matfold_label_differs=float(1 - agree[spg == g].mean()),
                     n_matfold_moves_in=int(((L.sga01_number.to_numpy() == g) & (spg != g)).sum()))
             for g in elig}
    losgo = holdout_sets(spg, form, elig)
    locso = holdout_sets(crys, form, sorted(set(crys)))
    M = pd.concat([manifest(losgo, n, "losgo", extra), manifest(locso, n, "locso")], ignore_index=True)
    arrays = {}
    for kind, sets in (("losgo", losgo), ("locso", locso)):
        for g, (test, polyex) in sets.items():
            arrays[f"{kind}|{g}|test"] = test
            arrays[f"{kind}|{g}|polymorph_drop"] = polyex
    out = os.path.join(dirs["splits"], "holdouts_v1.npz")
    atomic_save_npz(out, **arrays)
    atomic_write_csv(os.path.join(dirs["splits"], "holdouts_manifest.csv"), M)
    # K under MatFold's labels (sensitivity only; the audit decision uses the database labels)
    vc = pd.Series(L.sga01_number).value_counts()
    rec = dict(created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), code_hash=code_hash(),
               frozen_table_sha256=sha256_file(os.path.join(dirs["data"], "jarvis_frozen.csv.gz")),
               holdouts_file=dict(path=out, sha256=sha256_file(out)), n_rows=n,
               n_losgo=len(losgo), n_locso=len(locso), matfold_check=chk,
               label_agreement_overall=float(agree.mean()),
               crystal_system_agreement_overall=float((L.sga01_crys.str.lower() == F.crys.str.lower()).mean()),
               K_eligible_under_matfold_labels=int((vc >= int(meta["audit_settings"]["min_group_size"])).sum()),
               polymorph_share_test_formulas=dict(
                   median=float(np.median([len(p) for _, p in losgo.values()]) / n),
                   note="median fraction of all rows dropped from training in the polymorph-exclusion variant"),
               minutes=(time.time() - t0) / 60)
    atomic_write_json(os.path.join(dirs["splits"], "splits_v1.json"), rec)
    print(M[M.kind == "losgo"].drop(columns="kind").round(3).to_string(index=False))
    print(M[M.kind == "locso"].drop(columns="kind").to_string(index=False))
    print(json.dumps({k: v for k, v in rec.items() if k not in ("holdouts_file",)}, indent=2, default=str))
    print(f"[splits6] written {out} in {(time.time() - t0) / 60:.1f} min", flush=True)
    return rec
