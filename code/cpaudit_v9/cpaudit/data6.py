"""Phase 6 data freeze: the JARVIS-DFT analysis table, rebuilt from the pinned release on Drive with
exactly the Task 6.1 audit code, then checked against the audit's own outputs.

The frozen table (one row per de-duplicated entry of the primary scope chosen by the audit) fixes
the row order used by every later Phase 6 step (features, splits, models). It holds identifiers,
labels and the three targets; no target value is summarised here.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np
import pandas as pd

from .audit6 import classify_oxides, dedup, entry_table, group_counts, load_release, scopes
from .config import atomic_write_bytes, atomic_write_json, code_hash, config_hash, get_paths, sha256_file
from .phase5 import _gz

TARGET_COLS = ("formation_energy_peratom", "optb88vdw_bandgap", "mbj_bandgap")


def jarvis_structure(atoms: dict):
    """pymatgen Structure from a JARVIS `atoms` dict (as jarvis-tools' Atoms.pymatgen_converter)."""
    from pymatgen.core import Lattice, Structure
    return Structure(Lattice(atoms["lattice_mat"]), atoms["elements"], atoms["coords"],
                     coords_are_cartesian=bool(atoms["cartesian"]))


def phase6_dirs(cfg: dict) -> dict:
    paths = get_paths(cfg)
    d = {k: os.path.join(paths.out, v) for k, v in
         dict(audit="audit", data="data", features="features", splits="splits", report="report",
              inputs=cfg["jarvis"].get("inputs_subdir", "inputs")).items()}
    for k in ("data", "features", "splits", "report"):
        os.makedirs(d[k], exist_ok=True)
    d["root"] = paths.out
    return d


def audit_decision(dirs: dict) -> dict:
    f = os.path.join(dirs["audit"], "AUDIT_DECISION.json")
    if not os.path.exists(f):
        raise FileNotFoundError(f"{f} not found: run audit6 (cpaudit_v5) first.")
    with open(f) as fh:
        return json.load(fh)


def release_info(dirs: dict, meta: dict) -> dict:
    """The primary release as recorded by the audit, with its SHA-256 re-checked on Drive."""
    rel = meta["releases"][0]
    path = os.path.join(dirs["inputs"], rel["file"] + ".zip")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found (audit6 stores it there).")
    sha = sha256_file(path)
    if sha != rel["sha256"]:
        raise RuntimeError(f"{path}: SHA-256 {sha} differs from the audit's {rel['sha256']}")
    return dict(rel, path=path)


def load_entries_by_jid(info: dict) -> dict:
    return {e["jid"]: e for e in load_release(info)}


def frozen_path(dirs: dict) -> str:
    return os.path.join(dirs["data"], "jarvis_frozen.csv.gz")


def load_frozen(dirs: dict) -> pd.DataFrame:
    f = frozen_path(dirs)
    if not os.path.exists(f):
        raise FileNotFoundError(f"{f} not found: run freeze6 first.")
    return read_table(f)


def read_table(f: str) -> pd.DataFrame:
    # keep_default_na=False: formulas such as "NaN" (sodium nitride) or "NA" must stay strings;
    # only empty fields (how missing targets are written) become NaN
    F = pd.read_csv(f, keep_default_na=False, na_values=[""], float_precision="round_trip")
    for c in ("jid", "reduced_formula", "chemsys", "spg_symbol", "crys", "dimensionality", "oxide_status"):
        if c in F:
            F[c] = F[c].astype(str)
    return F


def run_freeze6(cfg: dict):
    t0 = time.time()
    dirs = phase6_dirs(cfg)
    meta = audit_decision(dirs)
    dec = meta["decision"]
    scope = dec["primary_scope"]
    info = release_info(dirs, meta)
    print(f"[freeze6] release {info['file']} (sha256 {info['sha256'][:12]}), primary scope: {scope}", flush=True)
    entries = load_release(info)
    T = classify_oxides(entry_table(entries), meta["audit_settings"]["other_anion_formers"])
    del entries
    S = dedup(scopes(T)[scope]).copy()
    min_n = int(meta["audit_settings"]["min_group_size"])
    G = group_counts(S, min_n)
    elig_audit = pd.read_csv(os.path.join(dirs["audit"], "audit_eligible_groups_primary.csv"))
    # the freeze must reproduce the audit exactly
    sc = pd.read_csv(os.path.join(dirs["audit"], "audit_scopes.csv"))
    n_audit = int(sc[(sc.release == dec["release"]) & (sc.scope == scope)].n_dedup.iloc[0])
    if len(S) != n_audit:
        raise AssertionError(f"de-duplicated row count {len(S)} differs from the audit ({n_audit})")
    e1 = sorted(G[G.eligible].spg_number.astype(int).tolist())
    e2 = sorted(elig_audit.spg_number.astype(int).tolist())
    if e1 != e2:
        raise AssertionError(f"eligible space groups differ from the audit: {set(e1) ^ set(e2)}")
    n1 = G.set_index("spg_number").loc[e2, "n"].to_numpy()
    if not np.array_equal(n1, elig_audit.set_index("spg_number").loc[e2, "n"].to_numpy()):
        raise AssertionError("eligible group sizes differ from the audit")

    S = S.sort_values("jid_num", kind="mergesort").reset_index(drop=True)
    S.insert(0, "row", np.arange(len(S)))
    S["eligible"] = S.spg_number.isin(e2)
    keep = ["row", "jid", "reduced_formula", "chemsys", "n_elements", "nat", "spg_number", "spg_symbol", "crys",
            "dimensionality", "eligible", "oxide_status", "oxide_oxidation_rule", "oxide_element_rule",
            *TARGET_COLS, "ehull"]
    F = S[keep]
    data = _gz(F)
    atomic_write_bytes(frozen_path(dirs), data)
    sha = sha256_file(frozen_path(dirs))
    back = load_frozen(dirs)
    if not (back.reduced_formula.tolist() == F.reduced_formula.astype(str).tolist()
            and back.jid.tolist() == F.jid.astype(str).tolist()
            and np.allclose(back[list(TARGET_COLS)].to_numpy(float), F[list(TARGET_COLS)].to_numpy(float),
                            equal_nan=True, rtol=0, atol=0)):
        raise AssertionError("the frozen table does not read back identically")
    rec = dict(created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), code_hash=code_hash(),
               config_hash=config_hash(cfg), release=dict(file=info["file"], url=info["url"], sha256=info["sha256"]),
               audit_decision_sha256=sha256_file(os.path.join(dirs["audit"], "AUDIT_DECISION.json")),
               scope=scope, n_rows=int(len(F)), n_eligible_groups=len(e2), n_rows_in_eligible=int(F.eligible.sum()),
               n_with_mbj=int(F.mbj_bandgap.notna().sum()),
               n_missing=dict(formation=int(F.formation_energy_peratom.isna().sum()),
                              optb88vdw=int(F.optb88vdw_bandgap.isna().sum())),
               frozen_table=dict(path=frozen_path(dirs), sha256=sha),
               checks="row count, eligible space groups and their sizes equal the Task 6.1 audit")
    atomic_write_json(os.path.join(dirs["data"], "freeze.json"), rec)
    print(f"[freeze6] {len(F)} entries ({int(F.eligible.sum())} in {len(e2)} eligible space groups); "
          f"TBmBJ present for {rec['n_with_mbj']}; reproduces the audit exactly. "
          f"Written {frozen_path(dirs)} (sha256 {sha[:12]}) in {(time.time() - t0) / 60:.1f} min", flush=True)
    return rec
