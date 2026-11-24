"""Data loading, de-duplication (v2 rule), composition groups and Fd-3m sub-families.

* De-duplication (v2): entries with the same input crystal structure (pymatgen StructureMatcher on
  the Vegard geometry files) AND both targets within `target_tol` are one record; the lowest
  NOMAD id is kept. Same-input entries whose targets differ by more are CONFLICTING records:
  kept, flagged in `Dataset.conflict`, and reported separately. See geometry.py for the reasoning.
  The result is computed once by `setup` and stored in phase0/dedup_v2.json; every later step
  reads that file, so all steps use the identical dataset.
* Composition group = (lattice, Al fraction, Ga fraction) rounded to 0.1, exactly as in v1/v5/v7.
  Every member of an identical-input cluster has the same lattice and composition, hence the same
  composition group, so duplicates can never straddle a group-disjoint split.
"""
from __future__ import annotations

import json
import os
import shutil
import urllib.request
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .config import Paths, atomic_write_json, sha256_file

LATTICE_NAME = {12: "C2/m", 33: "Pna21", 167: "R-3c", 194: "P63/mmc", 206: "Ia-3", 227: "Fd-3m"}
LATTICE_CODE = {v: k for k, v in LATTICE_NAME.items()}
LATTICE_ORDER = ["C2/m", "Pna21", "R-3c", "P63/mmc", "Ia-3", "Fd-3m"]
TARGETS = {"formation": "formation_energy_ev_natom", "bandgap": "bandgap_energy_ev"}
SUBFAMILY_NAME = {"227a": "Fd-3m (gamma~30)", "227b": "Fd-3m (gamma~120)"}


# ----------------------------------------------------------------------------- inputs
def load_train_csv(cfg: dict, paths: Paths) -> tuple[pd.DataFrame, str]:
    local = os.path.join(paths.inputs, "train.csv")
    if not os.path.exists(local):
        p1 = os.path.join(paths.phase2, "inputs", "train.csv")
        if os.path.exists(p1):
            shutil.copyfile(p1, local)
        else:
            url = cfg["paths"]["train_csv_url"]
            print(f"  downloading train.csv from {url}")
            urllib.request.urlretrieve(url, local + ".part")
            os.replace(local + ".part", local)
    return pd.read_csv(local), sha256_file(local)


def feature_cache_path(cfg: dict, paths: Paths) -> str:
    return os.path.join(paths.cache, "soap_float32.npy")


def ensure_feature_cache(cfg: dict, paths: Paths) -> str:
    """Copy the SOAP matrix from Drive to a local float32 .npy once per Colab session.
    Workers memory-map it instead of each decompressing 354 MB from Drive."""
    dst = feature_cache_path(cfg, paths)
    if os.path.exists(dst):
        return dst
    src = paths.drive_file(cfg["paths"]["soap_npz"])
    if not os.path.exists(src):
        raise FileNotFoundError(f"SOAP feature file not found: {src}")
    print(f"  caching SOAP features locally (one time per session): {src}", flush=True)
    with np.load(src) as z:
        X = np.asarray(z["features"], dtype=np.float32)
    exp = (cfg["data"]["expected_rows"], cfg["data"]["expected_features"])
    if X.shape != exp:
        raise AssertionError(f"SOAP matrix shape {X.shape} != expected {exp}")
    if not np.isfinite(X).all():
        raise AssertionError("SOAP matrix contains non-finite values")
    np.save(dst + ".part.npy", X)
    os.replace(dst + ".part.npy", dst)
    return dst


def load_features(cfg: dict, paths: Paths) -> np.ndarray:
    return np.load(ensure_feature_cache(cfg, paths), mmap_mode="r")


# ----------------------------------------------------------------------------- annotations
def composition_groups(df: pd.DataFrame, decimals: int = 1) -> np.ndarray:
    """'SG{sg}_Al{a}_Ga{g}' with np.round, reproducing the v1/v5/v7 group labels exactly."""
    al = np.round(df["percent_atom_al"].to_numpy(np.float64), decimals)
    ga = np.round(df["percent_atom_ga"].to_numpy(np.float64), decimals)
    sg = df["spacegroup"].to_numpy(int)
    return np.array([f"SG{s}_Al{a}_Ga{g}" for s, a, g in zip(sg, al, ga)])


def fd3m_subfamily(df: pd.DataFrame, gamma_split: float) -> np.ndarray:
    """Lattice label with Fd-3m split by cell setting: '227a' (gamma < split) / '227b'."""
    sg = df["spacegroup"].to_numpy(int)
    gam = df["lattice_angle_gamma_degree"].to_numpy(float)
    out = sg.astype(str).astype(object)
    m = sg == 227
    out[m & (gam < gamma_split)] = "227a"
    out[m & (gam >= gamma_split)] = "227b"
    return out.astype(str)


def id_table(cfg: dict, paths: Paths) -> pd.DataFrame:
    """All 2,400 entries indexed by NOMAD id, with lattice, composition group and sub-family."""
    df, _ = load_train_csv(cfg, paths)
    t = pd.DataFrame({"id": df["id"].to_numpy(int), "sg": df["spacegroup"].to_numpy(int),
                      "group": composition_groups(df, cfg["data"]["composition_decimals"]),
                      "sub": fd3m_subfamily(df, cfg["data"]["fd3m_gamma_split_deg"])})
    return t.set_index("id", drop=False)


# ----------------------------------------------------------------------------- de-duplication
def dedup_v2(cfg: dict, paths: Paths, df: pd.DataFrame) -> dict:
    """Identical-input clusters, duplicate groups and conflicting records (cached in phase0)."""
    f = os.path.join(paths.phase0, "dedup_v2.json")
    if os.path.exists(f):
        with open(f) as fh:
            return json.load(fh)
    from .geometry import (fetch_geometries, identical_input_clusters, read_xyz, same_structure,
                           split_by_targets)
    d = cfg["data"]["dedup"]
    geo = fetch_geometries(df["id"].to_numpy(int), cfg, paths)
    print("  matching input structures (pymatgen StructureMatcher) ...", flush=True)
    clusters = identical_input_clusters(df, geo, cfg)
    by_id = df.set_index("id", drop=False)
    dup_groups, conflicts, check = [], [], []
    for c in clusters:
        subs = split_by_targets(c, by_id, d["target_tol"])
        dup_groups += [s for s in subs if len(s) > 1]
        if len(subs) > 1:
            conflicts.append(subs)
        g0 = read_xyz(geo[c[0]])
        for j in c[1:]:
            ok, mis = same_structure(g0, read_xyz(geo[j]), tol=d["check_tol_ang"])
            check.append(dict(a=c[0], b=j, same_setting_match=bool(ok), mismatch_ang=float(mis)))
    res = dict(rule=d["rule"], matcher=d["matcher"], target_tol=d["target_tol"],
               identical_input_clusters=clusters, duplicate_groups=dup_groups,
               removed_ids=sorted(i for g in dup_groups for i in g[1:]),
               conflicting_clusters=conflicts,
               conflicting_ids=sorted(i for c in conflicts for s in c for i in s),
               independent_check=check)
    atomic_write_json(f, res)
    return res


# ----------------------------------------------------------------------------- dataset
@dataclass
class Dataset:
    df: pd.DataFrame            # de-duplicated, index 0..n-1
    orig_row: np.ndarray        # row of each entry in train.csv / the SOAP matrix
    sg: np.ndarray
    sub: np.ndarray             # lattice with Fd-3m split into 227a / 227b
    group: np.ndarray           # composition group
    y: dict                     # target name -> array
    conflict: np.ndarray        # True for members of a conflicting identical-input cluster
    dedup: dict                 # full de-duplication record

    @property
    def n(self) -> int:
        return len(self.df)

    @property
    def ids(self) -> np.ndarray:
        return self.df["id"].to_numpy(int)


def build_dataset(cfg: dict, paths: Paths) -> Dataset:
    df, h = load_train_csv(cfg, paths)
    if len(df) != cfg["data"]["expected_rows"]:
        raise AssertionError(f"train.csv has {len(df)} rows, expected {cfg['data']['expected_rows']}")
    want = cfg["data"].get("train_csv_sha256")
    if want and h != want:
        raise AssertionError(f"train.csv SHA-256 {h[:12]} differs from the Phase 1 file {want[:12]}")
    if not (df["id"].to_numpy() == np.arange(1, len(df) + 1)).all():
        raise AssertionError("train.csv ids are not 1..2400 in order; SOAP row alignment would break")
    for c in TARGETS.values():
        if not np.isfinite(df[c]).all():
            raise AssertionError(f"non-finite values in {c}")
    dd = dedup_v2(cfg, paths, df)
    keep = np.setdiff1d(np.arange(len(df)), np.asarray(dd["removed_ids"], int) - 1)
    d = df.iloc[keep].reset_index(drop=True)
    return Dataset(
        df=d, orig_row=keep.astype(np.int64), sg=d["spacegroup"].to_numpy(int),
        sub=fd3m_subfamily(d, cfg["data"]["fd3m_gamma_split_deg"]),
        group=composition_groups(d, cfg["data"]["composition_decimals"]),
        y={k: d[v].to_numpy(np.float64) for k, v in TARGETS.items()},
        conflict=np.isin(d["id"].to_numpy(int), dd["conflicting_ids"]), dedup=dd)
