"""Data loading, de-duplication, composition groups, Fd-3m sub-families and repeated splits.

Design (Phase 1, Tasks 1.1-1.2 of the revision plan)
----------------------------------------------------
* Exact duplicates (same lattice, composition, cell parameters and DFT targets AND an
  identical SOAP vector) are reduced to one entry (lowest NOMAD id kept). SOAP equality is
  a rotation-, translation- and permutation-invariant check that the two geometries are the
  same structure, which is stronger than comparing tabular columns.
* Composition group = (lattice, Al fraction, Ga fraction) rounded to 0.1, exactly as in the
  v5/v7 leakage analysis. Every split below is grouped on it, so no composition group is ever
  on both sides of a train/test, proper/calibration or CV+ fold boundary.
* All splits are stratified by lattice so each lattice keeps its share in every partition.
"""
from __future__ import annotations

import os
import urllib.request
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

from .config import Paths, atomic_save_npz, sha256_file

LATTICE_NAME = {12: "C2/m", 33: "Pna21", 167: "R-3c", 194: "P63/mmc", 206: "Ia-3", 227: "Fd-3m"}
TARGETS = {"formation": "formation_energy_ev_natom", "bandgap": "bandgap_energy_ev"}
_LATTICE_COLS = ["lattice_vector_1_ang", "lattice_vector_2_ang", "lattice_vector_3_ang",
                 "lattice_angle_alpha_degree", "lattice_angle_beta_degree", "lattice_angle_gamma_degree"]
_DUP_KEY = (["spacegroup", "percent_atom_al", "percent_atom_ga", "percent_atom_in"]
            + _LATTICE_COLS + list(TARGETS.values()))


# ----------------------------------------------------------------------------- inputs
def load_train_csv(cfg: dict, paths: Paths) -> tuple[pd.DataFrame, str]:
    local = os.path.join(paths.inputs, "train.csv")
    if not os.path.exists(local):
        url = cfg["paths"]["train_csv_url"]
        print(f"  downloading train.csv from {url}")
        urllib.request.urlretrieve(url, local + ".part")
        os.replace(local + ".part", local)
    df = pd.read_csv(local)
    return df, sha256_file(local)


def feature_cache_path(cfg: dict, paths: Paths) -> str:
    return os.path.join(paths.cache, "soap_float32.npy")


def ensure_feature_cache(cfg: dict, paths: Paths) -> str:
    """Copy the SOAP matrix from Drive to a local float32 .npy once per Colab session.
    Workers then memory-map it instead of each decompressing 354 MB from Drive."""
    dst = feature_cache_path(cfg, paths)
    if os.path.exists(dst):
        return dst
    src = paths.drive_file(cfg["paths"]["soap_npz"])
    if not os.path.exists(src):
        raise FileNotFoundError(f"SOAP feature file not found: {src}")
    print(f"  caching SOAP features locally (one time per session): {src}")
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
    """'SG{sg}_Al{a}_Ga{g}' with np.round, reproducing the v5/v7 group labels exactly."""
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


MAX_DUP_SOAP_RELDIFF = {"value": 0.0}   # largest SOAP difference inside accepted duplicate groups


def find_duplicates(df: pd.DataFrame, X: np.ndarray) -> list[list[int]]:
    """Groups of row positions that are tabular duplicates AND have identical SOAP vectors.
    Tolerance rtol=1e-4 absorbs floating-point noise from a different atom order or orientation
    of the same structure; distinct structures differ by orders of magnitude more."""
    groups = []
    for _, sub in df.groupby(_DUP_KEY, sort=False):
        if len(sub) < 2:
            continue
        rows = sub.index.to_list()
        ref = np.asarray(X[rows[0]], dtype=np.float64)
        same = [rows[0]]
        for r in rows[1:]:
            v = np.asarray(X[r], dtype=np.float64)
            if np.allclose(v, ref, rtol=1e-4, atol=1e-6):
                same.append(r)
                rel = float(np.max(np.abs(v - ref) / (np.abs(ref) + 1e-6)))
                MAX_DUP_SOAP_RELDIFF["value"] = max(MAX_DUP_SOAP_RELDIFF["value"], rel)
        if len(same) > 1:
            groups.append(sorted(same, key=lambda r: df.loc[r, "id"]))
    return groups


# ----------------------------------------------------------------------------- dataset
@dataclass
class Dataset:
    df: pd.DataFrame            # de-duplicated, index 0..n-1
    orig_row: np.ndarray        # row of each entry in train.csv / the SOAP matrix
    sg: np.ndarray
    sub: np.ndarray             # lattice with Fd-3m split into 227a / 227b
    group: np.ndarray           # composition group
    y: dict                     # target name -> array
    dup_groups: list            # duplicate groups as NOMAD ids (removed entries included)

    @property
    def n(self) -> int:
        return len(self.df)


def build_dataset(cfg: dict, paths: Paths, X: np.ndarray | None = None) -> Dataset:
    df, _ = load_train_csv(cfg, paths)
    if len(df) != cfg["data"]["expected_rows"]:
        raise AssertionError(f"train.csv has {len(df)} rows, expected {cfg['data']['expected_rows']}")
    if X is None:
        X = load_features(cfg, paths)
    for c in TARGETS.values():
        if not np.isfinite(df[c]).all():
            raise AssertionError(f"non-finite values in {c}")
    dups = find_duplicates(df, X)
    drop = sorted(r for g in dups for r in g[1:])
    keep = np.setdiff1d(np.arange(len(df)), drop)
    d = df.iloc[keep].reset_index(drop=True)
    return Dataset(
        df=d,
        orig_row=keep.astype(np.int64),
        sg=d["spacegroup"].to_numpy(int),
        sub=fd3m_subfamily(d, cfg["data"]["fd3m_gamma_split_deg"]),
        group=composition_groups(d, cfg["data"]["composition_decimals"]),
        y={k: d[v].to_numpy(np.float64) for k, v in TARGETS.items()},
        dup_groups=[[int(df.loc[r, "id"]) for r in g] for g in dups],
    )


# ----------------------------------------------------------------------------- splits
def _first_fold(n_splits: int, strat, groups, seed: int):
    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    tr, te = next(iter(sgkf.split(np.zeros(len(strat)), strat, groups)))
    return np.sort(tr), np.sort(te)


def _fold_ids(n_splits: int, strat, groups, seed: int) -> np.ndarray:
    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    f = np.full(len(strat), -1, dtype=np.int64)
    for k, (_, te) in enumerate(sgkf.split(np.zeros(len(strat)), strat, groups)):
        f[te] = k
    assert (f >= 0).all()
    return f


def make_split(ds: Dataset, cfg: dict, r: int) -> dict:
    """Repeat r: outer train/test, proper/calibration and CV+ fold labels (positions in ds)."""
    s = cfg["splits"]
    base = int(s["seed_base"])
    train, test = _first_fold(s["outer_folds"], ds.sg, ds.group, base + r)
    p_rel, c_rel = _first_fold(s["cal_folds"], ds.sg[train], ds.group[train], base + 1000 + r)
    proper, cal = train[p_rel], train[c_rel]
    cv_random = _fold_ids(cfg["cvplus"]["k_random"], ds.sg[train], ds.group[train], base + 2000 + r)
    lattices = np.unique(ds.sg[train])
    cv_lattice = np.searchsorted(lattices, ds.sg[train]).astype(np.int64)
    # --- guarantees the analysis relies on ---
    assert not set(ds.group[train]) & set(ds.group[test]), "composition group on both sides of train/test"
    assert not set(ds.group[proper]) & set(ds.group[cal]), "composition group in proper and calibration"
    for k in np.unique(cv_random):
        assert not set(ds.group[train][cv_random == k]) & set(ds.group[train][cv_random != k])
    assert len(np.intersect1d(train, test)) == 0 and len(np.intersect1d(proper, cal)) == 0
    return dict(repeat=r, train=train, test=test, proper=proper, cal=cal,
                cv_random=cv_random, cv_lattice=cv_lattice)


def save_split_manifest(ds: Dataset, split: dict, path: str) -> None:
    ids = ds.df["id"].to_numpy(int)
    atomic_save_npz(path, train_id=ids[split["train"]], test_id=ids[split["test"]],
                    proper_id=ids[split["proper"]], cal_id=ids[split["cal"]],
                    cv_random=split["cv_random"], cv_lattice=split["cv_lattice"])


def split_summary(ds: Dataset, split: dict) -> pd.DataFrame:
    rows = []
    for lat in sorted(np.unique(ds.sg)):
        row = {"repeat": split["repeat"], "lattice": LATTICE_NAME[lat], "sg": lat}
        for part in ("train", "test", "proper", "cal"):
            row[f"n_{part}"] = int((ds.sg[split[part]] == lat).sum())
        rows.append(row)
    return pd.DataFrame(rows)
