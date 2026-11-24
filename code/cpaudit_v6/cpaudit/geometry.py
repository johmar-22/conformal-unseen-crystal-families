"""Identical-input detection for de-duplication (v2).

What the geometry files are
---------------------------
The NOMAD 2018 Kaggle geometry.xyz files hold the idealised Vegard's-law starting structures
(prototype fractional coordinates, composition-interpolated cell); the two targets come from
the relaxed DFT calculations (Sutton et al., npj Comput. Mater. 5, 111, 2019). Two entries
whose geometry files describe the same crystal therefore have the same model input: every
descriptor computed from those files (SOAP included) is the same up to numerical noise.

Why v1 under-counted
--------------------
v1 required equal tabular records AND equal SOAP vectors (rtol 1e-4). Equal-input entries are
often written in a different but symmetry-equivalent setting; their SOAP vectors then differ by
numerical noise above that tolerance, and entries whose targets differ by ~1e-4 eV were not
candidates at all. v1 removed 9 entries.

Rule used here
--------------
1. Identical input: pymatgen StructureMatcher (Ong et al., Comput. Mater. Sci. 68, 314, 2013)
   with tight tolerances (the files are idealised, so true matches are exact), applied within
   blocks of equal (space group, composition, number of atoms). Niggli/primitive reduction makes
   the test independent of the cell setting and atom order.
2. Within an identical-input cluster, entries whose formation energy AND band gap agree within
   `target_tol` are duplicate records of one result: the lowest id is kept.
3. Identical-input entries whose targets differ by more than `target_tol` are CONFLICTING
   records. They are kept (they are part of the public benchmark), flagged, and every summary
   is also reported without them.
An independent check (integer lattice automorphisms + translation, same cell setting only) is
run on every pair found in step 1 and its agreement is reported.
"""
from __future__ import annotations

import itertools
import os
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np

_CAND = None


# ----------------------------------------------------------------------------- reading
def read_xyz(path: str):
    L, cart_rows, species = [], [], []
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p or p[0].startswith("#"):
                continue
            if p[0] == "lattice_vector":
                L.append([float(v) for v in p[1:4]])
            elif p[0] == "atom":
                cart_rows.append([float(v) for v in p[1:4]])   # Cartesian, Angstrom
                species.append(p[4])
    L = np.asarray(L, float)
    frac = np.asarray(cart_rows, float) @ np.linalg.inv(L)       # rows of L are lattice vectors
    return L, frac, np.asarray(species)


def _local_path(entry_id: int, paths) -> str | None:
    for p in (os.path.join(paths.drive_root, "train", str(entry_id), "geometry.xyz"),
              os.path.join(paths.cache, "geometry", f"{entry_id}.xyz")):
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    return None


def fetch_geometries(ids, cfg: dict, paths, workers: int = 16) -> dict:
    """id -> path. Uses NOMAD_Uncer_Quan/train/<id>/geometry.xyz when present, otherwise downloads
    once per Colab session into the local cache (never into Drive)."""
    out, missing = {}, []
    for i in ids:
        p = _local_path(int(i), paths)
        if p:
            out[int(i)] = p
        else:
            missing.append(int(i))
    if missing:
        d = os.path.join(paths.cache, "geometry")
        os.makedirs(d, exist_ok=True)
        print(f"  downloading {len(missing)} geometry files (one time per session)", flush=True)

        def get(i):
            dst = os.path.join(d, f"{i}.xyz")
            try:
                urllib.request.urlretrieve(cfg["paths"]["geometry_url"].format(id=i), dst + ".part")
                os.replace(dst + ".part", dst)
                return i, dst
            except Exception as e:
                return i, f"ERROR {e}"
        with ThreadPoolExecutor(workers) as ex:
            for i, p in ex.map(get, missing):
                if p.startswith("ERROR"):
                    raise RuntimeError(f"could not obtain geometry for id {i}: {p}")
                out[i] = p
    return out


# ----------------------------------------------------------------------------- identical inputs
def identical_input_clusters(df, geo_paths: dict, cfg: dict) -> list[list[int]]:
    """Clusters (sorted NOMAD ids) of entries with the same input crystal structure."""
    import warnings

    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.core import Lattice, Structure
    m = cfg["data"]["dedup"]["matcher"]
    sm = StructureMatcher(ltol=m["ltol"], stol=m["stol"], angle_tol=m["angle_tol"],
                          primitive_cell=True, scale=False, attempt_supercell=False)
    block = ["spacegroup", "percent_atom_al", "percent_atom_ga", "percent_atom_in", "number_of_total_atoms"]
    clusters = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _, g in df.groupby(block, sort=True):
            ids = sorted(int(i) for i in g["id"])
            if len(ids) < 2:
                continue
            structs = []
            for i in ids:
                L, f, s = read_xyz(geo_paths[i])
                st = Structure(Lattice(L), s.tolist(), f)
                st.properties["nomad_id"] = i
                structs.append(st)
            for grp in sm.group_structures(structs):
                if len(grp) > 1:
                    clusters.append(sorted(int(x.properties["nomad_id"]) for x in grp))
    return sorted(clusters)


def split_by_targets(cluster, df_by_id, tol: float):
    """Sub-clusters of an identical-input cluster whose two targets agree within tol."""
    subs = []
    for i in cluster:
        r = df_by_id.loc[i]
        for s in subs:
            ref = df_by_id.loc[s[0]]
            if (abs(r.formation_energy_ev_natom - ref.formation_energy_ev_natom) <= tol
                    and abs(r.bandgap_energy_ev - ref.bandgap_energy_ev) <= tol):
                s.append(i)
                break
        else:
            subs.append([i])
    return subs


# ----------------------------------------------------------------------------- independent check
def _candidates():
    global _CAND
    if _CAND is None:
        mats = np.array(list(itertools.product((-1, 0, 1), repeat=9)), dtype=float).reshape(-1, 3, 3)
        det = np.round(np.linalg.det(mats)).astype(int)
        _CAND = mats[np.abs(det) == 1]
    return _CAND


def lattice_automorphisms(L: np.ndarray, rtol: float = 1e-4) -> np.ndarray:
    """Integer unimodular A with A G A^T = G (G = L L^T): symmetry operations of the cell metric."""
    G = L @ L.T
    C = _candidates()
    GG = np.einsum("nij,jk,nlk->nil", C, G, C)
    ok = np.all(np.abs(GG - G) <= rtol * np.abs(G).max(), axis=(1, 2))
    return C[ok]


def _max_mismatch(fa, sa, fb, sb, L, budget):
    worst = 0.0
    for sp in np.unique(sa):
        d = fa[sa == sp][:, None, :] - fb[sb == sp][None, :, :]
        d -= np.round(d)
        worst = max(worst, float(np.linalg.norm(d @ L, axis=2).min(axis=1).max()))
        if worst > budget:
            return worst
    return worst


def same_structure(ga, gb, tol: float = 1e-3, metric_rtol: float = 1e-4) -> tuple[bool, float]:
    """Same cell setting required. Returns (identical, best max-atom mismatch in Angstrom)."""
    La, fa, sa = ga
    Lb, fb, sb = gb
    if La.shape != Lb.shape or np.abs(La - Lb).max() > tol:
        return False, float("inf")
    if sorted(sa.tolist()) != sorted(sb.tolist()):
        return False, float("inf")
    counts = {s: int((sa == s).sum()) for s in np.unique(sa)}
    rare = min(counts, key=counts.get)
    a0 = fa[sa == rare][0]
    best = float("inf")
    for A in lattice_automorphisms(La, metric_rtol):
        fa_rot, a0_rot = fa @ A, a0 @ A
        for b in fb[sb == rare]:
            best = min(best, _max_mismatch(fa_rot + (b - a0_rot), sa, fb, sb, La, budget=tol))
            if best <= tol:
                return True, best
    return False, best
