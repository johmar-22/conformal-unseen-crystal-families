"""Phase 6, Task 6.2: features for the JARVIS-DFT analysis table (no target is used).

Blocks (all from matminer, Ward et al. 2018, Comput. Mater. Sci. 152, 60):
  composition  the full Magpie attribute set of Ward et al. (2016, npj Comput. Mater. 2, 16028):
               ElementProperty('magpie') 132 + Stoichiometry 6 + ValenceOrbital(frac) 4 +
               IonProperty 3 = 145 features, computed once per reduced formula.
  density      DensityFeatures: density, volume per atom, packing fraction (3). If the packing
               fraction cannot be computed (an element without an atomic radius, e.g. a noble gas),
               density and volume per atom are computed directly and the packing fraction is NaN.
  crystalnn    SiteStatsFingerprint.from_preset('CrystalNNFingerprint_ops'): mean and standard
               deviation over sites of the CrystalNN coordination / local-order fingerprint
               (Zimmermann & Jain 2020, RSC Adv. 10, 6063), 122 features.
GlobalSymmetryFeatures are NOT used: in a leave-one-space-group-out split the held-out group's
symmetry values never occur in training and trees extrapolate them arbitrarily (plan, Section 3 D).
Missing values stay NaN (LightGBM handles them); nothing is imputed.

The structure blocks run in parallel worker processes, in chunks that are saved to Drive as they
finish, so a Colab disconnect costs at most one chunk (resume by running the step again).
"""
from __future__ import annotations

import json
import os
import time
import warnings

import numpy as np
import pandas as pd

from .config import atomic_save_npz, atomic_write_json, code_hash, sha256_file
from .data6 import audit_decision, jarvis_structure, load_frozen, phase6_dirs, release_info, load_release

N_COMPOSITION, N_DENSITY, N_CRYSTALNN = 145, 3, 122


# ----------------------------------------------------------------------------- featurizers
def composition_featurizers():
    from matminer.featurizers.composition import ElementProperty, IonProperty, Stoichiometry, ValenceOrbital
    return [ElementProperty.from_preset("magpie", impute_nan=False), Stoichiometry(),
            ValenceOrbital(props=["frac"], impute_nan=False), IonProperty(fast=True, impute_nan=False)]


def composition_labels():
    out = []
    for f in composition_featurizers():
        out += [f"{type(f).__name__}|{n}" for n in f.feature_labels()]
    return out


def composition_block(formulas) -> np.ndarray:
    from pymatgen.core import Composition
    fz = composition_featurizers()
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for f in formulas:
            c = Composition(f)
            r = []
            for z in fz:
                try:
                    r += list(z.featurize(c))
                except Exception:
                    r += [np.nan] * len(z.feature_labels())
            rows.append(r)
    X = np.asarray(rows, dtype=float)
    if X.shape[1] != N_COMPOSITION:
        raise AssertionError(f"composition block has {X.shape[1]} columns, expected {N_COMPOSITION}")
    return X


def structure_labels():
    from matminer.featurizers.structure import DensityFeatures, SiteStatsFingerprint
    return ([f"DensityFeatures|{n}" for n in DensityFeatures().feature_labels()]
            + [f"CrystalNN_ops|{n}" for n in SiteStatsFingerprint.from_preset("CrystalNNFingerprint_ops").feature_labels()])


_W = {}


def _init_worker():
    warnings.simplefilter("ignore")
    from matminer.featurizers.structure import DensityFeatures, SiteStatsFingerprint
    _W["density"] = DensityFeatures()
    _W["cnn"] = SiteStatsFingerprint.from_preset("CrystalNNFingerprint_ops")


def featurize_structure(atoms: dict):
    """(row of 3 + 122 values, error message or ''). Called in a worker process."""
    if not _W:
        _init_worker()
    err = []
    try:
        st = jarvis_structure(atoms)
    except Exception as e:  # pragma: no cover - malformed record
        return [np.nan] * (N_DENSITY + N_CRYSTALNN), f"structure: {type(e).__name__}: {e}"
    try:
        d = list(_W["density"].featurize(st))
    except Exception as e:
        d = [float(st.density), float(st.volume / len(st)), np.nan]
        err.append(f"packing_fraction: {type(e).__name__}")
    try:
        c = list(_W["cnn"].featurize(st))
    except Exception as e:
        c = [np.nan] * N_CRYSTALNN
        err.append(f"crystalnn: {type(e).__name__}: {str(e)[:80]}")
    return d + c, "; ".join(err)


# ----------------------------------------------------------------------------- driver
def _workers(cfg) -> int:
    w = cfg["features"].get("workers", "auto")
    return max(1, (os.cpu_count() or 2) if w == "auto" else int(w))


def run_features6(cfg: dict):
    t0 = time.time()
    dirs = phase6_dirs(cfg)
    F = load_frozen(dirs)
    meta = audit_decision(dirs)
    info = release_info(dirs, meta)
    fdir = dirs["features"]
    chunk_dir = os.path.join(fdir, "structure_chunks")
    os.makedirs(chunk_dir, exist_ok=True)
    n = len(F)
    print(f"[features6] {n} entries; composition block (145, per reduced formula) ...", flush=True)

    # composition block: once per reduced formula, cached on Drive
    comp_file = os.path.join(fdir, "composition_block.npz")
    if os.path.exists(comp_file):
        with np.load(comp_file, allow_pickle=False) as z:
            forms, Xc_u = z["formulas"].astype(str), z["X"]
    else:
        forms = np.array(sorted(F.reduced_formula.unique()))
        t1 = time.time()
        Xc_u = composition_block(forms)
        atomic_save_npz(comp_file, formulas=forms, X=Xc_u)
        print(f"[features6]   {len(forms)} formulas in {(time.time() - t1) / 60:.1f} min", flush=True)
    pos = {f: i for i, f in enumerate(forms)}
    Xc = Xc_u[[pos[f] for f in F.reduced_formula]]

    # structure blocks: parallel, chunked, resumable
    cs = int(cfg["features"].get("chunk_size", 2000))
    chunks = [(i, F.jid.to_numpy()[i * cs:(i + 1) * cs]) for i in range((n + cs - 1) // cs)]
    todo = [(i, j) for i, j in chunks if not os.path.exists(os.path.join(chunk_dir, f"chunk_{i:04d}.npz"))]
    W = _workers(cfg)
    print(f"[features6] structure blocks (density 3 + CrystalNN 122): {len(chunks)} chunks of {cs}, "
          f"{len(chunks) - len(todo)} done earlier, {W} worker processes", flush=True)
    if todo:
        import multiprocessing as mp
        by_jid = {e["jid"]: e["atoms"] for e in load_release(info)}
        ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
        with ctx.Pool(W, initializer=_init_worker) as pool:
            for k, (i, jids) in enumerate(todo, 1):
                t1 = time.time()
                res = pool.map(featurize_structure, [by_jid[j] for j in jids], chunksize=max(1, len(jids) // (4 * W)))
                X = np.asarray([r[0] for r in res], dtype=float)
                errs = np.array([r[1] for r in res])
                atomic_save_npz(os.path.join(chunk_dir, f"chunk_{i:04d}.npz"), jid=jids.astype(str), X=X, err=errs)
                el = time.time() - t1
                print(f"[features6]   chunk {i + 1}/{len(chunks)}: {len(jids)} structures in {el:.0f} s "
                      f"({int((errs != '').sum())} with a warning); about {el * (len(todo) - k) / 60:.0f} min left",
                      flush=True)
        del by_jid

    Xs, errs = [], []
    for i, jids in chunks:
        with np.load(os.path.join(chunk_dir, f"chunk_{i:04d}.npz"), allow_pickle=False) as z:
            if not np.array_equal(z["jid"].astype(str), jids.astype(str)):
                raise AssertionError(f"chunk {i}: entries differ from the frozen table (delete the chunk and rerun)")
            Xs.append(z["X"])
            errs.append(z["err"].astype(str))
    Xs, errs = np.vstack(Xs), np.concatenate(errs)
    names = composition_labels() + structure_labels()
    X = np.hstack([Xc, Xs])
    if X.shape != (n, N_COMPOSITION + N_DENSITY + N_CRYSTALNN) or len(names) != X.shape[1]:
        raise AssertionError(f"feature matrix shape {X.shape} / {len(names)} names")
    block = np.array(["composition"] * N_COMPOSITION + ["density"] * N_DENSITY + ["crystalnn"] * N_CRYSTALNN)
    out = os.path.join(fdir, "features_v1.npz")
    atomic_save_npz(out, X=X, names=np.array(names), block=block, jid=F.jid.to_numpy().astype(str),
                    row=F.row.to_numpy())
    qa = feature_qa(X, names, block, F, errs)
    rec = dict(created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), code_hash=code_hash(),
               frozen_table_sha256=sha256_file(os.path.join(dirs["data"], "jarvis_frozen.csv.gz")),
               features_file=dict(path=out, sha256=sha256_file(out)), shape=list(X.shape),
               blocks=dict(composition=N_COMPOSITION, density=N_DENSITY, crystalnn=N_CRYSTALNN),
               matminer_version=_version("matminer"), pymatgen_version=_version("pymatgen"),
               minutes=(time.time() - t0) / 60, qa=qa)
    atomic_write_json(os.path.join(fdir, "features_v1.json"), rec)
    print(json.dumps(qa, indent=2, default=str))
    print(f"[features6] written {out} {X.shape} in {(time.time() - t0) / 60:.1f} min", flush=True)
    return rec


def _version(mod):
    try:
        import importlib.metadata as im
        return im.version(mod)
    except Exception:
        return "unknown"


def feature_qa(X, names, block, F, errs) -> dict:
    """Missing values, constant columns, warnings, and feature-identical entries (inputs the model
    cannot tell apart; with different targets they are irreducible noise, as on NOMAD)."""
    nan_col = np.isnan(X).mean(axis=0)
    finite = np.where(np.isfinite(X), X, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        const = np.nanmax(finite, axis=0) == np.nanmin(finite, axis=0)
    rows_nan = {b: int(np.isnan(X[:, block == b]).any(axis=1).sum()) for b in np.unique(block)}
    key = pd.util.hash_pandas_object(pd.DataFrame(np.round(np.nan_to_num(X, nan=-9.99e9), 10)), index=False)
    dup = key.duplicated(keep=False).to_numpy()
    dgroups = pd.Series(np.arange(len(X)))[dup].groupby(key[dup].to_numpy())
    spg = F.spg_number.to_numpy()
    n_cross_spg = sum(len(set(spg[g.to_numpy()])) > 1 for _, g in dgroups)
    return dict(n_columns_all_nan=int((nan_col == 1).sum()), columns_with_nan=int((nan_col > 0).sum()),
                max_nan_fraction=float(nan_col.max()),
                worst_columns={names[i]: round(float(nan_col[i]), 4) for i in np.argsort(-nan_col)[:5] if nan_col[i] > 0},
                n_constant_columns=int(const.sum()), constant_columns=[names[i] for i in np.where(const)[0]][:20],
                rows_with_nan_by_block=rows_nan,
                n_structure_warnings=int((errs != "").sum()),
                n_crystalnn_failed=int((np.char.find(errs.astype(str), "crystalnn") >= 0).sum()) if errs.size else 0,
                n_rows_in_feature_identical_sets=int(dup.sum()), n_feature_identical_sets=int(dgroups.ngroups),
                n_feature_identical_sets_spanning_space_groups=int(n_cross_spg))
