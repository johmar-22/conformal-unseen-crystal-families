"""Phase 5: few-shot, family-specific recalibration (no model is refitted).

For each held-out family (6 LOSO lattices and the leave-two-out P63/mmc refit) the fixed Phase 2
model mu (mean of the five lattice-fold models, i.e. the predictor whose CV+ intervals were
evaluated in Phases 2 and 4) is recalibrated with m labelled structures of that family, and the
interval is evaluated on the family's other structures.

Sampling schemes (all random, repeated `draws` times per m):
  random         m random structures; coverage on all other structures of the family
                 (exact finite-population reference law known).
  clustered      whole composition groups in random order until m structures (as when several
                 configurations of one composition are computed).
  one_per_group  m groups, one structure each (single subsampling, Dunn et al. 2023).
  sub_a2b/b2a    Fd-3m only (exploratory): calibrate on one cell setting, test on the other.
For clustered and one_per_group two coverages are recorded: cov_new = mean over the untouched
composition groups of the covered fraction (new compositions; primary, and exactly k/(m+1) in
expectation for one_per_group), and cov_family = covered fraction of all non-calibration
structures. Permutation null (clustered, one_per_group): the same procedure after permuting the
residuals across the family's structures (same group sizes, no group structure).
"""
from __future__ import annotations

import json
import os
import time
import zlib

import numpy as np
import pandas as pd

from .config import atomic_write_bytes, atomic_write_csv, get_paths
from .data import LATTICE_NAME, id_table
from .fewshot import (METHODS, SET_FUNCS, Family, draw_clustered, draw_one_per_group, draw_random,
                      draw_subfamily, group_weighted, hull_width, in_set)
from .setup5 import check_preregistration

TARGETS = ("bandgap", "formation")


def parse_config(s: str):
    d, g = s.split(":")
    return d, int(g)


def config_label(design: str, g: int) -> str:
    return f"{LATTICE_NAME[g]} ({'LOSO' if design == 'loso' else 'refit without Fd-3m'})"


def load_family(paths, ids_tab, design: str, g: int, target: str) -> Family:
    f = os.path.join(paths.phase2_loso, f"{target}_{design}{g}.npz")
    z = np.load(f, allow_pickle=False)
    mu = z["mu_te_l"].mean(axis=0) if design == "loso" else z["mu_te"].mean(axis=0)
    ids = z["test_id"].astype(int)
    tab = ids_tab.loc[ids]
    if not (tab["sg"].to_numpy() == g).all():
        raise AssertionError(f"{f}: test set contains structures outside lattice {g}")
    y = z["y_test"].astype(float)
    return Family(z=y - mu, groups=tab["group"].to_numpy(), sub=tab["sub"].to_numpy(), ids=ids, y=y)


def _seed(base: int, *parts) -> int:
    return (int(base) + zlib.crc32("|".join(map(str, parts)).encode())) % (2 ** 32 - 1)


def _evaluate(z_cal, z_new, z_family, alpha, new_fn):
    """Coverage of the new compositions (primary for clustered / one_per_group: group-weighted over
    untouched groups) and of the whole remaining family (point-weighted)."""
    rows = []
    b = float(np.mean(z_cal))
    for meth in METHODS:
        S = SET_FUNCS[meth](z_cal, alpha)
        w = hull_width(S)
        cov_fam = float(in_set(z_family, S).mean())
        cov_new = new_fn(in_set(z_new, S)) if z_new is not None else np.nan
        rows.append(dict(method=meth, cov=(cov_new if z_new is not None else cov_fam), cov_new=cov_new,
                         cov_family=cov_fam, width=w, infinite=bool(not np.isfinite(w)), n_comp=len(S), b_hat=b,
                         mae_raw=float(np.mean(np.abs(z_family))),
                         mae_offset=float(np.mean(np.abs(z_family - b)))))
    return rows


def run_family(fam: Family, cfg: dict, design: str, g: int, target: str) -> pd.DataFrame:
    f = cfg["fewshot"]
    alpha = float(cfg["alpha_primary"])
    D, mtg = int(f["draws"]), int(f["min_test_groups"])
    out = []
    schemes = [("random", False), ("clustered", False), ("clustered", True),
               ("one_per_group", False), ("one_per_group", True)]
    if design == "loso" and g == 227 and f.get("subfamily", True):
        schemes += [("sub_a2b", False), ("sub_b2a", False)]
    allidx = np.arange(fam.N)
    for scheme, null in schemes:
        grid = f[{"random": "m_random", "clustered": "m_clustered", "one_per_group": "m_one_per_group"}.get(
            scheme, "m_subfamily")]
        for m in grid:
            rng = np.random.default_rng(_seed(cfg["seed_base"], design, g, target, scheme, m, null))
            for d in range(D):
                z = fam.z[rng.permutation(fam.N)] if null else fam.z
                new, tg = None, None
                if scheme == "random":
                    cal, rest = draw_random(fam, m, rng)
                elif scheme in ("clustered", "one_per_group"):
                    r = (draw_clustered if scheme == "clustered" else draw_one_per_group)(fam, m, rng, mtg)
                    if r is None:
                        break
                    cal, new, tg = r
                    rest = np.setdiff1d(allidx, cal, assume_unique=True)
                else:
                    src, dst = ("227a", "227b") if scheme == "sub_a2b" else ("227b", "227a")
                    r = draw_subfamily(fam, m, rng, src, dst)
                    if r is None:
                        break
                    cal, rest = r
                new_fn = ((lambda c, _t=new, _g=tg: group_weighted(c, _t, fam, _g)) if new is not None else None)
                rows = _evaluate(z[cal], (z[new] if new is not None else None), z[rest], alpha, new_fn)
                for row in rows:
                    out.append(dict(design=design, heldout=LATTICE_NAME[g], sg=g, target=target, scheme=scheme,
                                    null=null, m=int(m), draw=d, n_test=len(rest),
                                    n_test_groups=(len(tg) if tg is not None else np.nan), **row))
    return pd.DataFrame(out)


def run_phase5(cfg: dict, allow_deviation: bool = False):
    paths = get_paths(cfg)
    check_preregistration(cfg, paths, allow_deviation)
    ids_tab = id_table(cfg, paths)
    draws_dir = os.path.join(paths.out, "draws")
    os.makedirs(draws_dir, exist_ok=True)
    meta = []
    t0 = time.time()
    confs = [parse_config(c) for c in cfg["configs"]]
    print(f"Phase 5 | {len(confs)} families x 2 targets | {cfg['fewshot']['draws']} draws per (scheme, m)",
          flush=True)
    for i, (design, g) in enumerate(confs, 1):
        for target in TARGETS:
            fam = load_family(paths, ids_tab, design, g, target)
            meta.append(dict(design=design, heldout=LATTICE_NAME[g], sg=g, target=target, N=fam.N, G=fam.G,
                             sd_y=float(np.std(fam.y, ddof=1)), bias_full=float(fam.z.mean()),
                             mae_full=float(np.mean(np.abs(fam.z))),
                             n_227a=int(np.sum(fam.sub == "227a")), n_227b=int(np.sum(fam.sub == "227b"))))
            fn = os.path.join(draws_dir, f"{target}_{design}{g}.csv.gz")
            if os.path.exists(fn):
                print(f"  [{i}/{len(confs)}] {target} {config_label(design, g)}: done earlier, skipped", flush=True)
                continue
            t1 = time.time()
            df = run_family(fam, cfg, design, g, target)
            atomic_write_bytes(fn, _gz(df))
            print(f"  [{i}/{len(confs)}] {target} {config_label(design, g)}: N={fam.N}, groups={fam.G}, "
                  f"{len(df)} rows ({time.time() - t1:.0f} s)", flush=True)
    atomic_write_csv(os.path.join(paths.out, "phase5_families.csv"), pd.DataFrame(meta))
    with open(os.path.join(paths.out, "phase5_run.json"), "w") as fh:
        json.dump(dict(minutes=(time.time() - t0) / 60), fh)
    print(f"Phase 5 draws written ({(time.time() - t0) / 60:.1f} min). Next: report5")


def _gz(df: pd.DataFrame) -> bytes:
    import gzip
    import io
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
        gz.write(df.to_csv(index=False).encode())
    return buf.getvalue()
