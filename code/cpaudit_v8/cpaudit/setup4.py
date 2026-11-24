"""Setup for Phase 4: environment, inputs from Phase 2, and the immutable Phase 4 pre-registration
(written before any Phase 4 score is computed)."""
from __future__ import annotations

import glob
import json
import os
import platform
import shutil
import sys
import time

from .config import atomic_write_json, atomic_write_text, code_hash, config_hash, get_paths, sha256_file


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def environment() -> dict:
    import importlib
    pk = {}
    for m in ("numpy", "scipy", "pandas", "sklearn", "yaml", "matplotlib"):
        try:
            pk[m] = getattr(importlib.import_module(m), "__version__", None)
        except Exception:
            pk[m] = None
    return dict(time_utc=_now(), python=sys.version.split()[0], platform=platform.platform(),
                cpus=os.cpu_count(), packages=pk)


def link_phase2(cfg) -> dict:
    """Use exactly the Phase 2 dataset: copy its de-duplication record and check its units."""
    paths = get_paths(cfg)
    src = os.path.join(paths.phase2, "phase0", "dedup_v2.json")
    dst = os.path.join(paths.phase0, "dedup_v2.json")
    if not os.path.exists(src):
        raise FileNotFoundError(f"Phase 2 de-duplication record not found: {src}")
    if not os.path.exists(dst):
        shutil.copyfile(src, dst)
    if sha256_file(src) != sha256_file(dst):
        raise AssertionError("phase0/dedup_v2.json differs from the Phase 2 record")
    units = sorted(os.path.basename(f) for f in glob.glob(os.path.join(paths.phase2_loso, "*.npz")))
    need = [f"{t}_{d}{g}.npz" for t in ("bandgap", "formation")
            for d, gs in (("loso", (227, 12, 33, 167, 194, 206)), ("l2o", (12, 33, 167, 194, 206))) for g in gs]
    missing = [u for u in need if u not in units]
    if missing:
        raise FileNotFoundError(f"Phase 2 units missing: {missing}")
    info = dict(time_utc=_now(), dedup_sha256=sha256_file(dst),
                units={u: sha256_file(os.path.join(paths.phase2_loso, u)) for u in need})
    atomic_write_json(os.path.join(paths.phase0, "phase2_inputs.json"), info)
    atomic_write_json(os.path.join(paths.phase0, "environment.json"), environment())
    print(f"  Phase 2 inputs linked: {len(need)} units, de-duplication record {info['dedup_sha256'][:12]}")
    return info


PREREG_TEXT = """# Pre-registration: Phase 4 (detecting failure before it matters)

Registered (UTC): {time}
Config hash: {chash}
Code hash:   {khash}
Phase 2 inputs: phase0/phase2_inputs.json (SHA-256 of all 22 units)

Written once, before any Phase 4 score is computed; never overwritten. A later change to the
scientific settings is refused unless explicitly allowed and is then logged in DEVIATIONS.json.

## What was known when this was written
Phase 1, 2 and 3 results, including: LOSO error is mostly a lattice-wide offset; a domain
classifier separates every lattice from the rest (AUROC 1.000); for LOSO-Fd-3m band gap the
median kNN distance of the held-out structures (97-109) lies inside the range of the out-of-fold
lattices (85-131). No Phase 4 score, AUROC or flag rate had been computed.

## Data and configurations
n = 2,343 (Phase 2 de-duplication). 6 nested-LOSO holdouts and 5 leave-two-out holdouts, using
the Phase 2 models without refitting. Primary analysis: the 6 LOSO holdouts pooled (every
structure appears once).

## Outcomes (per structure, per target)
* Primary: miscovered by the reference interval, {ref} at alpha = {a} (the submitted design).
* Secondary: |error| above the in-distribution {q}-quantile (random-fold out-of-fold |residuals|
  on the same training lattices).

## Scores (higher = more doubt)
width ({ref} interval width); knn (mean distance to the {k} nearest training structures,
standardized SOAP, PCA-{pca}); kde (Schultz et al. 2025: Epanechnikov KDE on all standardized
features, bandwidth by the nearest-neighbour rule with quantile {bq}, d = 1 - KDE/max_train KDE);
committee (SD of the lattice-fold models); subsample (SD of the random-fold models; LOSO only);
classifier (out-of-fold probability of belonging to the screened set; logistic, PCA-{cpca}).

## Analyses
1. Pooled structure-level AUROC and average precision, 95 % cluster-bootstrap CI ({B} resamples of
   composition groups).
2. Within-lattice AUROC for every holdout with >= {mc} structures in each class, and their
   n-weighted mean.
3. Selective prediction: interval coverage and MAE of the retained structures when abstaining
   on the highest-score fraction (grid {grid0}-{grid1}).
4. Lattice level (11 configurations, descriptive): Spearman correlation of the mean score with the
   miscoverage rate and with |bias|, permutation p-value ({perm} permutations). The 11
   configurations reuse 6 test lattices, so they are not independent.
5. Flag: kNN-distance conformal p-values (Bates et al. 2023), flag at p <= {fa}; report the
   in-distribution false-flag rate and the flagged fraction of every held-out lattice.

## Decision rules (primary outcome, each target separately)
* H4 supported if max(AUROC_knn, AUROC_kde) - AUROC_width >= {dmin}, with the 95 % cluster-
  bootstrap CI of that difference (max taken inside each resample) above 0, and
  max(AUROC_knn, AUROC_kde) >= {amin}. Otherwise H4 is not supported.
* G4 triggered if every score has pooled AUROC < {g4} AND every score's n-weighted within-lattice
  AUROC is < {g4}. Then the recommendation narrows to: obtain labelled calibration data for every
  family screened, and Phase 5 carries the practical message.
* Expected for the flag (stated now): false-flag rate near {fa} and close to 100 % of every held-out
  lattice flagged, i.e. the flag detects "new family", not "will fail". The paper reports it as such.
* Every score, outcome, lattice and configuration is reported whichever way it falls.
"""


def prereg_paths(paths):
    return (os.path.join(paths.phase0, "PREREGISTRATION_phase4.json"),
            os.path.join(paths.phase0, "PREREGISTRATION_phase4.md"))


def preregister(cfg) -> str:
    paths = get_paths(cfg)
    js, md = prereg_paths(paths)
    if os.path.exists(js):
        pre = json.load(open(js))
        print(f"  pre-registration already exists ({pre['time_utc']}); config unchanged: "
              f"{pre['config_hash'] == config_hash(cfg)}")
        return md
    if os.path.exists(os.path.join(paths.scores, "phase4_structures.csv")) or glob.glob(
            os.path.join(paths.scores, "features_*.npz")):
        raise RuntimeError("Phase 4 scores already exist; a pre-registration written now would not be one.")
    link_phase2(cfg)
    d, r = cfg["detect"], cfg["decisions"]
    info = dict(time_utc=_now(), config_hash=config_hash(cfg), code_hash=code_hash(),
                config_path=cfg["_config_path"], config={k: v for k, v in cfg.items() if not k.startswith("_")})
    text = PREREG_TEXT.format(
        time=info["time_utc"], chash=info["config_hash"], khash=info["code_hash"], ref=cfg["reference_interval"],
        a=cfg["alpha_primary"], q=d["large_error_quantile"], k=d["knn"]["k"], pca=d["knn"]["pca_components"],
        bq=d["kde"]["bandwidth_quantile"], cpca=d["classifier"]["pca_components"], B=d["bootstrap"],
        mc=d["min_class"], grid0=d["abstention_grid"][0], grid1=d["abstention_grid"][-1],
        perm=d["lattice_permutations"], fa=d["flag"]["alpha"], dmin=r["h4_min_diff"], amin=r["h4_min_auroc"],
        g4=r["g4_auroc"])
    atomic_write_text(md, text)
    atomic_write_json(js, info)
    print(f"  wrote {md}")
    return md


def check_preregistration(cfg, paths, allow_deviation: bool):
    js, _ = prereg_paths(paths)
    if not os.path.exists(js):
        raise RuntimeError("No Phase 4 pre-registration found. Run the `setup` step first.")
    pre = json.load(open(js))
    now = config_hash(cfg)
    if pre["config_hash"] != now:
        msg = f"Config changed since pre-registration ({pre['config_hash'][:12]} -> {now[:12]})."
        if not allow_deviation:
            raise RuntimeError(msg + " Re-run with --allow-deviation to proceed; the change will be logged.")
        log = os.path.join(paths.phase0, "DEVIATIONS.json")
        items = json.load(open(log)) if os.path.exists(log) else []
        items.append({"time_utc": _now(), "from": pre["config_hash"], "to": now, "config": cfg["_config_path"]})
        atomic_write_json(log, items)
        print("  WARNING: " + msg + " Logged in phase0/DEVIATIONS.json")
