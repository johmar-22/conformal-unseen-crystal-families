"""Setup for Phase 5: environment, links to the Phase 2 inputs, and the immutable Phase 5
pre-registration (written before any Phase 5 quantity is computed on the real data)."""
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
    for m in ("numpy", "scipy", "pandas", "sklearn", "yaml", "matplotlib", "lightgbm"):
        try:
            pk[m] = getattr(importlib.import_module(m), "__version__", None)
        except Exception:
            pk[m] = None
    return dict(time_utc=_now(), python=sys.version.split()[0], platform=platform.platform(),
                cpus=os.cpu_count(), packages=pk)


def needed_inputs(cfg) -> list[str]:
    files = []
    for t in ("bandgap", "formation"):
        for c in cfg["configs"]:
            d, g = c.split(":")
            files.append(f"{t}_{d}{g}.npz")
        for c in cfg["train"]["configs"]:
            d, g = c.split(":")
            files.append(f"{t}_loso{g}.hp.json")       # leave-two-out refits reused the LOSO-X hyperparameters
    return sorted(set(files))


def link_phase2(cfg) -> dict:
    paths = get_paths(cfg)
    src = os.path.join(paths.phase2, "phase0", "dedup_v2.json")
    dst = os.path.join(paths.phase0, "dedup_v2.json")
    if not os.path.exists(src):
        raise FileNotFoundError(f"Phase 2 de-duplication record not found: {src}")
    if not os.path.exists(dst):
        shutil.copyfile(src, dst)
    if sha256_file(src) != sha256_file(dst):
        raise AssertionError("phase0/dedup_v2.json differs from the Phase 2 record")
    need = needed_inputs(cfg)
    missing = [u for u in need if not os.path.exists(os.path.join(paths.phase2_loso, u))]
    if missing:
        raise FileNotFoundError(f"Phase 2 files missing in {paths.phase2_loso}: {missing}")
    info = dict(time_utc=_now(), dedup_sha256=sha256_file(dst),
                files={u: sha256_file(os.path.join(paths.phase2_loso, u)) for u in need})
    atomic_write_json(os.path.join(paths.phase0, "phase2_inputs.json"), info)
    atomic_write_json(os.path.join(paths.phase0, "environment.json"), environment())
    print(f"  Phase 2 inputs linked: {len(need)} files, de-duplication record {info['dedup_sha256'][:12]}")
    return info


PREREG_TEXT = """# Pre-registration: Phase 5 (few-shot, family-specific recalibration)

Registered (UTC): {time}
Config hash: {chash}
Code hash:   {khash}
Phase 2 inputs: phase0/phase2_inputs.json (SHA-256 of every file used)

Written once, before any Phase 5 quantity is computed on the real data; never overwritten. A
later change to the scientific settings is refused unless explicitly allowed, and is then
logged in phase0/DEVIATIONS.json.

## What was known when this was written
All Phase 1-4 results. In particular: under LOSO the error of a held-out lattice is mostly one
offset for the whole lattice (|bias|/MAE >= 0.84 in 9 of the 12 LOSO lattice-target pairs; band gap
P63/mmc 0.49 and Pna21 0.00, formation P63/mmc 0.66; leave-two-out P63/mmc 0.97 and 0.86);
the Fd-3m cell-setting sub-families have offsets +0.625/+0.673 eV (band gap) and -0.099/-0.147
eV/atom (formation); no score detects failure (Phase 4, G4 triggered); composition groups, not
structures, are the exchangeability unit in-distribution (post1). No few-shot coverage, width or
dispersion had been computed on the real data; the code was tested on synthetic data only.

## Question
How many labelled structures of a new family restore valid, useful intervals, and does removing
the family offset make them narrower?

## Families and model
{nconf} families: {confs}. For each, the fixed Phase 2 model mu = mean of the five lattice-fold
models (the predictor whose CV+ intervals were evaluated in Phases 2 and 4), trained without the
family; nothing is refitted. Residual z = y - mu. Targets: band gap and formation energy.
alpha = {a}.

## Methods (each draw uses the same m labelled structures)
* split: split conformal with |z| scores, k = ceil((1 - alpha)(m + 1)).
* offset_full: full conformal for y = mu + b + e with b the mean residual of the augmented data
  (exact marginal coverage k/(m + 1) under exchangeability; closed form in cpaudit/fewshot.py).
  For this model the full-conformal set is always one interval and equals the jackknife+
  interval with the leave-one-out mean offset (Barber et al. 2021); the identity is verified in
  the self-test, so the exact full-conformal guarantee applies (not only jackknife+'s 1 - 2 alpha).

## Sampling schemes ({D} draws per scheme and m)
* random: m random structures, m in {mr}; coverage of all other structures of the family.
  Reference: exact finite-population law (negative hypergeometric; mean k/(m + 1), tends to
  Beta(k, m + 1 - k)).
* clustered: whole composition groups in random order until m structures, m in {mc} (as when
  several configurations of one composition are computed).
* one_per_group: m composition groups, one structure each (single subsampling, Dunn et al. 2023),
  m in {mo}.
  For clustered and one_per_group, at least {mtg} composition groups must remain untouched.
  Primary coverage = cov_new: mean over the untouched groups of the covered fraction (new
  compositions, the screening target). Also recorded: cov_family, the covered fraction of all
  non-calibration structures.
* Permutation null for clustered and one_per_group: the same procedure after permuting the
  residuals across the family's structures (same group sizes, no group structure).
* sub_a2b / sub_b2a (Fd-3m only, exploratory): calibrate on one cell setting (227a, gamma ~30;
  227b, gamma ~120), test on the other; m in {ms}.

## Why the plan's H5 is re-specified
With random sampling inside a family the mean coverage is fixed by construction: for random
labels the family coverage of split and offset_full is exactly k/(m + 1), and for one_per_group
labels so is the coverage of new compositions (Dunn et al. 2023). Fd-3m heterogeneity cannot move
these means. What composition structure can change is the spread across draws, i.e. the chance
that a particular set of labels gives a badly calibrated interval (calibration-conditional
coverage; Vovk 2012; Angelopoulos & Bates 2023, Sec. 3.2), and, for clustered labels, the mean as
well. H5 is therefore tested on the spread, against a permutation null with the same group sizes
and no group structure.

## Decision rules (each target separately)
D5.0 Code check (random scheme, split and offset_full): in every family and m,
     |mean coverage - k/(m + 1)| <= {nse} Monte Carlo SE (split: exact-law SD / sqrt(draws);
     offset_full: empirical SD / sqrt(draws)). A cell outside this band is re-run once with
     10 x the draws and a new seed; if it is still outside, that is a code error: stop and do
     not interpret Phase 5.
H5   (primary; plan H5 for the recommended protocol) one_per_group, split, m = {h5m}:
     rho = SD(cov_new) / SD(cov_new under the permutation null).
     Supported if rho <= {rs} in >= {smin} of {nconf} families; rejected if rho > {rr} in >= {rmin}
     families; otherwise partial. Rejection means the labels needed depend on which compositions
     are labelled, beyond exchangeable theory.
H5c  (clustered labels) clustered, split, m = {h5m}: a family counts as harmed if its rho > {rs} or
     its mean cov_new < k/(m + 1) - {cdrop}. Supported (clustering harms) if >= {cmin} families are
     harmed; otherwise not supported.
H5b  (offset correction) Cells fixed in advance from Phase 2 E1 (|bias|/MAE >= 0.80):
     band gap {cbg}; formation {cfe}. Statistic at one_per_group, m = {bm}: width ratio =
     median width (offset_full) / median width (split), and mean cov_new of offset_full.
     Supported if ratio <= {wr} and mean coverage >= {cmin_b} in >= {share:.0%} of the cells;
     rejected if ratio > {wrej} in at least half of the cells; otherwise partial.
Practical endpoint (reported for every family, scheme and method; stated expectation, not a
     rule): m* = smallest m with mean coverage >= {pmean}, P(coverage < {floor}) <= {prisk} and no
     unbounded interval. Expected under exchangeable theory: m* = 30 (Beta P(coverage < 0.80) =
     0.069 at m = 20 and 0.044 at m = 30), later if the test set is small.
Exploratory, with expectations stated now: (i) clustered calibration is more dispersed than
     one_per_group at equal m; (ii) sub-family transfer for Fd-3m: band gap near nominal (offsets
     +0.625 vs +0.673 eV), formation 227a -> 227b below nominal (offsets -0.099 vs -0.147 eV/atom).
Also reported: interval width in units of the family's SD next to the LOSO CV+ width (Phase 2),
     and MAE before and after subtracting the estimated offset.

## Optional training arm (`train5`, descriptive, no decision rule)
Families {tconfs}. Per draw, {tf:.0%} of the composition groups are held out as test; from the
other groups, n_cal = {ncal} structures are drawn for calibration and m in {tm} more for training.
Arms: refit with the m structures added, calibrated on n_cal (train); the fixed Phase 2 model
calibrated on m + n_cal (calibrate); the fixed Phase 2 model calibrated on n_cal (fixed); a refit
without the m structures calibrated on n_cal (baseline; separates the refit from the labels).
{tdraws} draws per m, Phase 2 hyperparameters, split and offset_full. Question: are new DFT
labels better spent on the model or on calibration?

Every family, target, scheme, m and method is reported whichever way it falls.
"""


def prereg_paths(paths):
    return (os.path.join(paths.phase0, "PREREGISTRATION_phase5.json"),
            os.path.join(paths.phase0, "PREREGISTRATION_phase5.md"))


def _cells(cfg, target):
    from .data import LATTICE_NAME
    out = []
    for c in cfg["decisions"]["h5b_cells"][target]:
        d, g = c.split(":")
        out.append(f"{LATTICE_NAME[int(g)]} ({'LOSO' if d == 'loso' else 'refit'})")
    return ", ".join(out)


def preregister(cfg) -> str:
    paths = get_paths(cfg)
    js, md = prereg_paths(paths)
    if os.path.exists(js):
        pre = json.load(open(js))
        print(f"  pre-registration already exists ({pre['time_utc']}); config unchanged: "
              f"{pre['config_hash'] == config_hash(cfg)}")
        return md
    if glob.glob(os.path.join(paths.out, "draws", "*.csv.gz")) or glob.glob(os.path.join(paths.out, "train", "*")):
        raise RuntimeError("Phase 5 results already exist; a pre-registration written now would not be one.")
    link_phase2(cfg)
    from .data import LATTICE_NAME
    f, r, t = cfg["fewshot"], cfg["decisions"], cfg["train"]

    def lab(c):
        d, g = c.split(":")
        return f"{LATTICE_NAME[int(g)]} ({'LOSO' if d == 'loso' else 'refit without Fd-3m'})"
    info = dict(time_utc=_now(), config_hash=config_hash(cfg), code_hash=code_hash(),
                config_path=cfg["_config_path"], config={k: v for k, v in cfg.items() if not k.startswith("_")})
    text = PREREG_TEXT.format(
        time=info["time_utc"], chash=info["config_hash"], khash=info["code_hash"],
        nconf=len(cfg["configs"]), confs="; ".join(lab(c) for c in cfg["configs"]), a=cfg["alpha_primary"],
        D=f["draws"], mr=f["m_random"], mc=f["m_clustered"], mo=f["m_one_per_group"], ms=f["m_subfamily"],
        mtg=f["min_test_groups"], nse=r["sanity_n_se"], h5m=r["h5_m"], rs=r["h5_ratio_support"],
        rr=r["h5_ratio_reject"], smin=r["h5_support_min"], rmin=r["h5_reject_min"],
        cbg=_cells(cfg, "bandgap"), cfe=_cells(cfg, "formation"), bm=r["h5b_m"], wr=r["h5b_width_ratio"],
        cmin_b=r["h5b_min_coverage"], share=r["h5b_share"], wrej=r["h5b_width_reject"],
        cdrop=r["h5c_mean_drop"], cmin=r["h5c_min_families"],
        pmean=r["practical_mean"], floor=r["coverage_floor"], prisk=r["practical_risk"],
        tconfs="; ".join(lab(c) for c in t["configs"]), tf=t["test_group_fraction"], ncal=t["n_cal"],
        tm=t["m"], tdraws=t["draws"])
    atomic_write_text(md, text)
    atomic_write_json(js, info)
    print(f"  wrote {md}")
    return md


def check_preregistration(cfg, paths, allow_deviation: bool):
    js, _ = prereg_paths(paths)
    if not os.path.exists(js):
        raise RuntimeError("No Phase 5 pre-registration found. Run the `setup` step first.")
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
