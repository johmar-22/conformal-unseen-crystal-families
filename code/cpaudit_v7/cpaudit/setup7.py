"""setup7: design, input hashes, timing benchmark and the Phase 6 pre-registration (v7).

The pre-registration is written once, before any model is fitted, and records:
  * config hash (all scientific settings, including the decision rules) and code hash;
  * content hashes of every input (frozen table, labels, features, holdouts);
  * the design of every target (content hash and per-fold summary);
  * the decision rules in words.
tune7 / fit7 refuse to run if the pre-registration is missing or if the config, the code or an
input changed since, unless --allow-deviation is given (then the change is logged in
deviations.json next to the pre-registration).
"""
from __future__ import annotations

import json
import os
import platform
import time

import numpy as np

from .config import atomic_write_json, atomic_write_text, code_hash, config_hash

RULES_TEXT = """
Decision rules (alpha = 0.10; each target judged separately; formation energy and the non-metal
OptB88vdW gap are primary, TBmBJ is descriptive only). The primary coverage of a held-out space
group is its CV+ coverage (inner folds by space group, absolute residual), averaged over the two
repeats.

* G5: at least {g5[share]:.0%} of eligible groups have primary coverage >= {g5[coverage]}: "most held-out
  groups are covered" (the paper then shows which fail and whether distance predicts them).
* H1-J (Regime A): supported if the marginal CV+ random coverage is in [{h1[marginal_lo]}, {h1[marginal_hi]}]
  AND at least {h1[share_support]:.0%} of eligible groups have CV+ random coverage < {h1[low]} AND Mondrian
  gives >= {h1[mondrian_coverage]} in at least {h1[mondrian_share]:.0%} of groups; heterogeneity absent if
  fewer than {h1[share_absent]:.0%} of groups are below {h1[low]}.
* G3-J: among groups with repeat-1 CV+ coverage < {g3[failing]}, normalized CV+ or CQR-CV+ lifts at least
  {g3[share]:.0%} to >= {g3[lifted]}: "adaptive scores help"; otherwise negative, as on NOMAD.
* H3-J: the valid group-level design (HCP, Lee et al. eq. 6; single subsampling, Dunn et al. 2023)
  is supported if the mean group-weighted coverage is >= {h3[support]} for both, rejected if < {h3[reject]}
  for either.
* H4-J (evaluated in v8): best distance AUROC - width AUROC >= {h4[delta]} with 95 % cluster-bootstrap CI
  above 0 AND best distance AUROC >= {h4[min_distance_auroc]}; G4-J if every score has pooled and
  within-group AUROC < {h4[g4_threshold]}.
* H5-J (evaluated in v8): in the 3 worst-covered groups, mean coverage at m = {h5[m]} (one structure per
  composition) >= {h5[mean_coverage]} in all 3 AND dispersion ratio vs the permutation null <= {h5[dispersion_ratio]}
  in at least {h5[min_groups]} of 3.
* H2 is descriptive only (calibration-tail share by group; ablation of the 10 % of training groups
  with the largest out-of-fold residuals).
"""

DEVIATIONS_TEXT = """
Deviations from the revision plan, stated in advance:
* Outer design: 10 held-out folds of space groups (MatFold nested K = 10, L = 5) instead of
  leave-one-space-group-out; the MatFold paper itself uses K = 10 for space groups.
* Inner folds: 5 folds by space group instead of leave-one-family-out.
* Tuning: cross-half nested tuning (h0 is tuned without half 0's groups and used for half 0's
  folds; h1 likewise) instead of tuning per holdout. LOCSO reuses h0, whose tuning data overlap
  the held-out crystal system (secondary analysis; disclosed).
* Band gap restricted to non-metals (gap > 0): 72.6 % of entries have a gap of 0.
* TBmBJ descriptive only (10 eligible groups among non-metals).
* No ALIGNN (optional in the plan).
* Group-weighted CV+ is our extension of Lee, Barber & Willett (arXiv:2306.06342 v4), who prove
  the weighting for split conformal (Theorem 1) and leave-one-group-out jackknife+ (Theorem 5).
"""


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def environment() -> dict:
    import importlib.metadata as im
    from .models import n_cpu, plan_parallel
    v = {}
    for p in ("numpy", "scipy", "pandas", "scikit-learn", "lightgbm", "optuna", "pymatgen", "matminer", "MatFold"):
        try:
            v[p] = im.version(p)
        except Exception:
            v[p] = "missing"
    return dict(python=platform.python_version(), platform=platform.platform(), cpus=n_cpu(), packages=v)


def prereg_files(d):
    return (os.path.join(d["phase0"], "PREREGISTRATION_phase6.json"),
            os.path.join(d["phase0"], "PREREGISTRATION_phase6.md"))


def benchmark(cfg, ctx) -> dict:
    """Time one 1,000-round fit at the size of a formation and of a band-gap unit, on a synthetic
    target (no real result is computed before the pre-registration)."""
    from .fit7 import _lgb_task, _mock_target, _run
    from .models import plan_parallel
    out = {}
    hp = dict(learning_rate=0.05, num_leaves=63, max_depth=8, min_child_samples=30, subsample=0.8,
              colsample_bytree=0.8, reg_alpha=0.01, reg_lambda=0.01, n_rounds=1000)
    rng = np.random.default_rng(0)
    for name, n in (("formation", 65000), ("bandgap", 18000)):
        rows = np.sort(rng.choice(ctx.X.shape[0], size=n, replace=False))
        y = _mock_target(np.asarray(ctx.X[rows]), "formation")
        t0 = time.time()
        _run([_lgb_task(ctx, rows, y, {"p": rows[:2000]}, hp, 1)], cfg)
        out[f"{name}_fit_seconds_1000_rounds"] = time.time() - t0
    w, th = plan_parallel(cfg)
    out.update(parallel_fits=w, threads_per_fit=th)
    return out


def projected_hours(m, bench, hp_rounds=1000) -> dict:
    """Rough projection, assuming about 1,000 boosting rounds per fit (tuning sets the real number)."""
    w = bench["parallel_fits"]
    tf, tb = bench["formation_fit_seconds_1000_rounds"], bench["bandgap_fit_seconds_1000_rounds"]
    n_out = 2 * m["outer_folds_per_half"]
    L = m["inner_folds"]
    per = {"primary": n_out * L * 3, "randinner": n_out * L, "polyex": n_out * L, "rep2": n_out * L,
           "locso": 7 * 6, "regimeA": 10 * 8, "h3": m["h3"]["n_partitions"],
           "tuning": 3 * m["tuning"]["n_trials"] * m["tuning"]["inner_folds"] * 0.6}
    out = {}
    for t, tt in (("formation", tf), ("bandgap", tb)):
        out[t] = {b: round(n * tt / w / 3600, 2) for b, n in per.items()}
        out[t]["total"] = round(sum(out[t].values()), 1)
    out["tbmbj"] = {"total": round((n_out * L + 2 * 45 * 0.6) * tb * 0.45 / w / 3600, 2)}
    return out


def run_setup7(cfg):
    from .fit7 import MOCK, Ctx, build_designs, content_hashes
    ctx = Ctx(cfg)
    d = ctx.d
    pj, pm = prereg_files(d)
    if os.path.exists(pj):
        with open(pj) as f:
            pre = json.load(f)
        print(f"[setup7] pre-registration already written at {pre['time_utc']} (not rewritten)")
        require_preregistration(cfg, allow_deviation=False)
        print("[setup7] config, code and inputs unchanged since the pre-registration")
        return pre
    print("[setup7] building and checking the designs ...", flush=True)
    designs = build_designs(cfg, ctx)
    print("[setup7] hashing inputs ...", flush=True)
    hashes = content_hashes(d)
    print("[setup7] timing one 1,000-round fit per data size (synthetic target) ...", flush=True)
    bench = benchmark(cfg, ctx)
    proj = projected_hours(ctx.m, bench)
    m = ctx.m
    pre = dict(time_utc=_now(), mock=MOCK, config_hash=config_hash(cfg), code_hash=code_hash(),
               config_path=cfg["_config_path"], models7=m, lgbm=cfg["lgbm"], inputs=hashes,
               designs=designs, environment=environment(), benchmark=bench, projected_hours=proj)
    rules = RULES_TEXT.format(**m["decisions"])
    md = [f"# Phase 6 pre-registration (JARVIS-DFT, Tasks 6.4-6.5){' - MOCK RUN' if MOCK else ''}", "",
          f"Written {pre['time_utc']} by setup7, before any model is fitted on JARVIS.", "",
          f"* config hash {pre['config_hash']}", f"* code hash {pre['code_hash']}",
          "* input content hashes: " + ", ".join(f"{k} {v[:12]}" for k, v in hashes.items()), "",
          "## Design", "",
          "Regime B: MatFold criterion (space group re-detected with SpacegroupAnalyzer, symprec 0.1),",
          f"nested K = {2 * m['outer_folds_per_half']} outer folds (two halves x {m['outer_folds_per_half']}) and L = {m['inner_folds']} inner folds by",
          "space group; groups with fewer than "
          f"{m['min_group_size']} rows in the target's scope always train. Methods: CV+, group-weighted",
          "CV+, normalized CV+ (kNN sigma), CQR-CV+ (alpha 0.10), ensemble-SD baseline; contrasts: random",
          "inner folds, polymorph exclusion, repeat 2, LOCSO. Regime A: 10 formula-grouped splits with",
          "CV+ random, split, Mondrian, normalized split, CQR. H3: 10 partitions (pool = one half),",
          f"{m['h3']['n_reshuffles']} calibration/test reshuffles each, HCP and single subsampling.", ""]
    for t, s in designs.items():
        md.append(f"* {t}: {s['n_rows']} rows, {s['n_eligible']} eligible groups ({s['rows_in_eligible']} rows); "
                  f"halves {s['half_sizes_groups']} groups; groups per outer fold {s['outer_groups_per_fold_r1']}; "
                  f"design hash {s['content_hash'][:12]}")
    md += ["", "## Decision rules", rules, "## Deviations", DEVIATIONS_TEXT, "## Environment and time", "",
           f"* {json.dumps(pre['environment'])}", f"* benchmark {json.dumps(bench)}",
           f"* projected hours (about 1,000 rounds per fit) {json.dumps(proj)}", ""]
    atomic_write_json(pj, pre)
    atomic_write_text(pm, "\n".join(md))
    print("\n".join(md))
    print(f"[setup7] pre-registration written: {pj}")
    return pre


def require_preregistration(cfg, allow_deviation: bool):
    from .fit7 import content_hashes, dirs7
    d = dirs7(cfg)
    pj, _ = prereg_files(d)
    if not os.path.exists(pj):
        raise RuntimeError("No Phase 6 pre-registration found: run setup7 first.")
    with open(pj) as f:
        pre = json.load(f)
    problems = []
    if pre["config_hash"] != config_hash(cfg):
        problems.append(f"config changed ({pre['config_hash'][:12]} -> {config_hash(cfg)[:12]})")
    if pre["code_hash"] != code_hash():
        problems.append(f"code changed ({pre['code_hash'][:12]} -> {code_hash()[:12]})")
    now = content_hashes(d)
    for k, v in pre["inputs"].items():
        if now.get(k) != v:
            problems.append(f"input {k} changed")
    if problems:
        msg = "; ".join(problems)
        if not allow_deviation:
            raise RuntimeError(f"Since the pre-registration: {msg}. Re-run with --allow-deviation to "
                               "proceed; the change is then logged in deviations.json.")
        f = os.path.join(d["phase0"], "deviations.json")
        items = json.load(open(f)) if os.path.exists(f) else []
        items.append(dict(time_utc=_now(), problems=problems, config_hash=config_hash(cfg), code_hash=code_hash()))
        atomic_write_json(f, items)
        print(f"[prereg] DEVIATION logged: {msg}", flush=True)
    return pre
