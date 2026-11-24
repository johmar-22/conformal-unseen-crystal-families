"""setup8: checks that v7 is complete and unchanged, and writes the Task 6.5 pre-registration
addendum (phase08_v1/phase0/PREREGISTRATION_phase6_v8.md/.json) before any v8 result exists.

What is new here and what is not
--------------------------------
* H4-J and H5-J (and G4-J) were pre-registered in v7 (2026-11-18 20:16 UTC). setup8 stops if
  their rules (models7.decisions) differ from the v7 pre-registration.
* The addendum fixes the operational details that v7 left open (listed in ADDENDUM_TEXT), and
  records the three worst-covered groups per target chosen by the v7 rule. It is written AFTER
  the v7 results (G5, H1-J, G3-J, H3-J) were seen and BEFORE any detection, few-shot or
  hierarchical jackknife+ result is computed; the paper must say so.
* scores8, fewshot8, hjplus8 and report8 refuse to run if the config, the code or a v7 input
  changed since the addendum, unless --allow-deviation is given (then logged in deviations.json).
"""
from __future__ import annotations

import hashlib
import json
import os
import time

from .config import atomic_write_json, atomic_write_text, code_hash

ADDENDUM_TEXT = """
Operational details fixed before any Task 6.5 result is computed (v7 results already seen):

* Structure-level analyses (H4-J, H5-J, H2, hierarchical jackknife+) use the repeat-1 units of
  Regime B; group-level summaries (coverage, width, MAE, interval score, distance) average the two
  repeats where repeat 2 exists (TBmBJ: repeat 1 only).
* Distance to training = NOMAD Phase 4 kNN score: mean Euclidean distance to the {distance[k]} nearest
  outer-training structures (rare groups + all calibration groups of the outer fold), standardized
  features with training-median imputation, first {distance[pca_components]} principal components. Group distance =
  median over the group's held-out structures.
* Across groups: shares of groups below 0.90 and 0.80 (primary CV+, alpha 0.10); Spearman rho of
  group coverage with group distance; 95 % percentile bootstrap CI over groups ({bootstrap} resamples)
  and a Monte Carlo permutation p-value ({n_perm} permutations).
* H4-J (rule pre-registered in v7, unchanged): outcome = held-out structure outside its primary
  CV+ interval (alpha 0.10, repeat 1). Scores: CV+ width, kNN distance, KDE dissimilarity
  (Schultz et al. 2025; Epanechnikov, bandwidth rule quantile {kde[bandwidth_quantile]}; reference set and bandwidth
  sample drawn at random when the training set exceeds {kde[n_reference]} / {kde[n_bandwidth]} rows), committee
  (SD of the 5 inner models) and an out-of-fold domain classifier (PCA-{classifier[pca_components]}, logistic, {classifier[cv_folds]} folds
  grouped by reduced formula). "Cluster bootstrap over groups" = resampling held-out space groups
  for the pooled AUROCs and their difference; within-group AUROCs (groups with >= {min_class} structures
  in each class), n-weighted, resample reduced formulas within each group. Abstention curves and
  the conformal kNN flag (p <= {flag[alpha]}, Bates et al. 2023; training rows split 3:1:1 into fit,
  calibration and in-distribution test parts, formula-grouped) are descriptive. The three worst
  groups' within-group AUROCs are reported separately (plan: "failure detector for the three
  worst-covered groups").
* H5-J (rule pre-registered in v7, unchanged): the 3 groups with the lowest primary coverage
  (mean over repeats; ties: lower repeat-1 coverage). Fixed model = mean of the 5 repeat-1 inner
  models (never saw the group); composition = reduced formula; m in {fewshot[m]}; {fewshot[draws]} draws per
  scheme and m; schemes random, clustered and one structure per composition, with permutation
  null; methods split and offset-corrected full conformal (NOMAD Phase 5 code). The rule uses
  the split method and the one-per-composition scheme at m = 20. Stated before the run: the
  JARVIS table keeps one entry per (reduced formula, space group), so within a held-out group
  almost every composition occurs once (at most 1.08 entries per formula in any eligible group).
  The one-per-composition, clustered and random schemes therefore nearly coincide, the dispersion
  ratio against the permutation null is close to 1 by construction, and H5-J rests in practice
  on the mean-coverage criterion. NOMAD's clustering contrast (H5c) cannot be replicated here.
* H3 (valid design run in v7): additionally reported as how often a held-out group is covered at
  the nominal level (within-group coverage >= 1 - alpha) over all partitions and reshuffles,
  and the share of unbounded thresholds.
* Hierarchical jackknife+ (plan: "cheaper alternative using all groups"; Lee, Barber & Willett,
  arXiv:2306.06342 v4, Appendix C.2.2, eq. 20, Theorem 5: coverage >= 1 - 2 alpha under
  hierarchical exchangeability and a hierarchically symmetric algorithm). For every repeat-1
  outer fold, one LightGBM model per training space group k, fitted on all other training
  groups and the rare groups (the half's tuned hyperparameters, one seed for all models, rows in
  a fixed order). Every calibration group therefore trains the models of the others and is
  scored by a model that never saw it. Targets: {hjplus[targets]}. No decision rule: reported
  descriptively against the theorem's bound and against the v7 CV+ and group-weighted CV+
  intervals on the same held-out structures. TBmBJ has 10 groups (fewer than about 20), so only
  this version is reported for it; with 9 calibration groups the weight 1/(K+1) = 0.10 equals
  alpha = 0.10 (upper end = largest term) and alpha = 0.05 is unbounded.
* H2 (descriptive): calibration-tail enrichment per group (share of residuals above the
  {h2[tail_quantile]}-quantile relative to the group's share of calibration points) and a calibration-only
  ablation that removes the {h2[drop_share]:.0%} of calibration groups with the largest out-of-fold MAE and
  recomputes CV+ with the same models (as NOMAD's no-Fd-3m calibration ablation).
* Oxides (descriptive, no extra fits): held-out structures classed as oxides by the
  oxidation-state rule of the audit; groups with >= {oxide[min_rows]} such structures.
* JARVIS is never pooled with NOMAD: H6 is a side-by-side table of separate results.
"""


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def files8(d8):
    return (os.path.join(d8["phase0"], "PREREGISTRATION_phase6_v8.json"),
            os.path.join(d8["phase0"], "PREREGISTRATION_phase6_v8.md"))


def config_hash8(cfg) -> str:
    """Scientific settings that v8 depends on: models7 (incl. the pre-registered rules), models8, lgbm."""
    from .fit7 import m7
    from .run8 import m8
    keep = dict(models7=m7(cfg), models8=m8(cfg), lgbm=cfg["lgbm"])
    return hashlib.sha256(json.dumps(keep, sort_keys=True).encode()).hexdigest()


def _v7_checks(cfg):
    """v7 pre-registration present, its rules unchanged, inputs unchanged, units complete."""
    from .fit7 import content_hashes, dirs7, m7
    from .setup7 import prereg_files
    d7 = dirs7(cfg)
    pj, _ = prereg_files(d7)
    if not os.path.exists(pj):
        raise RuntimeError("v7 pre-registration not found: v7 must be complete before v8.")
    pre7 = json.load(open(pj))
    now7 = json.loads(json.dumps(m7(cfg)))
    if now7 != pre7["models7"]:
        diff = sorted(k for k in set(now7) | set(pre7["models7"]) if now7.get(k) != pre7["models7"].get(k))
        raise RuntimeError(f"models7 differs from the v7 pre-registration in {diff}: the H4-J / H5-J rules "
                           "must be the registered ones.")
    inp = content_hashes(d7)
    bad = [k for k, v in pre7["inputs"].items() if inp.get(k) != v]
    if bad:
        raise RuntimeError(f"inputs changed since the v7 pre-registration: {bad}")
    m = m7(cfg)
    n_out = 2 * m["outer_folds_per_half"]
    missing = []
    for t in m["targets"]:
        for blk, rep in (("primary", 1), ("rep2", 2)):
            if blk not in m["blocks"][t]:
                continue
            for o in range(n_out):
                s = os.path.join(d7["unitsB"], f"{t}_{blk}_r{rep}_o{o:02d}")
                if not (os.path.exists(s + ".npz") and os.path.exists(s + ".csv")):
                    missing.append(os.path.basename(s))
        if "h3" in m["blocks"][t]:
            for p in range(int(m["h3"]["n_partitions"])):
                if not os.path.exists(os.path.join(d7["unitsH3"], f"{t}_H3_p{p:02d}.csv")):
                    missing.append(f"{t}_H3_p{p:02d}")
        for which in ("h0", "h1"):
            if not os.path.exists(os.path.join(d7["tune"], f"hp_{t}_{which}.json")):
                missing.append(f"hp_{t}_{which}")
    if missing:
        raise RuntimeError(f"v7 outputs missing ({len(missing)}): {missing[:8]} ...")
    return pre7, inp


def run_setup8(cfg):
    from .run8 import DETECT_TARGETS, MOCK, Ctx8, dirs8, m8, worst_groups
    d8 = dirs8(cfg)
    pj, pm = files8(d8)
    if os.path.exists(pj):
        pre = json.load(open(pj))
        print(f"[setup8] addendum already written at {pre['time_utc']} (not rewritten)")
        require_prereg8(cfg, allow_deviation=False)
        print("[setup8] config, code and v7 inputs unchanged since the addendum")
        return pre
    print("[setup8] checking v7 (pre-registration, rules, inputs, units) ...", flush=True)
    pre7, inp = _v7_checks(cfg)
    ctx = Ctx8(cfg)
    m = m8(cfg)
    worst = {}
    for t in DETECT_TARGETS:
        W = worst_groups(ctx, t, int(m["fewshot"]["n_worst"]))
        worst[t] = [dict(group=int(r.group), symbol=ctx.sgsym.get(int(r.group), str(r.group)),
                         primary_coverage=round(float(r.coverage), 4)) for _, r in W.iterrows()]
    pre = dict(time_utc=_now(), mock=MOCK, config_hash8=config_hash8(cfg), code_hash=code_hash(),
               v7_preregistration=dict(time_utc=pre7["time_utc"], config_hash=pre7["config_hash"],
                                       code_hash=pre7["code_hash"]),
               registered_rules_v7=dict(h4=pre7["models7"]["decisions"]["h4"], h5=pre7["models7"]["decisions"]["h5"]),
               v7_inputs=inp, models8=m, worst_groups=worst)
    text = ADDENDUM_TEXT.format(**m)
    md = [f"# Phase 6 Task 6.5 pre-registration addendum (v8){' - MOCK RUN' if MOCK else ''}", "",
          f"Written {pre['time_utc']} by setup8, after the v7 results and before any v8 result.", "",
          f"* v7 pre-registration {pre7['time_utc']} (config hash {pre7['config_hash'][:12]}); its rules "
          "are checked unchanged, including H4-J and H5-J:",
          f"  - H4-J {json.dumps(pre['registered_rules_v7']['h4'])}",
          f"  - H5-J {json.dumps(pre['registered_rules_v7']['h5'])}",
          f"* v8 config hash {pre['config_hash8'][:12]}, code hash {pre['code_hash'][:12]}",
          "* v7 input content hashes unchanged: " + ", ".join(f"{k} {v[:12]}" for k, v in inp.items()), "",
          "## Worst-covered groups (H5-J; v7 rule: lowest primary coverage, mean of 2 repeats)", ""]
    for t, L in worst.items():
        md.append(f"* {t}: " + "; ".join(f"{w['symbol']} (#{w['group']}, {w['primary_coverage']:.3f})" for w in L))
    md += ["", "## Operational details", text]
    atomic_write_json(pj, pre)
    atomic_write_text(pm, "\n".join(md))
    print("\n".join(md))
    print(f"[setup8] addendum written: {pj}")
    return pre


def require_prereg8(cfg, allow_deviation: bool):
    from .fit7 import content_hashes, dirs7
    from .run8 import dirs8
    d8 = dirs8(cfg)
    pj, _ = files8(d8)
    if not os.path.exists(pj):
        raise RuntimeError("No v8 addendum found: run setup8 first.")
    pre = json.load(open(pj))
    problems = []
    if pre["config_hash8"] != config_hash8(cfg):
        problems.append(f"config changed ({pre['config_hash8'][:12]} -> {config_hash8(cfg)[:12]})")
    if pre["code_hash"] != code_hash():
        problems.append(f"code changed ({pre['code_hash'][:12]} -> {code_hash()[:12]})")
    now = content_hashes(dirs7(cfg))
    problems += [f"v7 input {k} changed" for k, v in pre["v7_inputs"].items() if now.get(k) != v]
    if problems:
        msg = "; ".join(problems)
        if not allow_deviation:
            raise RuntimeError(f"Since the v8 addendum: {msg}. Re-run with --allow-deviation to proceed; "
                               "the change is then logged in deviations.json.")
        f = os.path.join(d8["phase0"], "deviations.json")
        items = json.load(open(f)) if os.path.exists(f) else []
        items.append(dict(time_utc=_now(), problems=problems, config_hash8=config_hash8(cfg), code_hash=code_hash()))
        atomic_write_json(f, items)
        print(f"[prereg8] DEVIATION logged: {msg}", flush=True)
    return pre
