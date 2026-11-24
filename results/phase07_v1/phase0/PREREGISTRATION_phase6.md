# Phase 6 pre-registration (JARVIS-DFT, Tasks 6.4-6.5)

Written 2026-11-18T20:16:44Z by setup7, before any model is fitted on JARVIS.

* config hash baea14d0ddf1bcd464ec07c3a649521e4d33f4ece13e04a7b1758ffc145970fe
* code hash 43ee6e98c5adc2fe38ccba6cb20b23951a6b18623cd732c443dbb28e6181df20
* input content hashes: frozen_table 050448539876, labels 7b0cf364fe2a, features_X c3218ec310bc, features_jid 2ee9ee808e05, holdouts_v2 743fd8e3a66b

## Design

Regime B: MatFold criterion (space group re-detected with SpacegroupAnalyzer, symprec 0.1),
nested K = 10 outer folds (two halves x 5) and L = 5 inner folds by
space group; groups with fewer than 200 rows in the target's scope always train. Methods: CV+, group-weighted
CV+, normalized CV+ (kNN sigma), CQR-CV+ (alpha 0.10), ensemble-SD baseline; contrasts: random
inner folds, polymorph exclusion, repeat 2, LOCSO. Regime A: 10 formula-grouped splits with
CV+ random, split, Mondrian, normalized split, CQR. H3: 10 partitions (pool = one half),
200 calibration/test reshuffles each, HCP and single subsampling.

* formation: 81884 rows, 59 eligible groups (75197 rows); halves [30, 29] groups; groups per outer fold [6, 6, 6, 6, 6, 6, 5, 6, 6, 6]; design hash 2326bfffc75d
* bandgap: 22475 rows, 31 eligible groups (15902 rows); halves [15, 16] groups; groups per outer fold [3, 3, 3, 3, 3, 3, 3, 4, 3, 3]; design hash 87e03d80be0a
* tbmbj: 7829 rows, 10 eligible groups (3083 rows); halves [5, 5] groups; groups per outer fold [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; design hash 15d3bdd4aa84

## Decision rules

Decision rules (alpha = 0.10; each target judged separately; formation energy and the non-metal
OptB88vdW gap are primary, TBmBJ is descriptive only). The primary coverage of a held-out space
group is its CV+ coverage (inner folds by space group, absolute residual), averaged over the two
repeats.

* G5: at least 75% of eligible groups have primary coverage >= 0.85: "most held-out
  groups are covered" (the paper then shows which fail and whether distance predicts them).
* H1-J (Regime A): supported if the marginal CV+ random coverage is in [0.88, 0.92]
  AND at least 10% of eligible groups have CV+ random coverage < 0.85 AND Mondrian
  gives >= 0.88 in at least 90% of groups; heterogeneity absent if
  fewer than 5% of groups are below 0.85.
* G3-J: among groups with repeat-1 CV+ coverage < 0.8, normalized CV+ or CQR-CV+ lifts at least
  50% to >= 0.85: "adaptive scores help"; otherwise negative, as on NOMAD.
* H3-J: the valid group-level design (HCP, Lee et al. eq. 6; single subsampling, Dunn et al. 2023)
  is supported if the mean group-weighted coverage is >= 0.88 for both, rejected if < 0.85
  for either.
* H4-J (evaluated in v8): best distance AUROC - width AUROC >= 0.05 with 95 % cluster-bootstrap CI
  above 0 AND best distance AUROC >= 0.6; G4-J if every score has pooled and
  within-group AUROC < 0.7.
* H5-J (evaluated in v8): in the 3 worst-covered groups, mean coverage at m = 20 (one structure per
  composition) >= 0.88 in all 3 AND dispersion ratio vs the permutation null <= 1.25
  in at least 2 of 3.
* H2 is descriptive only (calibration-tail share by group; ablation of the 10 % of training groups
  with the largest out-of-fold residuals).

## Deviations

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

## Environment and time

* {"python": "3.13.18", "platform": "Linux-6.6.159+-x86_64-with-glibc2.39", "cpus": 8, "packages": {"numpy": "2.4.6", "scipy": "1.16.3", "pandas": "2.2.3", "scikit-learn": "1.6.1", "lightgbm": "4.6.0", "optuna": "5.1.0", "pymatgen": "2026.10.24", "matminer": "0.10.1", "MatFold": "1.3.3"}}
* benchmark {"formation_fit_seconds_1000_rounds": 20.993663787841797, "bandgap_fit_seconds_1000_rounds": 7.416641473770142, "parallel_fits": 2, "threads_per_fit": 4}
* projected hours (about 1,000 rounds per fit) {"formation": {"primary": 0.44, "randinner": 0.15, "polyex": 0.15, "rep2": 0.15, "locso": 0.12, "regimeA": 0.23, "h3": 0.03, "tuning": 0.24, "total": 1.5}, "bandgap": {"primary": 0.15, "randinner": 0.05, "polyex": 0.05, "rep2": 0.05, "locso": 0.04, "regimeA": 0.08, "h3": 0.01, "tuning": 0.08, "total": 0.5}, "tbmbj": {"total": 0.05}}
