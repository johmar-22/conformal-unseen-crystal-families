# Phase 6 Task 6.5 pre-registration addendum (v8)

Written 2026-11-19T06:44:03Z by setup8, after the v7 results and before any v8 result.

* v7 pre-registration 2026-11-18T20:16:44Z (config hash baea14d0ddf1); its rules are checked unchanged, including H4-J and H5-J:
  - H4-J {"delta": 0.05, "min_distance_auroc": 0.6, "g4_threshold": 0.7}
  - H5-J {"m": 20, "mean_coverage": 0.88, "dispersion_ratio": 1.25, "min_groups": 2}
* v8 config hash 452f4cb9e312, code hash f260c5e34f8d
* v7 input content hashes unchanged: frozen_table 050448539876, labels 7b0cf364fe2a, features_X c3218ec310bc, features_jid 2ee9ee808e05, holdouts_v2 743fd8e3a66b

## Worst-covered groups (H5-J; v7 rule: lowest primary coverage, mean of 2 repeats)

* formation: F-43m (#216, 0.582); Pm (#6, 0.704); P-6m2 (#187, 0.771)
* bandgap: P4/nmm (#129, 0.784); P4mm (#99, 0.788); Fm-3m (#225, 0.805)

## Operational details

Operational details fixed before any Task 6.5 result is computed (v7 results already seen):

* Structure-level analyses (H4-J, H5-J, H2, hierarchical jackknife+) use the repeat-1 units of
  Regime B; group-level summaries (coverage, width, MAE, interval score, distance) average the two
  repeats where repeat 2 exists (TBmBJ: repeat 1 only).
* Distance to training = NOMAD Phase 4 kNN score: mean Euclidean distance to the 25 nearest
  outer-training structures (rare groups + all calibration groups of the outer fold), standardized
  features with training-median imputation, first 50 principal components. Group distance =
  median over the group's held-out structures.
* Across groups: shares of groups below 0.90 and 0.80 (primary CV+, alpha 0.10); Spearman rho of
  group coverage with group distance; 95 % percentile bootstrap CI over groups (2000 resamples)
  and a Monte Carlo permutation p-value (100000 permutations).
* H4-J (rule pre-registered in v7, unchanged): outcome = held-out structure outside its primary
  CV+ interval (alpha 0.10, repeat 1). Scores: CV+ width, kNN distance, KDE dissimilarity
  (Schultz et al. 2025; Epanechnikov, bandwidth rule quantile 0.3; reference set and bandwidth
  sample drawn at random when the training set exceeds 20000 / 5000 rows), committee
  (SD of the 5 inner models) and an out-of-fold domain classifier (PCA-50, logistic, 5 folds
  grouped by reduced formula). "Cluster bootstrap over groups" = resampling held-out space groups
  for the pooled AUROCs and their difference; within-group AUROCs (groups with >= 10 structures
  in each class), n-weighted, resample reduced formulas within each group. Abstention curves and
  the conformal kNN flag (p <= 0.05, Bates et al. 2023; training rows split 3:1:1 into fit,
  calibration and in-distribution test parts, formula-grouped) are descriptive. The three worst
  groups' within-group AUROCs are reported separately (plan: "failure detector for the three
  worst-covered groups").
* H5-J (rule pre-registered in v7, unchanged): the 3 groups with the lowest primary coverage
  (mean over repeats; ties: lower repeat-1 coverage). Fixed model = mean of the 5 repeat-1 inner
  models (never saw the group); composition = reduced formula; m in [10, 20, 30]; 500 draws per
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
  scored by a model that never saw it. Targets: ['bandgap', 'tbmbj', 'formation']. No decision rule: reported
  descriptively against the theorem's bound and against the v7 CV+ and group-weighted CV+
  intervals on the same held-out structures. TBmBJ has 10 groups (fewer than about 20), so only
  this version is reported for it; with 9 calibration groups the weight 1/(K+1) = 0.10 equals
  alpha = 0.10 (upper end = largest term) and alpha = 0.05 is unbounded.
* H2 (descriptive): calibration-tail enrichment per group (share of residuals above the
  0.9-quantile relative to the group's share of calibration points) and a calibration-only
  ablation that removes the 10% of calibration groups with the largest out-of-fold MAE and
  recomputes CV+ with the same models (as NOMAD's no-Fd-3m calibration ablation).
* Oxides (descriptive, no extra fits): held-out structures classed as oxides by the
  oxidation-state rule of the audit; groups with >= 50 such structures.
* JARVIS is never pooled with NOMAD: H6 is a side-by-side table of separate results.
