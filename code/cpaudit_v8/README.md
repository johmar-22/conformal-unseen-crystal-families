# cpaudit_v8: Phase 6 on JARVIS-DFT (Tasks 6.2-6.3 data and splits; 6.4 models and UQ; 6.5 analyses)

Upload `cpaudit_v7.zip` and `cpaudit_run_v7.py` to `MyDrive/NOMAD_Uncer_Quan/code_bundles/`
(next to v1-v6.1; nothing is replaced). One short Colab cell per step:

```python
from google.colab import drive; drive.mount("/content/drive", force_remount=True)
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v7.py" setup
```

| Step | What it does | Time (8 CPUs) |
|---|---|---|
| `setup` | self-test (all earlier tools plus the v6 data, feature and split tools) | 3-5 min |
| `freeze6` | analysis table rebuilt from the release on Drive; stops unless it equals the Task 6.1 audit | ~1 min |
| `features6` | Magpie 145 + density 3 + CrystalNN 122 features; parallel; checkpoint every 2,000 structures | 15-25 min |
| `splits6` | LOSGO (59 space groups) and LOCSO (7 crystal systems) under MatFold labels (primary) and JARVIS labels (sensitivity), polymorph-exclusion variant | 5-8 min |
| `matfold6` | MatFold itself on all 81,884 structures; its labels and leave-one-out splits must equal ours row by row | ~10 min |
| `report6` | one summary | seconds |
| `all` | the five steps above | |

No target is used and no model is fitted. Outputs: `phase06_v1/data/`, `features/`, `splits/`,
`splits_matfold/`, `report/`. Earlier outputs are only read.

## Design

* **Freeze.** The pinned release on Drive (SHA-256 checked against the audit) is processed with the
  audit's own functions; the step stops unless the row count, the 59 eligible space groups and their
  sizes equal the audit. Row order = JARVIS id, used by every later step. Formulas such as "NaN"
  (sodium nitride) are kept as text when the table is read back, and floats round-trip exactly.
* **Features** (matminer; Ward et al. 2018): the full Magpie set of Ward et al. (2016)
  (ElementProperty 132, Stoichiometry 6, ValenceOrbital 4, IonProperty 3), DensityFeatures (3) and
  SiteStatsFingerprint `CrystalNNFingerprint_ops` (122; Zimmermann & Jain 2020). No
  GlobalSymmetryFeatures (their values for a held-out space group never occur in training).
  NaN is kept, not imputed. QA: missing values, constant columns, featurizer warnings and
  feature-identical entries.
* **Splits** (MatFold protocol; Witman & Schindler 2025, checked against the MatFold 1.3.3 source).
  Labels are re-detected from each structure with SpacegroupAnalyzer(symprec = 0.1), exactly as
  MatFold does; these are the primary labels. The JARVIS database labels are a sensitivity family
  (2.2 % of entries change space group, mostly P1, Pm, P2/m and Cm; per-group test-set Jaccard median
  0.97, minimum 0.78). Leave-one-label-out outer folds; space groups with fewer than 200 entries are
  never held out and always train (MatFold `fraction_lower_limit`); 59 groups under either label
  system. Leave-one-crystal-system-out over all 7 systems. The polymorph-exclusion variant (training
  drops every entry whose reduced formula occurs in the held-out group) is our extension, not part of
  MatFold. `matfold6` runs MatFold itself on every structure, keeps its CSV and JSON split files
  (recreatable with `MatFold.from_json`) and stops unless its labels, compositions and every fold
  equal ours.

## Task 6.4 (v7): models and conformal intervals, pre-registered

Steps: `setup7` (design, input content hashes, timing, PRE-REGISTRATION), `tune7`, `fit7`,
`report7`. Outputs in `phase07_v1/`. Full design: `PHASE6_v7_DESIGN.md` (project); the
pre-registration written by `setup7` is the binding version.

* Targets: formation energy (81,884 rows, 59 eligible groups), OptB88vdW gap of non-metals (gap > 0;
  22,475 rows, 31 groups), TBmBJ gap of non-metals (7,829 rows, 10 groups; descriptive only).
* Regime B: MatFold nested (K = 10, L = 5): two halves of the eligible groups, 5 outer folds each,
  inner folds by space group; rare groups always train. Methods: CV+, group-weighted CV+ (Lee et al.
  weighting; our CV+ extension), normalized CV+ (kNN sigma), CQR-CV+ (alpha 0.10), ensemble SD.
  Contrasts: random inner folds, polymorph exclusion, repeat 2, LOCSO.
* Regime A: 10 formula-grouped splits; CV+ random, split, Mondrian per space group, normalized
  split, CQR.
* H3: 10 partitions (pool = one half), 200 reshuffles, HCP (Lee et al. eq. 6) and single
  subsampling (Dunn et al. 2023).
* Tuning: cross-half nested (h0 without half 0's groups, used for half 0's folds; h1 likewise),
  Regime A on split 0's proper set; Optuna TPE, 15 trials, early stopping, max 2,000 rounds.
* Fast exact CV+ (`cvfast.py`): bisection on the value with exact term counts; equal to the
  brute-force definition bit for bit in the self-tests and spot-checked inside every unit.
* Guards: tune7 / fit7 stop if the pre-registration is missing or the config, code or inputs
  changed (content hashes, independent of the zlib build); `--allow-deviation` logs the change.


## Task 6.5 (v8): analyses, after a pre-registration addendum

v8 reads phase07_v1 (v7 units, designs, tuned hyperparameters; never modified) and writes phase08_v1/
and paper_figures_v3/. Order: setup8, scores8, fewshot8, hjplus8, report8, figures8.

* setup8: checks that the v7 pre-registration, its rules (including H4-J and H5-J) and its inputs are
  unchanged and that all v7 units exist, then writes PREREGISTRATION_phase6_v8.md/.json with the
  operational details and the three worst-covered groups per target (v7 rule), before any v8 result.
* scores8: kNN distance to training (NOMAD Phase 4 score), KDE dissimilarity (Schultz et al. 2025),
  inner-model SD, domain classifier, conformal kNN flag, for every held-out structure (repeat 1);
  kNN distance only for repeat 2 and TBmBJ.
* fewshot8: H5-J with the NOMAD Phase 5 functions on the 3 worst groups per target.
* hjplus8: hierarchical jackknife+ of Lee, Barber & Willett (arXiv:2306.06342 v4, App. C.2.2, eq. 20,
  Theorem 5: >= 1 - 2 alpha), one leave-one-group-out model per training space group per outer fold,
  one seed and a fixed row order for all models (symmetric algorithm); checked against eq. 20 written
  literally inside every unit.
* report8: per-group table (coverage, width, MAE, interval score, distance), shares below 0.90 / 0.80,
  Spearman with a bootstrap over groups, H3 coverage at the nominal level, HJ+ against CV+ and
  group-weighted CV+, H4-J / G4-J, H5-J, H2 (descriptive), oxide subset, H6 side by side (never pooled).
* Verification: self-tests (t_phase8), a full mock run on the real features with synthetic targets,
  independent recomputation of every reported quantity from the saved files, guard tests.
