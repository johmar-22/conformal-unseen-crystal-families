# cpaudit_v9: final paper figures (Fig. 1-8, S1-S7) from the results on Drive; earlier steps unchanged

## v9 (2026-11-19): paper figures

Upload `cpaudit_v9.zip` and `cpaudit_run_v9.py` to `MyDrive/NOMAD_Uncer_Quan/code_bundles/` (nothing is
replaced). One figure at a time (new time-stamped folder under `paper_figures_v5_drafts/` each run):

```python
from google.colab import drive; drive.mount("/content/drive", force_remount=True)
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v9.py" figures9 --only fig1
```

The full set (`figures9` without `--only`) goes to `paper_figures_v5/`, refused if that folder exists and is
not empty. Nothing is fitted or recomputed: every number comes from a result file of Phases 1-8, and every
file read is listed with its SHA-256 in `FIGURES9_MANIFEST.json` (the run stops if any file outside the Drive
folder was read). All figures pass the automatic overlap check (`figcheck.py`) before they are saved.

| Figure | Content | Drawn from |
|---|---|---|
| Fig. 1 | study design: split designs, family sizes, families available vs size threshold, workflow (every count in a box read from a result file; methods listed only if results exist) | phase01, phase02 reports; phase06 audit, data, features; phase07 units; phase08 report |
| Fig. 2 | regime A on NOMAD (a-d) and group-level Mondrian (e-g) | phase01 report, phase02 post1 |
| Fig. 3 | regime B on NOMAD: coverage of the held-out lattice by method | phase02 report |
| Fig. 4 | mechanism (a-d) and no covariate overlap for weighted CP (e, f) | phase02 loso, overlap |
| Fig. 5 | failure detection, NOMAD (a-c) and JARVIS-DFT (d-f), never pooled | phase04 report; phase08 report, posthoc |
| Fig. 6 | NOMAD few-shot recalibration (a-f) and the dispersion test (g, h) | phase05 report, posthoc |
| Fig. 7 | JARVIS-DFT unseen space groups (a-d) and group-level guarantees (e, f) | phase07 unitsB, unitsH3; phase08 |
| Fig. 8 | JARVIS-DFT few-shot recalibration of the 3 worst groups per target | phase08 fewshot, posthoc |
| S1-S7 | calibration tail; formation transfer; detection AUROC matrix; lattice-level detection; sub-family transfer and training arm; JARVIS few-shot per group; hierarchical jackknife+ vs CV+ | as above |

Composite figures reuse the drawing code of the single figures: each part is drawn into a matplotlib
sub-figure, its panel letters are re-assigned (also the letters inside legends, e.g. "(a, b)"), and the whole
page goes through the overlap check once.

---

# Earlier steps (v8.1 and before)

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

## v8.1 (2026-11-19): post hoc analyses and figure fixes (NOT pre-registered)

Designed after the v8 results were seen. Reads `phase07_v1/` and `phase08_v1/`; writes only
`phase08_v1/posthoc/` and `paper_figures_v4/`. A rerun recomputes and stops if any table would change.

```python
from google.colab import drive; drive.mount("/content/drive", force_remount=True)
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v81.py" posthoc8
```

* `posthoc8`
  1. Few-shot vs no model: the fewshot8 label sets are re-created with the registered seeds and
     checked draw by draw against the saved draws; on each, a no-model interval is built from the
     same labels (offset-corrected full conformal on y, the posthoc5 definition). Also the coverage
     spread (5 % / 95 %) over draws.
  2. Inner-model SD at group level (Spearman with group bootstrap and permutation p; CV+ width and
     kNN distance for comparison) and abstention inside the three worst groups.
  3. Oxides split into ABO3 compositions (composition rule only) and other oxides.
  4. Degenerate group-level intervals: atom 1/(K+1) vs alpha (HJ+, group-weighted CV+, HCP) and the
     single-subsampling rank.
* `figures81`: Fig. 8 with leader lines on every label; Fig. 9 legend (0.70 = G4-J threshold,
  0.50 = chance, 0.80 = HJ+ bound); Fig. 10 redesigned (coverage band over draws; width relative to
  the no-model interval). Refuses to draw into a folder that already holds a complete set.

The v8 steps setup8 / scores8 / fewshot8 / hjplus8 / report8 stop with "code changed" under v8.1:
the registered v8 results are final.
