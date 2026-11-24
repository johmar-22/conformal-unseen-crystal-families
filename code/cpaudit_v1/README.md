# cpaudit_v1: Phase 0 and Phase 1 of the revision plan

Conformal-prediction reliability audit for the NOMAD 2018 TCO benchmark. This bundle
implements Phase 0 (freeze, integrity checks, pre-registration) and Phase 1 (a clean
in-distribution baseline over 25 repeated, composition-grouped splits).

## How to run (Colab Pro, A100 runtime)

1. Create `MyDrive/NOMAD_Uncer_Quan/code_bundles/` and upload **both** `cpaudit_v1.zip` and
   `cpaudit_run_v1.py` into it.
2. One short cell per step:

```python
from google.colab import drive; drive.mount("/content/drive")
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v1.py" setup
```

| Step | What it does | Time |
|---|---|---|
| `setup` | self-test, Phase 0 (freeze + integrity checks), pre-registration | 5-10 min |
| `phase1` | all 25 splits; `phase1 0-4` runs splits 0-4 only | several hours, see below |
| `report` | tables, pre-registered endpoint, figure | ~1 min |

The launcher unpacks the zip to `/content/cpaudit_v1` (local), re-unpacks it when you upload a
newer zip, and installs missing packages once per session. Phase 1 resumes where it stopped:
every finished (split, target) unit is saved atomically and HPO trials are saved after each
trial, so after a disconnect you simply run the same cell again.

## Nothing in your Drive folder is modified

* The code is unzipped to `/content/cpaudit_v1` (local), never into Drive.
* Inputs are only **read**: `soap_features_corrected_v2.npz`, `data_prep_v5.npz`
  (alignment check), `duplicate_entries_v7.csv` (cross-check), and `figstyle.py` (figure style).
* Every output goes to a **new** folder `NOMAD_Uncer_Quan/phase01_v1/`:

```
phase01_v1/
  inputs/train.csv                    downloaded once from the Sutton et al. repository
  phase0/freeze_manifest_v7.json      SHA-256 of every existing file (Task 0.2)
  phase0/integrity_report.json        split, alignment, duplicate and group checks
  phase0/environment.json             package versions, CPUs, RAM, GPU
  phase0/PREREGISTRATION.md / .json   hypotheses, endpoint, decision rules, hashes (Task 0.3)
  splits/split_XX.npz, _sizes.csv     NOMAD ids of every partition of every split
  hpo/hp_<target>_splitXX.json        tuned hyperparameters (+ per-trial log)
  runs/rXX_<target>_<tag>.csv/.npz    per-split metrics and all intervals / predictions
  report/PHASE1_SUMMARY.md            tables, endpoint decision, HPO sensitivity
  report/phase1_summary.csv           aggregated over splits
  report/fig_phase1_lattice_coverage  .png (600 dpi) + .pdf, in your figstyle
```

## What Phase 1 does

**Data (Tasks 1.1-1.2).** The 31 duplicate entries (15 groups) are reduced to 15, keeping
the lowest NOMAD id: n = 2,384. A pair counts as a duplicate only if the tabular record and the
SOAP vector are both identical. Phase 0 checks the result against `duplicate_entries_v7.csv`.
The unit of exchangeability is the **composition group** (lattice, Al and Ga fraction rounded
to 0.1), the same labels as your v7 leakage table (231 groups in the v5 test set, reproduced).
There are 366 groups in total, so each lattice has about 52-64 groups.

**Splits.** 25 independent outer splits from `StratifiedGroupKFold` (fold 0 = test, 20 %),
grouped by composition and stratified by lattice. Inside each training pool, a 25 %
calibration set and the CV+ folds are also composition-grouped. No composition group is ever
on both sides of any boundary; this is asserted in code for every split.

**Model.** LightGBM, L2 objective, `max_bin = 63`, the same search space as manuscript
Table 1. There are two corrections:
* **Bagging now works.** v5 set `subsample` without `subsample_freq`, and LightGBM bags only
  when `bagging_freq > 0`, so the reported subsample ratios had no effect. Here
  `bagging_freq = 1`.
* **The refit uses the tuned number of trees.** v5 refit with a fixed 700 trees after tuning
  with early stopping. Here the refit uses the mean best iteration of the selected trial.

Hyperparameters are tuned with Optuna TPE (25 trials, MAE, 5-fold grouped CV) on the
**proper-training set of split 0 only**, so they never see calibration or test data of that
split. Splits 1 and 2 are re-tuned on their own data (`tag = hpo_own`) to measure how much
coverage depends on the tuning.

**Interval methods**, all from the same fits so differences reflect the calibration only:

| Method | Calibration | Reference |
|---|---|---|
| `cvplus_random` | CV+, 5 composition-grouped random folds (Regime A) | Barber et al., Ann. Stat. 2021 |
| `cvplus_lattice` | CV+, one fold per lattice (your submitted design), as a reference | same |
| `split_global` | split conformal, one quantile | Lei et al., JASA 2018 |
| `split_mondrian` | split conformal, one quantile per lattice | Vovk, Mach. Learn. 2013; Bostrom, PMLR 2022 |
| `split_normalized_knn` | split, scores scaled by a kNN distance in PCA-50 SOAP space | Papadopoulos et al., JAIR 2011 |
| `cqr` | conformalized quantile regression, alpha = 0.10 | Romano et al., NeurIPS 2019 |

Metrics per split, method, alpha (0.05/0.10/0.20) and lattice, plus the two Fd-3m cell
settings (gamma ~30 deg: 247 entries; gamma ~120 deg: 146 entries):
* coverage with a Clopper-Pearson CI;
* mean and median width, and the fraction of infinite intervals;
* interval score (Gneiting & Raftery, JASA 2007);
* MAE, and Spearman correlation between width and |error|;
* the number of test and calibration composition groups;
* for split methods, the probability that a *valid* procedure would fall below nominal on
  that split (Beta-binomial; Vovk 2013).

**Pre-registration.** `prereg` writes the hypotheses, the primary endpoint (mean Fd-3m
coverage under `cvplus_random` at alpha = 0.10) and the decision thresholds before any result
exists. `phase1` refuses to run without it, and refuses a changed config unless you pass
`--allow-deviation`, which logs the change.

## Verification already done (in this bundle)
* `selftest` (no data needed):
  * split, Mondrian, normalized and CQR coverage equal the finite-sample value k/(n+1) in
    simulation;
  * CV+ matches the Barber et al. definition by brute force;
  * split, normalized and Mondrian intervals are identical to **crepes 0.9.1**;
  * CV+ intervals are identical to **MAPIE 0.8.6 `plus`** (max difference 0.0).
* A full dry run (phase0 -> prereg -> phase1 -> report -> resume) was executed on a mock
  copy of your folder built from the real `train.csv`. That run confirmed:
  * the v5 split is reproduced exactly (row order verified against your
    `fig_parity_v6_bandgap.csv`);
  * the duplicate groups match the v7 table exactly.

## Compute and the A100
LightGBM runs on CPU. The A100 runtime helps because it comes with more CPU cores and RAM;
the GPU itself is idle. LightGBM's GPU mode (OpenCL) needs a custom build, uses single-precision
histograms (results differ from CPU), and is rarely faster at 1.5k rows. Phase 0 tests whether
your LightGBM build can use the GPU and reports it in `environment.json`. Leave
`lgbm.device: cpu` unless you have a reason to change it.

A 700-tree fit on 1,500 x 23,976 features took about 2-3 min on 2 cores here. With 3
parallel 4-thread fits on a 12-core host, expect roughly:
* HPO: about 30-45 min per target for split 0, plus the same for each check split;
* Phase 1: about 5-8 min per (split, target), so **about 6-10 h for all 25 splits**.

Run it in chunks (`--repeats 0-4`, `5-9`, ...) if your sessions are short. Keep the browser
tab open; Colab Pro stops idle sessions.

## Task 0.1 (no code)
Before Phase 1 results are interpreted, give an outside reader the central question, the
hypotheses table and `phase0/PREREGISTRATION.md`, and ask: is this worth answering, and which
result would surprise you?

## Interpreting the numbers honestly
* **Correlated test points.** Points in one composition group are correlated. The effective
  test size per lattice is about 10-13 groups per split, not about 75 points. Read coverage
  together with `n_test_comp_groups`; the Beta-binomial reference understates the spread.
* **The guarantee is marginal over random group assignments.** With about 12 calibration
  groups per lattice, a single split can be far from nominal; the 25-split distribution is
  what to report.
* **`cvplus_lattice` is not a valid in-distribution design.** It is kept only to show, on
  identical data, how much of the manuscript's over-coverage comes from lattice-held-out
  calibration.
