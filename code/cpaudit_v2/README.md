# cpaudit_v2: Phase 1 follow-up, Phase 2 (LOSO mechanism), Phase 3 (overlap)

Upload `cpaudit_v2.zip` and `cpaudit_run_v2.py` to `MyDrive/NOMAD_Uncer_Quan/code_bundles/`
(next to the v1 files; nothing is replaced). One short Colab cell per step:

```python
from google.colab import drive; drive.mount("/content/drive")
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v2.py" setup
```

| Step | What it does | Time (8 CPUs) |
|---|---|---|
| `setup` | self-test, v2 de-duplication + integrity report, Phase 2/3 pre-registration | 5-10 min |
| `post1` | exploratory checks on saved Phase 1 outputs (no refitting) | 5-15 min |
| `overlap` | Phase 3 overlap diagnostics (must run before `phase2`) | ~5 min |
| `phase2 bandgap` | nested LOSO (6 holdouts) + leave-two-out (5), band gap | ~7-9 h |
| `phase2 formation` | same, formation energy | ~9-11 h |
| `report2` | tables, pre-registered decisions, figure | ~1 min |

`phase2` alone runs both targets. Every unit and every HPO trial is saved as it finishes, so
after a disconnect run the same cell again. Inputs are only read (SOAP npz, `phase01_v1/`,
`figstyle.py`); all outputs go to the new folder `NOMAD_Uncer_Quan/phase02_v1/`:

```
phase02_v1/
  phase0/   environment.json, dedup_v2.json, dedup_table_v2.csv, integrity_report_v2.json,
            PREREGISTRATION_phase2.md/.json, DEVIATIONS.json (only if a change is allowed)
  post1/    POST1_SUMMARY.md, dedup_sensitivity_*.csv, cluster_mondrian_*.csv, fd3m_subfamily_phase1.csv
  overlap/  overlap_diagnostics.csv, overlap_decision.json, overlap_weights.npz
  hpo/      hp_<target>_loso<sg>.json + per-trial logs (nested, lattice-out inner folds)
  loso/     <target>_loso<sg>.csv/.npz/_calcomp.csv/_grouplevel.csv, <target>_l2o<sg>.csv/.npz
  report/   PHASE2_SUMMARY.md, phase2_*.csv, phase2_decisions.json, fig_phase2_loso_coverage.png/.pdf
```

## Changes relative to v1

* **De-duplication.** v1 required equal tabular record and equal SOAP vector (rtol 1e-4) and
  removed 9 entries. v2 compares the input structures themselves (pymatgen StructureMatcher on
  the Vegard geometry files, setting- and order-independent) and treats same-input entries with
  both targets within 0.002 eV(/atom) as one record. Same-input entries with larger target
  differences are *conflicting records*: kept, flagged, and every coverage is also reported
  without them. Geometry files are read from `NOMAD_Uncer_Quan/train/` if present, otherwise
  downloaded once per session to local disk.
* loky's "A worker stopped" warning is silenced (it is a memory-leak guard; results unaffected).

## Methods (Phase 2)

Nested LOSO: hyperparameters tuned (Optuna TPE, 15 trials) with leave-one-lattice-out inner
folds over the five training lattices. CV+ with the lattice folds (Barber et al. 2021); nested
normalized and CQR scores inside CV+ (Gupta, Kuchibhotla & Ramdas 2022; Romano et al. 2019;
distance-based sigma as in Hu et al. 2022); in-distribution random-fold CV+ as contrast;
calibration-only ablation without the Fd-3m fold; leave-two-out refit; group-level threshold
(Dunn, Wasserman & Ramdas 2023). Phase 3: domain-classifier AUROC, Kish ESS and test-point
weight share; JAW-style weighted CV+ (Tibshirani et al. 2019; Prinster et al. 2022) only where
the pre-registered overlap gate passes.
