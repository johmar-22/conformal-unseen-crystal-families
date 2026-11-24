# cpaudit_v4: Phase 5, few-shot family-specific recalibration

Upload `cpaudit_v4.zip` and `cpaudit_run_v4.py` to `MyDrive/NOMAD_Uncer_Quan/code_bundles/`
(next to v1-v3; nothing is replaced). One short Colab cell per step:

```python
from google.colab import drive; drive.mount("/content/drive")
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v4.py" setup
```

| Step | What it does | Time (8 CPUs) |
|---|---|---|
| `setup` | self-test, link the Phase 2 units (hashes recorded), Phase 5 pre-registration | ~5 min |
| `phase5` | few-shot draws: 7 families x 2 targets x 4 schemes x m grid x 500 draws; no refitting | ~3-6 min |
| `report5` | summaries, decisions D5.0, H5, H5c, H5b, practical m*, widths, sub-family table | ~1-3 min |
| `figures` | paper figures Fig. 1-7 and S1-S5 (read only) | ~1-2 min |
| `train5` | OPTIONAL: labels into training vs into calibration (LightGBM refits, resumable) | ~1-3 h |

Inputs are only read (`phase01_v1/`, `phase02_v1/`, `phase04_v1/`, the SOAP npz for `train5`);
outputs go to `NOMAD_Uncer_Quan/phase05_v1/` (`phase0/`, `draws/`, `train/`, `report/`) and
`paper_figures_v1/`.

## Design (pre-registered at `setup`; full text in phase05_v1/phase0/PREREGISTRATION_phase5.md)

* Families: the 6 LOSO lattices and the leave-two-out P63/mmc refit. Model: the fixed Phase 2
  predictor (mean of the five lattice-fold models). Residual z = y - mu.
* Methods: split conformal (|z|), and the offset-corrected interval = full conformal for
  y = mu + b + e. For this model full conformal equals jackknife+ with a leave-one-out offset
  (proof in cpaudit/fewshot.py, verified in the self-test), so it keeps the exact 1 - alpha
  guarantee.
* Schemes: random structures (exact finite-population reference law), clustered labels (whole
  composition groups), one structure per composition (Dunn et al. 2023), Fd-3m cell-setting
  transfer; permutation null for the group schemes.
* Decisions: D5.0 code check; H5 (spread of coverage for one-per-composition labels vs the
  permutation null); H5c (clustered labels harm coverage); H5b (offset correction narrows the
  intervals by >= 20 % in the cells where Phase 2 found an offset-dominated error).

Literature: Vovk (2012/2013) training-conditional validity; Angelopoulos & Bates (2023) Sec. 3.2;
Lei et al. (2018); Vovk, Gammerman & Shafer (2005) full conformal; Barber et al. (2021)
jackknife+; Dunn, Wasserman & Ramdas (2023) and Lee, Barber & Willett (2023) hierarchical data;
Fisch et al. (2021) few-shot conformal prediction.

## Paper figures

`cpaudit/paperfigs.py` in the figures4papers house style (Chen Liu, CC BY-NC 4.0), adapted in
`cpaudit/pubstyle.py`. Fig. 7 (few-shot recalibration), Fig. S4 (dispersion), Fig. S5 (cell-setting
transfer and, if run, the training arm) are new in v4.
