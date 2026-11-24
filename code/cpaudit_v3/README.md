# cpaudit_v3: Phase 4, detecting failure before it matters

Upload `cpaudit_v3.zip` and `cpaudit_run_v3.py` to `MyDrive/NOMAD_Uncer_Quan/code_bundles/`
(next to v1 and v2; nothing is replaced). One short Colab cell per step:

```python
from google.colab import drive; drive.mount("/content/drive")
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v3.py" setup
```

| Step | What it does | Time (8 CPUs) |
|---|---|---|
| `setup` | self-test, link the 22 Phase 2 units (hashes recorded), Phase 4 pre-registration | ~3-5 min |
| `phase4` | scores for every held-out structure in 6 LOSO + 5 leave-two-out configurations | ~20-40 min |
| `report4` | AUROC/AP with cluster-bootstrap CIs, abstention curves, lattice level, flags, decisions, figure | ~5-10 min |
| `figures` | paper figures (Fig. 1-6, S1-S3) from the saved Phase 1-4 outputs; read only | ~1-2 min |

No model is refitted. Inputs are only read (SOAP npz, `phase01_v1/`, `phase02_v1/`, `figstyle.py`); outputs go to
`NOMAD_Uncer_Quan/phase04_v1/` (`phase0/`, `scores/`, `report/`). Each configuration's feature
scores are saved as they finish, so after a disconnect run the same cell again.

## Design (pre-registered at `setup`)

* Outcome: miscoverage of the submitted design's interval (CV+ with lattice folds, alpha = 0.1)
  for the held-out lattice; secondary outcome |error| above the in-distribution 90th percentile.
* Scores: CV+ width; kNN distance (standardized SOAP, PCA-50, k = 25); KDE dissimilarity
  (Schultz et al., npj Comput. Mater. 2025: Epanechnikov kernel, nearest-neighbour bandwidth,
  all standardized features); lattice-fold committee SD; random-fold subsample-ensemble SD;
  domain-classifier probability.
* Analyses: pooled and within-lattice AUROC with a cluster bootstrap over composition groups;
  selective prediction (Geifman & El-Yaniv 2017); lattice-level Spearman with permutation p
  (descriptive); kNN-distance conformal p-value flag (Bates et al., Ann. Stat. 2023) with an
  in-distribution false-flag check.
* Decision rules H4 and G4, and the Phase 2 bias/one-sided-miss supplement, are in
  `phase04_v1/phase0/PREREGISTRATION_phase4.md` and `report/PHASE4_SUMMARY.md`.

## Paper figures (`figures` step, added in 3.1.0)

`cpaudit/paperfigs.py` draws the manuscript figures from files the pipeline has already written;
nothing is refitted and no Phase 4 code or setting is changed (the config hash is unchanged
because only `paths` keys were added). Style: the figures4papers house style by Chen Liu
(https://github.com/ChenLiu-1996/figures4papers, CC BY-NC 4.0), adapted in `cpaudit/pubstyle.py`
for a 180 mm journal page, Colab fonts (Liberation Sans) and no LaTeX.

| File | Content | Inputs |
|---|---|---|
| fig1_indistribution_coverage | per-split coverage by lattice (CV+ random, Mondrian, lattice folds); width and coverage, other lattices vs Fd-3m | phase01_v1 |
| fig2_group_level_mondrian | pooled vs one-per-group Mondrian; width cost | phase02_v1/post1 |
| fig3_unseen_lattice_coverage | LOSO coverage heatmap, 7 designs x 6 lattices, width in lattice SD | phase02_v1/report, loso |
| fig4_mechanism_offset_transfer | signed-error violins; fold-model transfer matrix; LOSO -> ablation -> refit | phase02_v1 |
| fig5_no_covariate_overlap | classifier log-odds ridgelines; test-point weight share and Kish ESS | phase02_v1/overlap |
| fig6_failure_detection | AUROC forest plots, abstention curves, lattice level, outlier flags | phase04_v1/report |
| figS1-S3 | calibration-tail composition; formation transfer; detection radar | phase02_v1, phase04_v1 |

Output: `NOMAD_Uncer_Quan/paper_figures_v1/` as PNG (600 dpi), PDF and SVG (editable text),
plus `FIGURES_MANIFEST.json` with the SHA-256 of every input. Figures whose inputs do not exist
yet are skipped, so the step can be run before Phase 4 (Fig. 1-5, S1, S2 only).
Requires matplotlib >= 3.7 (Colab has it).
