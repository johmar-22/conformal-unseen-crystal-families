# cpaudit_v5: Phase 5 post hoc tables, Phase 6 Task 6.1 (JARVIS-DFT audit), paper figures

Upload `cpaudit_v5.zip` and `cpaudit_run_v5.py` to `MyDrive/NOMAD_Uncer_Quan/code_bundles/`
(next to v1-v4; nothing is replaced). One short Colab cell per step:

```python
from google.colab import drive; drive.mount("/content/drive")
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v5.py" setup
```

| Step | What it does | Time |
|---|---|---|
| `setup` | self-test (conformal, few-shot, post hoc, audit rules, figure overlap check) | ~3-5 min |
| `posthoc5` | Phase 5 post hoc tables (NOT pre-registered): no-model baseline from the same labels, residual ICC | ~1 min |
| `audit6` | Phase 6 Task 6.1: JARVIS-DFT counts per space group and the scope decision | ~2-4 min |
| `figures` | Fig. 1-7, S1-S6 from the Drive files, each checked for overlaps | ~2 min |
| `all` | `posthoc5`, `audit6`, `figures` | |

Nothing refits a model. Existing Drive files are only read. New outputs:
`phase05_v1/posthoc/`, `phase06_v1/` (`inputs/` = the JARVIS zips, `audit/`), `paper_figures_v2/`.

## posthoc5 (not pre-registered)

* The Phase 5 label sets (one structure per composition) are re-created with the seeds frozen in
  `phase05_v1/phase0/PREREGISTRATION_phase5.json`; every re-created draw is compared with the
  saved Phase 5 draws (model widths and coverages must be identical) before anything is written.
* No-model interval: full conformal for y = b + e on the same labels (the offset-corrected
  interval with no model), the honest yardstick for "useful".
* ICC(1) of the residuals by composition group (ANOVA estimator, unequal group sizes) and its
  Spearman correlation (permutation p) with the pre-registered clustered dispersion ratio.

## audit6 (plan Task 6.1)

* Pinned release: `jdft_3d-9-24-2025.json` (the file jarvis-tools 2026.6.12 loads for "dft_3d",
  93,902 entries), plus the 2022 release (75,993) for comparison. Zips are stored on Drive and
  their SHA-256 recorded.
* Counts only; no target value is summarised. Duplicates by (reduced formula, space group) keep the
  lowest formation energy per atom.
* Oxides: primary rule = most probable charge-balanced oxidation states give a negative state to O
  only (pymatgen); sensitivity = element rule. Decision rule (fixed in `config/phase06.yaml`):
  oxides are primary only if both rules give >= 20 space groups with >= 200 structures; otherwise
  all chemistries are primary and oxides secondary.

## figures

Every figure is drawn from the Drive files only (the manifest lists each input with its SHA-256,
and the step stops if any input lies outside the Drive folder). After layout, `cpaudit/figcheck.py`
checks each figure for text on text, text on another panel, legends or labels covering data, and
labels covering data or spines; the result is printed and stored in `FIGURES_MANIFEST.json`.
Style: figures4papers (Chen Liu, CC BY-NC 4.0), adapted in `cpaudit/pubstyle.py`.
