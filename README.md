# Conformal prediction intervals for materials in unseen crystal families

Code, decision-rule files, result tables and figures for

> Johaimen M. Omar and Muhammed Tan, "Conformal prediction intervals for materials in unseen crystal families: failure, detection and few-shot recalibration", manuscript submitted for publication.

The study asks whether conformal prediction intervals still cover materials from crystal families
that are absent from training. It uses NOMAD 2018 (Al-Ga-In sesquioxides, 6 lattices) and
JARVIS-DFT (59 space groups) with LightGBM models of the formation energy and the band gap.

## Contents

| Path | Content |
|:--|:--|
| `results/phase01_v1` ... `results/phase08_v1` | every result file behind the paper, in the folder layout the code reads |
| `figures/` | Figures 1-8 and S1-S7, drawn from `results/` by `code/cpaudit_v9` (map: `docs/FIGURES.md`) |
| `code/cpaudit_v1` ... `code/cpaudit_v9` | the analysis package, one folder per bundle, exactly as run |
| `code/bundles/` | the bundle zips and launchers, byte-identical to the files that were run |
| `code/nomad_soap/` | the SOAP descriptor recipe for NOMAD |
| `tools/` | download of the third-party inputs; figure rendering from `results/` |
| `docs/PROVENANCE.md` | which bundle produced which folder; decision-rule files and their hashes |
| `docs/DATA.md` | data sources, licences, SHA-256 of all inputs, excluded files |
| `MANIFEST.sha256` | SHA-256 and size of every file in this repository |

## Data

The third-party data are not redistributed here. `tools/fetch_inputs.py` downloads them and
checks their SHA-256.

* NOMAD 2018 (Sutton et al., npj Comput. Mater. 5, 111, 2019): `train.csv` from
  https://github.com/csutton7/nomad_2018_kaggle_dataset.
* JARVIS-DFT (Choudhary et al., npj Comput. Mater. 6, 173, 2020): release
  `jdft_3d-9-24-2025.json` from figshare (file 64391379), the release that jarvis-tools loads as `dft_3d`.

Two large derived inputs are attached to the GitHub release (see `docs/DATA.md`):
`soap_features_corrected_v2.npz` (NOMAD SOAP features, 2,400 x 23,976) and
`features_v1.npz` (JARVIS-DFT features, 81,884 x 270). The JARVIS-DFT features are needed to rerun
the JARVIS-DFT models exactly, because recomputing them gives small floating-point differences
across platforms, and the code refuses inputs whose hash changed.

## How the results were produced

Every step ran in Google Colab (Python 3.13, 8 CPUs, CPU-only LightGBM) through a one-line
launcher, `python cpaudit_run_vN.py <step>`. Each step wrote to a new folder and refused to
overwrite earlier output.

| Bundle | Steps | Output | Paper |
|:--|:--|:--|:--|
| `cpaudit_v1` | setup, phase1, report | `phase01_v1` | Regime A on NOMAD: Fig. 2a-d, Table S2 |
| `cpaudit_v2` | setup, post1, overlap, phase2 (band gap, formation), report2 | `phase02_v1` | de-duplication, Regime B on NOMAD, covariate-shift checks: Figs. 2e-g, 3, 4, S1, S2; Tables S3-S5 |
| `cpaudit_v3` | setup, phase4, report4 | `phase04_v1` | failure detection on NOMAD: Fig. 6a-c, Figs. S3, S4, Table S6 |
| `cpaudit_v4` | setup, phase5, report5, train5 | `phase05_v1` | few-shot recalibration on NOMAD: Fig. 7, Fig. S5, Table S7 |
| `cpaudit_v5` | posthoc5, audit6 | `phase05_v1/posthoc, phase06_v1/audit` | post hoc NOMAD tables (Table S8); JARVIS-DFT data audit (Table S9a, b) |
| `cpaudit_v6` | freeze6, features6 (first splits6, report6) | `phase06_v1/data, features, splits (v1 files), report` | JARVIS-DFT data set and features |
| `cpaudit_v61` | splits6, matfold6, report6 | `phase06_v1/splits (v2 files), splits_matfold, report (v2 files)` | JARVIS-DFT splits used by all models (Table S9c) |
| `cpaudit_v7` | setup7, tune7, fit7, report7 | `phase07_v1` | Regime A and B on JARVIS-DFT, group-level guarantees: Fig. 5, Tables S10, S11 |
| `cpaudit_v8` | setup8, scores8, fewshot8, hjplus8, report8 | `phase08_v1` | detection, few-shot and HJ+ on JARVIS-DFT: Figs. 6d-f, 8, S6, S7; Tables S10, S12 |
| `cpaudit_v81` | posthoc8 | `phase08_v1/posthoc` | post hoc JARVIS-DFT analyses: Tables S12, S13 |
| `cpaudit_v9` | figures9 | `figures/` | all figures (drawn from results/ only) |

## Decision rules

The files `results/*/phase0/PREREGISTRATION*.md` and `.json` are the decision-rule files described
in the Methods and listed in Supplementary Table S14. The code wrote each file once, before the
results it judges, and stored the SHA-256 of the configuration and of the code in it. Later steps
refused to run if the configuration, the code or an input had changed. The file names come from
the code; the rules were not deposited in an external registry. Analyses added after results were
seen are in folders named `posthoc` (and `post1`) and are labelled post hoc or exploratory in the
paper. `docs/PROVENANCE.md` matches every recorded code hash to an archived bundle.

Comments in the code and in the first rule files call the lattice-fold CV+ reference "the
submitted design". This is the reference that, as stated in the Methods, had been computed on
NOMAD before any decision rule was written.

## Reproducing

```bash
pip install -r requirements.txt
pip install --no-deps matminer==0.10.1 MatFold==1.3.3   # their pins conflict with the stack above
python tools/fetch_inputs.py                             # NOMAD train.csv into results/*/inputs
python tools/render_figures.py --out /tmp/figures_check  # all figures from results/ (minutes)
```

A full rerun follows the order of the table above, starting from a folder that holds only the
inputs, with `paths.drive_root` in each bundle's configuration pointing to that folder.
`code/cpaudit_vN/README.md` describes the steps of each bundle. Phase 1 alone took about 19 h on
8 CPUs. A rerun writes new decision-rule files; compare their configuration and code hashes with
the archived ones in `docs/PROVENANCE.md`.

## Software

`requirements.txt` pins the versions recorded by the code at run time (`results/*/phase0/environment.json`).
Figures were checked with matplotlib 3.10 and 3.11. The NOMAD SOAP features were computed with
DScribe 2.1 (`code/nomad_soap/`).

## Licence

Code: MIT (`LICENSE`), except `code/*/cpaudit/pubstyle.py`, which adapts the figure style of
figures4papers by C. Liu and is shared under CC BY-NC 4.0. Results and figures: CC BY 4.0.
Third-party data keep their own terms (see `docs/DATA.md`).

## Citation

See `CITATION.cff`. Please also cite NOMAD 2018 and JARVIS-DFT.
