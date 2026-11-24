# Data

## Inputs not included in this repository

| File (Drive layout) | Content | Size | SHA-256 |
|:--|:--|--:|:--|
| `phase01_v1/inputs/train.csv` | NOMAD 2018 train.csv (Sutton et al.) | 0.2 MB | `2be6168b578c4be91b7380e7404a2b58912790d8e03c51b28f4f279f2d8bb872` |
| `soap_features_corrected_v2.npz` | NOMAD SOAP features, 2,400 x 23,976 (release asset) | 354.5 MB | `73121e8adab96213f4494ab8e83a370bff44f4982d241ddf7e29ea77cc226d35` |
| `phase06_v1/features/features_v1.npz` | JARVIS-DFT features (release asset) | 56.2 MB | `70e2adad586edec749e6aa9ce3c14fdbb0b02cb691aea7a16c30c901131034d9` |
| `phase06_v1/inputs/jdft_3d-9-24-2025.json.zip` | JARVIS-DFT release used (figshare 64391379) | 48.4 MB | `f9e0a3309f0000d5de1ec9e49c93963109ea45e63f451103f8f6595d2eabf7f5` |
| `phase06_v1/inputs/jdft_3d-12-12-2022.json.zip` | JARVIS-DFT 2022 release, counted only (figshare 38521619) | 40.8 MB | `d4c64660e9e1fa45c82bd8868a96ec10162195eed69636445972c05550d8d0d6` |

* NOMAD 2018 `train.csv`: https://raw.githubusercontent.com/csutton7/nomad_2018_kaggle_dataset/master/train.csv (copies were stored in `phase01_v1`, `phase02_v1`, `phase04_v1` and `phase05_v1/inputs`). Cite Sutton et al., npj Comput. Mater. 5, 111 (2019).
* JARVIS-DFT: https://ndownloader.figshare.com/files/64391379 (used) and https://ndownloader.figshare.com/files/38521619 (2022 release, counted only). Cite Choudhary et al., npj Comput. Mater. 6, 173 (2020).
* `soap_features_corrected_v2.npz` and `features_v1.npz`: GitHub release assets of this repository. Put them at `results/soap_features_corrected_v2.npz` and `results/phase06_v1/features/features_v1.npz`. The SOAP recipe is in `code/nomad_soap/make_soap_features.py`; a recomputation may differ in the last digits, so compare with the SHA-256 above before reusing archived models.

Derived per-structure tables of JARVIS-DFT (identifiers, space groups, targets, predictions) are included under `results/phase06_v1` to `phase08_v1`.

## Files left out of `results/`

| Rule | Reason |
|:--|:--|
| `*/inputs/*` | third-party input data (NOMAD train.csv, JARVIS-DFT release zips); download with tools/fetch_inputs.py, SHA-256 in docs/DATA.md |
| `*PROJECT_LOG*` | internal working notes |
| `phase06_v1/features/features_v1.npz` | 56 MB JARVIS-DFT feature matrix: release asset (recomputation is not bit-identical across platforms) |
| `phase06_v1/features/composition_block.npz` | intermediate checkpoint of features6 |
| `phase06_v1/features/structure_chunks/*` | intermediate checkpoints of features6 |
| `phase06_v1/splits_matfold/*/*.csv` | 68 MB of MatFold split tables; the MatFold JSON split definitions are kept, and matfold6 recreates the tables |
| `phase08_v1/hjplus/parts/*` | intermediate checkpoints of hjplus8; merged results are in phase08_v1/hjplus/ |
| `*__pycache__*` | Python cache |
| `*.pyc` | Python cache |
| `*.ipynb_checkpoints*` | notebook cache |
| `*desktop.ini` | system file |
| `*.DS_Store` | system file |
| `*Thumbs.db` | system file |
