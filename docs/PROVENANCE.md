# Provenance

## Bundles

`code/bundles/` holds each bundle zip and launcher byte-for-byte. The code hash is the SHA-256 that the code itself records (file name and bytes of `cpaudit/*.py`, sorted by name).

| Bundle | Zip SHA-256 | Launcher SHA-256 | Code hash | Steps | Output |
|:--|:--|:--|:--|:--|:--|
| `cpaudit_v1` | `eead95f41820` | `dc343e1e4ef0` | `fa22ab7a0f41` | setup, phase1, report | `phase01_v1` |
| `cpaudit_v2` | `c5dd69d5124e` | `7e2e685b44a1` | `88eab1e03e62` | setup, post1, overlap, phase2 (band gap, formation), report2 | `phase02_v1` |
| `cpaudit_v3` | `345ba12ee79d` | `4b7ab394952e` | `8f8e4da22378` | setup, phase4, report4 | `phase04_v1` |
| `cpaudit_v4` | `39e6895a33f8` | `f69872434234` | `fbebadfd1099` | setup, phase5, report5, train5 | `phase05_v1` |
| `cpaudit_v5` | `f2a13a09a89c` | `75b1a4c122c0` | `e2b26efc5edb` | posthoc5, audit6 | `phase05_v1/posthoc, phase06_v1/audit` |
| `cpaudit_v6` | `b00e5cd46279` | `a7685ee22b2e` | `84c94726ebaf` | freeze6, features6 (first splits6, report6) | `phase06_v1/data, features, splits (v1 files), report` |
| `cpaudit_v61` | `9fa52da7f0ae` | `5c046721c7f6` | `768545a70e7e` | splits6, matfold6, report6 | `phase06_v1/splits (v2 files), splits_matfold, report (v2 files)` |
| `cpaudit_v7` | `dada94b09401` | `7592438363ca` | `43ee6e98c5ad` | setup7, tune7, fit7, report7 | `phase07_v1` |
| `cpaudit_v8` | `7a34ff48d3fa` | `a9da138d3a43` | `f260c5e34f8d` | setup8, scores8, fewshot8, hjplus8, report8 | `phase08_v1` |
| `cpaudit_v81` | `3ee8caa43310` | `a0fecee49ab4` | `c51dede456d0` | posthoc8 | `phase08_v1/posthoc` |
| `cpaudit_v9` | `87bdc334d896` | `e639cc476ce9` | `48bd48cdd2f1` | figures9 | `figures/` |

Full SHA-256 values are in `MANIFEST.sha256`.

## Decision-rule files

The configuration hash was recomputed from each archived bundle's own configuration file with that bundle's code; the code hash was recomputed from the bundle's files.

| File | Written (UTC) | Configuration hash | Recomputed from | Code hash | Bundle with this code hash |
|:--|:--|:--|:--|:--|:--|
| `results/phase01_v1/phase0/PREREGISTRATION.json` | 2026-11-16T08:15:53Z | `0e3daf052c14` | cpaudit_v1/config/phase01.yaml | `fa22ab7a0f41` | cpaudit_v1 |
| `results/phase02_v1/phase0/PREREGISTRATION_phase2.json` | 2026-11-17T07:14:00Z | `8b86e24ced22` | cpaudit_v2/config/phase02.yaml | `88eab1e03e62` | cpaudit_v2 |
| `results/phase04_v1/phase0/PREREGISTRATION_phase4.json` | 2026-11-17T19:13:28Z | `f654359b5ecc` | cpaudit_v3/config/phase04.yaml | `e411f88d5433` | none of the archived bundles |
| `results/phase05_v1/phase0/PREREGISTRATION_phase5.json` | 2026-11-17T20:12:01Z | `db2e1156ab5e` | cpaudit_v4/config/phase05.yaml | `fbebadfd1099` | cpaudit_v4 |
| `results/phase07_v1/phase0/PREREGISTRATION_phase6.json` | 2026-11-18T20:16:44Z | `baea14d0ddf1` | cpaudit_v7/config/phase06.yaml | `43ee6e98c5ad` | cpaudit_v7 |
| `results/phase08_v1/phase0/PREREGISTRATION_phase6_v8.json` | 2026-11-19T06:44:03Z | `452f4cb9e312` | not recomputed | `f260c5e34f8d` | cpaudit_v8 |

**Phase 4 note.** The bundle that wrote the Phase 4 rule file (code hash `e411f88d5433`) was replaced on Drive, after the Phase 4 run, by release 3.1.1 of `cpaudit_v3`, which changed a figure axis label and the wording of a confidence interval in the figure code only (run notes kept at the time). The archived `cpaudit_v3` is this later release; the earlier build was not kept. The Phase 4 configuration hash was recomputed from the archived bundle and equals the one in the rule file.

## Code hashes recorded in result files

| File | Key | Code hash | Bundle |
|:--|:--|:--|:--|
| `results/phase05_v1/posthoc/phase5_posthoc.json` | `code_hash` | `e2b26efc5edb` | cpaudit_v5 |
| `results/phase06_v1/audit/AUDIT_DECISION.json` | `code_hash` | `e2b26efc5edb` | cpaudit_v5 |
| `results/phase06_v1/data/freeze.json` | `code_hash` | `84c94726ebaf` | cpaudit_v6 |
| `results/phase06_v1/features/features_v1.json` | `code_hash` | `84c94726ebaf` | cpaudit_v6 |
| `results/phase06_v1/splits/splits_v1.json` | `code_hash` | `84c94726ebaf` | cpaudit_v6 |
| `results/phase06_v1/splits/splits_v2.json` | `code_hash` | `768545a70e7e` | cpaudit_v61 |
| `results/phase06_v1/splits_matfold/matfold6.json` | `code_hash` | `768545a70e7e` | cpaudit_v61 |
| `results/phase08_v1/posthoc/POSTHOC8.json` | `code_hash` | `c51dede456d0` | cpaudit_v81 |
