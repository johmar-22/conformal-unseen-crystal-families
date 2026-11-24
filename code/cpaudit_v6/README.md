# cpaudit_v6: Phase 6 Tasks 6.2-6.3 (JARVIS-DFT freeze, features, splits)

Upload `cpaudit_v6.zip` and `cpaudit_run_v6.py` to `MyDrive/NOMAD_Uncer_Quan/code_bundles/`
(next to v1-v5; nothing is replaced). One short Colab cell per step:

```python
from google.colab import drive; drive.mount("/content/drive", force_remount=True)
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v6.py" setup
```

| Step | What it does | Time (8 CPUs) |
|---|---|---|
| `setup` | self-test (all earlier tools plus the v6 data, feature and split tools) | 3-5 min |
| `freeze6` | analysis table rebuilt from the release on Drive; stops unless it equals the Task 6.1 audit | ~1 min |
| `features6` | Magpie 145 + density 3 + CrystalNN 122 features; parallel; checkpoint every 2,000 structures | 15-25 min |
| `splits6` | LOSGO (59 space groups) and LOCSO (7 crystal systems), polymorph-exclusion variant, MatFold check | 3-5 min |
| `report6` | one summary | seconds |
| `all` | the four steps above | |

No target is used and no model is fitted. Outputs: `phase06_v1/data/`, `features/`, `splits/`,
`report/`. Earlier outputs are only read.

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
* **Splits** (MatFold protocol; Witman & Schindler 2025): leave-one-space-group-out over the
  eligible groups, leave-one-crystal-system-out, and the polymorph-exclusion variant (training drops
  every entry whose reduced formula occurs in the held-out group). Group labels are the JARVIS labels
  the audit counted; MatFold re-detects space groups with SpacegroupAnalyzer(symprec = 0.1). We run
  MatFold itself on 2,000 random structures to verify that our re-detection equals its labels, then
  report per group how often MatFold's label differs from the database label.
