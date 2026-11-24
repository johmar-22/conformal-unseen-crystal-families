# Phase 6 Task 6.1: JARVIS-DFT data audit (dft_3d, 93902 entries)

Release: jdft_3d-9-24-2025.json (https://ndownloader.figshare.com/files/64391379), SHA-256 f9e0a3309f0000d5de1ec9e49c93963109ea45e63f451103f8f6595d2eabf7f5.

## Scopes after de-duplication by (reduced formula, space group)

| release     | scope                 |   n_raw |   n_dedup |   n_removed |   n_spacegroups |   K_eligible |   share_of_entries_in_eligible |   n_mbj_dedup_within_subset |   K_eligible_mbj |   n_mbj_in_eligible_groups |   median_largest_chemsys_share_eligible |
|:------------|:----------------------|--------:|----------:|------------:|----------------:|-------------:|-------------------------------:|----------------------------:|-----------------:|---------------------------:|----------------------------------------:|
| dft_3d      | all_chemistries       |   93902 |     81884 |       12018 |             215 |           59 |                          0.918 |                       19877 |               24 |                      16986 |                                   0.006 |
| dft_3d      | oxides_oxidation_rule |   17085 |     13885 |        3200 |             191 |           19 |                          0.604 |                        3218 |                1 |                       1637 |                                   0.01  |
| dft_3d      | oxides_element_rule   |   15341 |     12709 |        2632 |             187 |           22 |                          0.66  |                        2836 |                0 |                       1619 |                                   0.012 |
| dft_3d_2022 | all_chemistries       |   75993 |     66154 |        9839 |             213 |           55 |                          0.906 |                       18239 |               22 |                      15458 |                                   0.007 |
| dft_3d_2022 | oxides_oxidation_rule |   14268 |     11475 |        2793 |             187 |           17 |                          0.569 |                        3112 |                1 |                       1482 |                                   0.01  |
| dft_3d_2022 | oxides_element_rule   |   12579 |     10274 |        2305 |             182 |           15 |                          0.533 |                        2697 |                0 |                       1241 |                                   0.011 |

## Eligible space groups (n >= 200) by threshold

| release     |   min_group_size |   all_chemistries |   oxides_element_rule |   oxides_oxidation_rule |
|:------------|-----------------:|------------------:|----------------------:|------------------------:|
| dft_3d      |              100 |                82 |                    32 |                      33 |
| dft_3d      |              150 |                74 |                    25 |                      27 |
| dft_3d      |              200 |                59 |                    22 |                      19 |
| dft_3d      |              300 |                50 |                    14 |                      13 |
| dft_3d      |              500 |                36 |                     3 |                       7 |
| dft_3d_2022 |              100 |                77 |                    28 |                      31 |
| dft_3d_2022 |              150 |                61 |                    21 |                      21 |
| dft_3d_2022 |              200 |                55 |                    15 |                      17 |
| dft_3d_2022 |              300 |                41 |                    10 |                      10 |
| dft_3d_2022 |              500 |                32 |                     3 |                       4 |

## Oxide classification (O-containing compounds)

| oxide_status                  |   dft_3d |   dft_3d_2022 |
|:------------------------------|---------:|--------------:|
| no_charge_balanced_assignment |     3526 |          2599 |
| other_anion                   |     4565 |          3906 |
| oxide                         |    17085 |         14268 |

## Decision (rule fixed in the config before the run)

```
{
  "release": "dft_3d",
  "rule": "oxides primary only if K >= 20 under both oxide rules",
  "K": {
    "oxides_oxidation_rule": 19,
    "oxides_element_rule": 22,
    "all_chemistries": 59
  },
  "oxide_definition_flips_the_decision": true,
  "primary_scope": "all_chemistries",
  "secondary_scope": "oxides_oxidation_rule",
  "K_primary": 59,
  "K_secondary": 19,
  "h3_design_primary": "valid group-level split (>= 10 calibration groups + 49 training groups)",
  "h3_design_secondary": "group-level jackknife+ only (Lee et al. 2023; fewer than 20 eligible groups)",
  "tbmbj_primary": {
    "n": 19877,
    "K_eligible": 24
  }
}
```

Other releases: {"dft_3d_2022": {"K": {"oxides_oxidation_rule": 17, "oxides_element_rule": 15, "all_chemistries": 55}, "primary_scope": "all_chemistries"}}

## Eligible space groups of the primary scope (all_chemistries)

|   spg_number | spg_symbol   | crys         |    n |   n_mbj |   n_chemsys |   largest_chemsys_share |   share_formula_in_other_spg |
|-------------:|:-------------|:-------------|-----:|--------:|------------:|------------------------:|-----------------------------:|
|          225 | Fm-3m        | cubic        | 8607 |    2656 |        7642 |                   0     |                        0.18  |
|          221 | Pm-3m        | cubic        | 4600 |     939 |        3728 |                   0.001 |                        0.358 |
|          216 | F-43m        | cubic        | 4518 |    1425 |        4064 |                   0.001 |                        0.132 |
|          123 | P4/mmm       | tetragonal   | 4228 |     831 |        3625 |                   0.001 |                        0.345 |
|          139 | I4/mmm       | tetragonal   | 3958 |    1056 |        3546 |                   0.001 |                        0.257 |
|           99 | P4mm         | tetragonal   | 3490 |     719 |        2779 |                   0.001 |                        0.414 |
|          187 | P-6m2        | hexagonal    | 3464 |     499 |        3267 |                   0.002 |                        0.439 |
|           12 | C2/m         | monoclinic   | 2777 |     539 |        2226 |                   0.007 |                        0.342 |
|           62 | Pnma         | orthorhombic | 2520 |     388 |        2343 |                   0.001 |                        0.344 |
|          194 | P6_3/mmc     | hexagonal    | 2425 |     454 |        2100 |                   0.002 |                        0.525 |
|          166 | R-3m         | trigonal     | 2218 |     529 |        1939 |                   0.003 |                        0.35  |
|          164 | P-3m1        | trigonal     | 1988 |     363 |        1811 |                   0.004 |                        0.278 |
|            2 | P-1          | triclinic    | 1763 |     323 |        1381 |                   0.01  |                        0.281 |
|           63 | Cmcm         | orthorhombic | 1713 |     338 |        1554 |                   0.003 |                        0.308 |
|           14 | P2_1/c       | monoclinic   | 1557 |     313 |        1387 |                   0.007 |                        0.301 |
|            1 | P1           | triclinic    | 1367 |     245 |         812 |                   0.03  |                        0.361 |
|           15 | C2/c         | monoclinic   | 1320 |     258 |        1179 |                   0.005 |                        0.292 |
|          129 | P4/nmm       | tetragonal   | 1253 |     341 |        1238 |                   0.002 |                        0.344 |
|            8 | Cm           | monoclinic   | 1124 |     162 |         914 |                   0.008 |                        0.397 |
|           71 | Immm         | orthorhombic | 1085 |     184 |        1027 |                   0.003 |                        0.262 |
|          119 | I-4m2        | tetragonal   |  982 |     405 |         843 |                   0.003 |                        0.457 |
|            6 | Pm           | monoclinic   |  941 |     162 |         843 |                   0.007 |                        0.612 |
|          191 | P6/mmm       | hexagonal    |  930 |     165 |         816 |                   0.005 |                        0.343 |
|          227 | Fd-3m        | cubic        |  925 |     204 |         860 |                   0.003 |                        0.364 |
|           38 | Amm2         | orthorhombic |  919 |     158 |         854 |                   0.012 |                        0.348 |
|          189 | P-62m        | hexagonal    |  819 |      94 |         769 |                   0.005 |                        0.242 |
|          156 | P3m1         | trigonal     |  795 |     218 |         708 |                   0.015 |                        0.303 |
|          160 | R3m          | trigonal     |  775 |     195 |         705 |                   0.006 |                        0.403 |
|           11 | P2_1/m       | monoclinic   |  768 |     157 |         734 |                   0.004 |                        0.262 |
|           65 | Cmmm         | orthorhombic |  733 |     115 |         663 |                   0.008 |                        0.445 |
|          148 | R-3          | trigonal     |  716 |     213 |         687 |                   0.004 |                        0.265 |
|          107 | I4mm         | tetragonal   |  568 |     174 |         567 |                   0.004 |                        0.287 |
|           74 | Imma         | orthorhombic |  565 |      76 |         534 |                   0.005 |                        0.407 |
|          186 | P6_3mc       | hexagonal    |  557 |     143 |         538 |                   0.004 |                        0.469 |
|           47 | Pmmm         | orthorhombic |  553 |     163 |         457 |                   0.005 |                        0.658 |
|            5 | C2           | monoclinic   |  540 |     100 |         439 |                   0.022 |                        0.391 |
|          140 | I4/mcm       | tetragonal   |  459 |     104 |         448 |                   0.004 |                        0.255 |
|           25 | Pmm2         | orthorhombic |  446 |      97 |         412 |                   0.007 |                        0.543 |
|          127 | P4/mbm       | tetragonal   |  404 |      53 |         397 |                   0.005 |                        0.134 |
|           10 | P2/m         | monoclinic   |  363 |      91 |         325 |                   0.017 |                        0.614 |
|            4 | P2_1         | monoclinic   |  359 |      75 |         303 |                   0.042 |                        0.373 |
|          193 | P6_3/mcm     | hexagonal    |  358 |      10 |         351 |                   0.006 |                        0.148 |
|          141 | I4_1/amd     | tetragonal   |  352 |     112 |         348 |                   0.006 |                        0.636 |
|           36 | Cmc2_1       | orthorhombic |  342 |      91 |         328 |                   0.006 |                        0.31  |
|           59 | Pmmn         | orthorhombic |  318 |      72 |         304 |                   0.009 |                        0.572 |
|          115 | P-4m2        | tetragonal   |  316 |      82 |         301 |                   0.009 |                        0.193 |
|           44 | Imm2         | orthorhombic |  313 |      72 |         303 |                   0.013 |                        0.326 |
|          167 | R-3c         | trigonal     |  308 |      88 |         304 |                   0.006 |                        0.422 |
|          136 | P4_2/mnm     | tetragonal   |  304 |      80 |         290 |                   0.007 |                        0.428 |
|          121 | I-42m        | tetragonal   |  304 |      84 |         296 |                   0.01  |                        0.191 |
|            9 | Cc           | monoclinic   |  283 |      56 |         251 |                   0.025 |                        0.438 |
|           82 | I-4          | tetragonal   |  274 |     102 |         272 |                   0.007 |                        0.307 |
|           31 | Pmn2_1       | orthorhombic |  247 |      71 |         230 |                   0.02  |                        0.409 |
|          229 | Im-3m        | cubic        |  247 |      54 |         242 |                   0.008 |                        0.304 |
|          122 | I-42d        | tetragonal   |  244 |     135 |         241 |                   0.008 |                        0.434 |
|           64 | Cmce         | orthorhombic |  238 |      26 |         229 |                   0.008 |                        0.387 |
|           55 | Pbam         | orthorhombic |  230 |      36 |         225 |                   0.009 |                        0.174 |
|          215 | P-43m        | cubic        |  227 |      45 |         216 |                   0.009 |                        0.093 |
|           87 | I4/m         | tetragonal   |  208 |      46 |         207 |                   0.01  |                        0.236 |

## Crystal systems (secondary analysis; only 7 groups)

| crystal_system   |   all_chemistries |   oxides_element_rule |   oxides_oxidation_rule |
|:-----------------|------------------:|----------------------:|------------------------:|
| triclinic        |              3130 |                   811 |                    1293 |
| monoclinic       |             10399 |                  2934 |                    3943 |
| orthorhombic     |             12656 |                  2491 |                    2826 |
| tetragonal       |             18358 |                  2124 |                    1865 |
| trigonal         |              7722 |                  1671 |                    1858 |
| hexagonal        |              9232 |                   834 |                     592 |
| cubic            |             20387 |                  1844 |                    1508 |
