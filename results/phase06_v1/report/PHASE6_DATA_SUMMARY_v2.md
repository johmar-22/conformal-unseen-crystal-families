# Phase 6 data report (Tasks 6.2-6.3; no target used)

## Freeze

* 81884 entries of scope `all_chemistries` (75205 in 59 eligible space groups); TBmBJ gap present for 18391. Release jdft_3d-9-24-2025.json (sha256 f9e0a3309f00). Frozen table sha256 c8f9a4777277.
* Check: row count, eligible space groups and their sizes equal the Task 6.1 audit.

## Features

* Matrix 81884 x 270: composition 145 (Magpie, Ward 2016), density 3, CrystalNN 122 (matminer 0.10.1, pymatgen 2026.10.24); 43 min.
* Missing values: 137 columns contain NaN (max share 0.0047); 0 all-NaN columns; rows with any NaN by block {'composition': 270, 'crystalnn': 27, 'density': 386}.
* Constant columns: 0 [].
* Structure warnings: 413 (CrystalNN failed for 27).
* Feature-identical entries: 0 rows in 0 sets (0 sets span more than one space group).

## Splits (MatFold protocol, Witman & Schindler 2025)

* Labels: SpacegroupAnalyzer(symprec = 0.1), as MatFold computes them; primary = MatFold labels (`losgo_mf`, `locso_mf`), JARVIS database labels as sensitivity (`losgo_db`, `locso_db`).
* LOSGO over 59 space groups with at least 200 entries (MatFold fraction_lower_limit 0.002436; rarer groups always in training); same groups under both label systems and the audit: True. LOCSO over 7 crystal systems.
* MatFold itself, run on all 81884 structures: label mismatches {'n': 81884, 'sgnum_mismatch': 0, 'crystalsys_mismatch': 0, 'composition_mismatch': 0}; its leave-one-out splits equal ours row by row: space groups True (59 folds), crystal systems True (7 folds).
* MatFold label = JARVIS label for 0.9783 of entries (crystal system 0.9812); per-group test-set Jaccard, MatFold vs JARVIS labels: median 0.969, minimum 0.776 (space group 8).
* Near-duplicate twins (same formula and MatFold space group, different JARVIS space group): 920 rows; with JARVIS labels 840 held-out rows have their twin in training; with MatFold labels each such pair stays on one side of every split.
* Polymorph exclusion (our extension, not part of MatFold) drops a median of 677 training entries per held-out group (range 52-3116).

|   label |   n_test |   n_train |   n_dropped_polymorph |   n_train_polymorph_excluded |   share_label_differs_other_system |   n_moves_in_vs_other_system | spg_symbol   | crys         |   jaccard_test_mf_vs_db |
|--------:|---------:|----------:|----------------------:|-----------------------------:|-----------------------------------:|-----------------------------:|:-------------|:-------------|------------------------:|
|       1 |     1129 |     80755 |                  1012 |                        79743 |                              0.017 |                           19 | P1           | triclinic    |                   0.801 |
|       2 |     1623 |     80261 |                  1124 |                        79137 |                              0.012 |                           19 | P-1          | triclinic    |                   0.9   |
|       4 |      337 |     81547 |                   402 |                        81145 |                              0.015 |                            5 | P2_1         | monoclinic   |                   0.912 |
|       5 |      543 |     81341 |                   675 |                        80666 |                              0.077 |                           42 | C2           | monoclinic   |                   0.861 |
|       6 |      790 |     81094 |                   970 |                        80124 |                              0.032 |                           25 | Pm           | monoclinic   |                   0.792 |
|       8 |     1052 |     80832 |                  1000 |                        79832 |                              0.096 |                          101 | Cm           | monoclinic   |                   0.776 |
|       9 |      273 |     81611 |                   384 |                        81227 |                              0.026 |                            7 | Cc           | monoclinic   |                   0.917 |
|      10 |      316 |     81568 |                   556 |                        81012 |                              0.044 |                           14 | P2/m         | monoclinic   |                   0.801 |
|      11 |      756 |     81128 |                   451 |                        80677 |                              0.019 |                           14 | P2_1/m       | monoclinic   |                   0.949 |
|      12 |     2663 |     79221 |                  2140 |                        77081 |                              0.039 |                          104 | C2/m         | monoclinic   |                   0.888 |
|      14 |     1543 |     80341 |                  1316 |                        79025 |                              0.007 |                           11 | P2_1/c       | monoclinic   |                   0.977 |
|      15 |     1328 |     80556 |                  1035 |                        79521 |                              0.019 |                           25 | C2/c         | monoclinic   |                   0.969 |
|      25 |      499 |     81385 |                   608 |                        80777 |                              0.156 |                           78 | Pmm2         | orthorhombic |                   0.803 |
|      31 |      249 |     81635 |                   348 |                        81287 |                              0.02  |                            5 | Pmn2_1       | orthorhombic |                   0.968 |
|      36 |      343 |     81541 |                   370 |                        81171 |                              0.026 |                            9 | Cmc2_1       | orthorhombic |                   0.952 |
|      38 |      985 |     80899 |                   887 |                        80012 |                              0.093 |                           92 | Amm2         | orthorhombic |                   0.883 |
|      44 |      336 |     81548 |                   339 |                        81209 |                              0.113 |                           38 | Imm2         | orthorhombic |                   0.849 |
|      47 |      527 |     81357 |                   863 |                        80494 |                              0.047 |                           25 | Pmmm         | orthorhombic |                   0.869 |
|      55 |      225 |     81659 |                   116 |                        81543 |                              0     |                            0 | Pbam         | orthorhombic |                   0.978 |
|      59 |      319 |     81565 |                   603 |                        80962 |                              0.031 |                           10 | Pmmn         | orthorhombic |                   0.942 |
|      62 |     2526 |     79358 |                  2066 |                        77292 |                              0.008 |                           19 | Pnma         | orthorhombic |                   0.987 |
|      63 |     1712 |     80172 |                  1701 |                        78471 |                              0.013 |                           23 | Cmcm         | orthorhombic |                   0.973 |
|      64 |      242 |     81642 |                   325 |                        81317 |                              0.033 |                            8 | Cmce         | orthorhombic |                   0.951 |
|      65 |      786 |     81098 |                   946 |                        80152 |                              0.084 |                           66 | Cmmm         | orthorhombic |                   0.901 |
|      71 |     1113 |     80771 |                   524 |                        80247 |                              0.063 |                           70 | Immm         | orthorhombic |                   0.903 |
|      74 |      582 |     81302 |                   792 |                        80510 |                              0.038 |                           22 | Imma         | orthorhombic |                   0.954 |
|      82 |      275 |     81609 |                   335 |                        81274 |                              0.007 |                            2 | I-4          | tetragonal   |                   0.989 |
|      87 |      212 |     81672 |                   199 |                        81473 |                              0.047 |                           10 | I4/m         | tetragonal   |                   0.927 |
|      99 |     3478 |     78406 |                  2302 |                        76104 |                              0.004 |                           13 | P4mm         | tetragonal   |                   0.989 |
|     107 |      580 |     81304 |                   266 |                        81038 |                              0.026 |                           15 | I4mm         | tetragonal   |                   0.969 |
|     115 |      318 |     81566 |                   204 |                        81362 |                              0.009 |                            3 | P-4m2        | tetragonal   |                   0.987 |
|     119 |      988 |     80896 |                   711 |                        80185 |                              0.017 |                           17 | I-4m2        | tetragonal   |                   0.972 |
|     121 |      306 |     81578 |                   177 |                        81401 |                              0.007 |                            2 | I-42m        | tetragonal   |                   0.993 |
|     122 |      244 |     81640 |                   334 |                        81306 |                              0.004 |                            1 | I-42d        | tetragonal   |                   0.992 |
|     123 |     4266 |     77618 |                  2418 |                        75200 |                              0.02  |                           87 | P4/mmm       | tetragonal   |                   0.968 |
|     127 |      404 |     81480 |                   119 |                        81361 |                              0.007 |                            3 | P4/mbm       | tetragonal   |                   0.985 |
|     129 |     1263 |     80621 |                   918 |                        79703 |                              0.01  |                           12 | P4/nmm       | tetragonal   |                   0.989 |
|     136 |      312 |     81572 |                   586 |                        80986 |                              0.029 |                            9 | P4_2/mnm     | tetragonal   |                   0.968 |
|     139 |     3985 |     77899 |                  1925 |                        75974 |                              0.019 |                           74 | I4/mmm       | tetragonal   |                   0.97  |
|     140 |      463 |     81421 |                   304 |                        81117 |                              0.011 |                            5 | I4/mcm       | tetragonal   |                   0.987 |
|     141 |      351 |     81533 |                   689 |                        80844 |                              0.014 |                            5 | I4_1/amd     | tetragonal   |                   0.969 |
|     148 |      719 |     81165 |                   555 |                        80610 |                              0.014 |                           10 | R-3          | trigonal     |                   0.977 |
|     156 |      824 |     81060 |                   504 |                        80556 |                              0.04  |                           33 | P3m1         | trigonal     |                   0.955 |
|     160 |      832 |     81052 |                   922 |                        80130 |                              0.072 |                           60 | R3m          | trigonal     |                   0.925 |
|     164 |     2039 |     79845 |                  1222 |                        78623 |                              0.026 |                           53 | P-3m1        | trigonal     |                   0.973 |
|     166 |     2267 |     79617 |                  1887 |                        77730 |                              0.03  |                           69 | R-3m         | trigonal     |                   0.961 |
|     167 |      315 |     81569 |                   408 |                        81161 |                              0.032 |                           10 | R-3c         | trigonal     |                   0.959 |
|     186 |      564 |     81320 |                   677 |                        80643 |                              0.023 |                           13 | P6_3mc       | hexagonal    |                   0.967 |
|     187 |     3487 |     78397 |                  3116 |                        75281 |                              0.007 |                           24 | P-6m2        | hexagonal    |                   0.993 |
|     189 |      822 |     81062 |                   485 |                        80577 |                              0.004 |                            3 | P-62m        | hexagonal    |                   0.996 |
|     191 |      951 |     80933 |                   910 |                        80023 |                              0.022 |                           21 | P6/mmm       | hexagonal    |                   0.978 |
|     193 |      359 |     81525 |                    99 |                        81426 |                              0.003 |                            1 | P6_3/mcm     | hexagonal    |                   0.997 |
|     194 |     2453 |     79431 |                  2865 |                        76566 |                              0.012 |                           30 | P6_3/mmc     | hexagonal    |                   0.987 |
|     215 |      229 |     81655 |                    52 |                        81603 |                              0.009 |                            2 | P-43m        | cubic        |                   0.991 |
|     216 |     4533 |     77351 |                  1259 |                        76092 |                              0.003 |                           15 | F-43m        | cubic        |                   0.997 |
|     221 |     4665 |     77219 |                  3034 |                        74185 |                              0.014 |                           65 | Pm-3m        | cubic        |                   0.986 |
|     225 |     8713 |     73171 |                  2978 |                        70193 |                              0.012 |                          108 | Fm-3m        | cubic        |                   0.987 |
|     227 |      934 |     80950 |                   983 |                        79967 |                              0.01  |                            9 | Fd-3m        | cubic        |                   0.99  |
|     229 |      249 |     81635 |                   319 |                        81316 |                              0.016 |                            4 | Im-3m        | cubic        |                   0.976 |

| label        |   n_test |   n_train |   n_dropped_polymorph |   n_train_polymorph_excluded |
|:-------------|---------:|----------:|----------------------:|-----------------------------:|
| cubic        |    20589 |     61295 |                  6281 |                        55014 |
| hexagonal    |     9317 |     72567 |                  6210 |                        66357 |
| monoclinic   |     9965 |     71919 |                  4918 |                        67001 |
| orthorhombic |    12876 |     69008 |                  6507 |                        62501 |
| tetragonal   |    18470 |     63414 |                  7411 |                        56003 |
| triclinic    |     2752 |     79132 |                  1698 |                        77434 |
| trigonal     |     7915 |     73969 |                  4232 |                        69737 |
