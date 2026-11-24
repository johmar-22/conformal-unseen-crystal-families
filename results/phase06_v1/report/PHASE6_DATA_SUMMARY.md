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

## Splits

* LOSGO over 59 eligible space groups (JARVIS labels); LOCSO over 7 crystal systems.
* MatFold check on 2000 random structures: our re-detection = MatFold's labels for 1.0000 (space group) and 1.0000 (crystal system).
* MatFold's label (symprec 0.1) equals the JARVIS label for 0.9783 of entries (crystal system 0.9812); eligible groups under MatFold labels: 59.
* Polymorph exclusion drops a median of 715 training entries per held-out group (range 54-3108).

|   label |   n_test |   n_train |   n_dropped_polymorph |   n_train_polymorph_excluded | spg_symbol   | crys         |   share_matfold_label_differs |   n_matfold_moves_in |
|--------:|---------:|----------:|----------------------:|-----------------------------:|:-------------|:-------------|------------------------------:|---------------------:|
|       1 |     1367 |     80517 |                  1354 |                        79163 | P1           | triclinic    |                         0.188 |                   19 |
|       2 |     1763 |     80121 |                  1453 |                        78668 | P-1          | triclinic    |                         0.09  |                   19 |
|       4 |      359 |     81525 |                   519 |                        81006 | P2_1         | monoclinic   |                         0.075 |                    5 |
|       5 |      540 |     81344 |                   776 |                        80568 | C2           | monoclinic   |                         0.072 |                   42 |
|       6 |      941 |     80943 |                  1279 |                        79664 | Pm           | monoclinic   |                         0.187 |                   25 |
|       8 |     1124 |     80760 |                  1219 |                        79541 | Cm           | monoclinic   |                         0.154 |                  101 |
|       9 |      283 |     81601 |                   411 |                        81190 | Cc           | monoclinic   |                         0.06  |                    7 |
|      10 |      363 |     81521 |                   715 |                        80806 | P2/m         | monoclinic   |                         0.168 |                   14 |
|      11 |      768 |     81116 |                   526 |                        80590 | P2_1/m       | monoclinic   |                         0.034 |                   14 |
|      12 |     2777 |     79107 |                  2432 |                        76675 | C2/m         | monoclinic   |                         0.079 |                  104 |
|      14 |     1557 |     80327 |                  1368 |                        78959 | P2_1/c       | monoclinic   |                         0.016 |                   11 |
|      15 |     1320 |     80564 |                  1053 |                        79511 | C2/c         | monoclinic   |                         0.013 |                   25 |
|      25 |      446 |     81438 |                   550 |                        80888 | Pmm2         | orthorhombic |                         0.056 |                   78 |
|      31 |      247 |     81637 |                   349 |                        81288 | Pmn2_1       | orthorhombic |                         0.012 |                    5 |
|      36 |      342 |     81542 |                   385 |                        81157 | Cmc2_1       | orthorhombic |                         0.023 |                    9 |
|      38 |      919 |     80965 |                   769 |                        80196 | Amm2         | orthorhombic |                         0.028 |                   92 |
|      44 |      313 |     81571 |                   295 |                        81276 | Imm2         | orthorhombic |                         0.048 |                   38 |
|      47 |      553 |     81331 |                   893 |                        80438 | Pmmm         | orthorhombic |                         0.092 |                   25 |
|      55 |      230 |     81654 |                   191 |                        81463 | Pbam         | orthorhombic |                         0.022 |                    0 |
|      59 |      318 |     81566 |                   600 |                        80966 | Pmmn         | orthorhombic |                         0.028 |                   10 |
|      62 |     2520 |     79364 |                  2106 |                        77258 | Pnma         | orthorhombic |                         0.005 |                   19 |
|      63 |     1713 |     80171 |                  1784 |                        78387 | Cmcm         | orthorhombic |                         0.014 |                   23 |
|      64 |      238 |     81646 |                   349 |                        81297 | Cmce         | orthorhombic |                         0.017 |                    8 |
|      65 |      733 |     81151 |                   885 |                        80266 | Cmmm         | orthorhombic |                         0.018 |                   66 |
|      71 |     1085 |     80799 |                   538 |                        80261 | Immm         | orthorhombic |                         0.039 |                   70 |
|      74 |      565 |     81319 |                   816 |                        80503 | Imma         | orthorhombic |                         0.009 |                   22 |
|      82 |      274 |     81610 |                   329 |                        81281 | I-4          | tetragonal   |                         0.004 |                    2 |
|      87 |      208 |     81676 |                   192 |                        81484 | I4/m         | tetragonal   |                         0.029 |                   10 |
|      99 |     3490 |     78394 |                  2317 |                        76077 | P4mm         | tetragonal   |                         0.007 |                   13 |
|     107 |      568 |     81316 |                   249 |                        81067 | I4mm         | tetragonal   |                         0.005 |                   15 |
|     115 |      316 |     81568 |                   200 |                        81368 | P-4m2        | tetragonal   |                         0.003 |                    3 |
|     119 |      982 |     80902 |                   745 |                        80157 | I-4m2        | tetragonal   |                         0.011 |                   17 |
|     121 |      304 |     81580 |                   177 |                        81403 | I-42m        | tetragonal   |                         0     |                    2 |
|     122 |      244 |     81640 |                   340 |                        81300 | I-42d        | tetragonal   |                         0.004 |                    1 |
|     123 |     4228 |     77656 |                  2509 |                        75147 | P4/mmm       | tetragonal   |                         0.012 |                   87 |
|     127 |      404 |     81480 |                   120 |                        81360 | P4/mbm       | tetragonal   |                         0.007 |                    3 |
|     129 |     1253 |     80631 |                   921 |                        79710 | P4/nmm       | tetragonal   |                         0.002 |                   12 |
|     136 |      304 |     81580 |                   575 |                        81005 | P4_2/mnm     | tetragonal   |                         0.003 |                    9 |
|     139 |     3958 |     77926 |                  1999 |                        75927 | I4/mmm       | tetragonal   |                         0.012 |                   74 |
|     140 |      459 |     81425 |                   240 |                        81185 | I4/mcm       | tetragonal   |                         0.002 |                    5 |
|     141 |      352 |     81532 |                   691 |                        80841 | I4_1/amd     | tetragonal   |                         0.017 |                    5 |
|     148 |      716 |     81168 |                   556 |                        80612 | R-3          | trigonal     |                         0.01  |                   10 |
|     156 |      795 |     81089 |                   510 |                        80579 | P3m1         | trigonal     |                         0.005 |                   33 |
|     160 |      775 |     81109 |                   881 |                        80228 | R3m          | trigonal     |                         0.004 |                   60 |
|     164 |     1988 |     79896 |                  1207 |                        78689 | P-3m1        | trigonal     |                         0.001 |                   53 |
|     166 |     2218 |     79666 |                  1941 |                        77725 | R-3m         | trigonal     |                         0.009 |                   69 |
|     167 |      308 |     81576 |                   406 |                        81170 | R-3c         | trigonal     |                         0.01  |                   10 |
|     186 |      557 |     81327 |                   672 |                        80655 | P6_3mc       | hexagonal    |                         0.011 |                   13 |
|     187 |     3464 |     78420 |                  3108 |                        75312 | P-6m2        | hexagonal    |                         0     |                   24 |
|     189 |      819 |     81065 |                   487 |                        80578 | P-62m        | hexagonal    |                         0     |                    3 |
|     191 |      930 |     80954 |                   879 |                        80075 | P6/mmm       | hexagonal    |                         0     |                   21 |
|     193 |      358 |     81526 |                   100 |                        81426 | P6_3/mcm     | hexagonal    |                         0     |                    1 |
|     194 |     2425 |     79459 |                  2848 |                        76611 | P6_3/mmc     | hexagonal    |                         0.001 |                   30 |
|     215 |      227 |     81657 |                    54 |                        81603 | P-43m        | cubic        |                         0     |                    2 |
|     216 |     4518 |     77366 |                  1238 |                        76128 | F-43m        | cubic        |                         0     |                   15 |
|     221 |     4600 |     77284 |                  3065 |                        74219 | Pm-3m        | cubic        |                         0     |                   65 |
|     225 |     8607 |     73277 |                  3011 |                        70266 | Fm-3m        | cubic        |                         0     |                  108 |
|     227 |      925 |     80959 |                   984 |                        79975 | Fd-3m        | cubic        |                         0     |                    9 |
|     229 |      247 |     81637 |                   322 |                        81315 | Im-3m        | cubic        |                         0.008 |                    4 |

| label        |   n_test |   n_train |   n_dropped_polymorph |   n_train_polymorph_excluded |
|:-------------|---------:|----------:|----------------------:|-----------------------------:|
| cubic        |    20387 |     61497 |                  6373 |                        55124 |
| hexagonal    |     9232 |     72652 |                  6197 |                        66455 |
| monoclinic   |    10399 |     71485 |                  5643 |                        65842 |
| orthorhombic |    12656 |     69228 |                  6459 |                        62769 |
| tetragonal   |    18358 |     63526 |                  7494 |                        56032 |
| triclinic    |     3130 |     78754 |                  2281 |                        76473 |
| trigonal     |     7722 |     74162 |                  4291 |                        69871 |
