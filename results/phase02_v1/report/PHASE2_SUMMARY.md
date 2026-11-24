# Phase 2 / Phase 3 summary (alpha = 0.1)

## Phase 3 overlap diagnostics (computed before LOSO)

| heldout                     |   sg |   n_g |   n_rest |   auroc_logistic |   auroc_knn |   auroc_max |   ess_rest |   ess_frac |   max_weight_share |   median_test_weight_share |   frac_test_unbounded |   overlap_ok |
|:----------------------------|-----:|------:|---------:|-----------------:|------------:|------------:|-----------:|-----------:|-------------------:|---------------------------:|----------------------:|-------------:|
| Fd-3m                       |  227 |   393 |     1950 |            1     |       1     |       1     |   1034.09  |      0.53  |              0.005 |                      0.988 |                     1 |            0 |
| C2/m                        |   12 |   343 |     2000 |            1     |       0.999 |       1     |    462.637 |      0.231 |              0.018 |                      0.984 |                     1 |            0 |
| Pna21                       |   33 |   427 |     1916 |            1     |       1     |       1     |      2.845 |      0.001 |              0.525 |                      0.934 |                     1 |            0 |
| R-3c                        |  167 |   350 |     1993 |            1     |       1     |       1     |   1524.64  |      0.765 |              0.003 |                      0.986 |                     1 |            0 |
| P63/mmc                     |  194 |   351 |     1992 |            1     |       1     |       1     |   1034.08  |      0.519 |              0.006 |                      0.987 |                     1 |            0 |
| Ia-3                        |  206 |   479 |     1864 |            1     |       1     |       1     |    912.402 |      0.489 |              0.012 |                      0.984 |                     1 |            0 |
| 227a vs 227b (within Fd-3m) |  227 |   146 |      247 |            0.888 |       0.749 |       0.888 |    nan     |    nan     |            nan     |                    nan     |                   nan |          nan |

## bandgap: LOSO coverage of the held-out lattice

| heldout   |   cvplus_lattice |   cvplus_random_in |   ncvplus_knn |   cqr_cvplus |   cvplus_lattice_noFd3m_cal |   ncvplus_knn_noFd3m_cal |   cqr_cvplus_noFd3m_cal |
|:----------|-----------------:|-------------------:|--------------:|-------------:|----------------------------:|-------------------------:|------------------------:|
| C2/m      |            0.997 |              0.493 |         0.997 |        0.997 |                       0.878 |                    0.898 |                   0.942 |
| Pna21     |            1     |              0.988 |         1     |        1     |                       1     |                    1     |                   0.993 |
| R-3c      |            0.974 |              0.534 |         0.971 |        0.911 |                       0.731 |                    0.866 |                   0.617 |
| P63/mmc   |            0.957 |              0.704 |         0.98  |        0.997 |                       0.957 |                    0.977 |                   0.997 |
| Ia-3      |            0.992 |              0.914 |         0.99  |        0.992 |                       0.979 |                    0.971 |                   0.977 |
| Fd-3m     |            0.338 |              0.178 |         0.326 |        0.504 |                     nan     |                  nan     |                 nan     |

## bandgap: LOSO median width

| heldout   |   cvplus_lattice |   cvplus_random_in |   ncvplus_knn |   cqr_cvplus |   cvplus_lattice_noFd3m_cal |   ncvplus_knn_noFd3m_cal |   cqr_cvplus_noFd3m_cal |
|:----------|-----------------:|-------------------:|--------------:|-------------:|----------------------------:|-------------------------:|------------------------:|
| C2/m      |            1.901 |              0.747 |         1.998 |        2.565 |                       1.437 |                    1.531 |                   2.254 |
| Pna21     |            2.277 |              0.828 |         2.297 |        2.212 |                       2.235 |                    2.288 |                   1.755 |
| R-3c      |            1.524 |              0.846 |         1.946 |        2.406 |                       1.181 |                    1.535 |                   2.057 |
| P63/mmc   |            1.872 |              0.695 |         2.17  |        2.532 |                       1.502 |                    1.783 |                   2.23  |
| Ia-3      |            1.836 |              0.842 |         1.551 |        2.533 |                       1.56  |                    1.243 |                   1.989 |
| Fd-3m     |            0.864 |              0.376 |         0.846 |        1.622 |                     nan     |                  nan     |                 nan     |

## bandgap: LOSO MAE (fold-model ensembles)

| heldout   |   mae_lattice_ensemble |   mae_random_ensemble |
|:----------|-----------------------:|----------------------:|
| C2/m      |                 0.3782 |                0.3538 |
| Pna21     |                 0.1183 |                0.1333 |
| R-3c      |                 0.3938 |                0.4172 |
| P63/mmc   |                 0.2559 |                0.286  |
| Ia-3      |                 0.2612 |                0.2832 |
| Fd-3m     |                 0.6501 |                0.6251 |

## bandgap: coverage with conflicting records excluded from evaluation

| heldout   |   cvplus_lattice |   cvplus_random_in |   ncvplus_knn |   cqr_cvplus |   cvplus_lattice_noFd3m_cal |   ncvplus_knn_noFd3m_cal |   cqr_cvplus_noFd3m_cal |
|:----------|-----------------:|-------------------:|--------------:|-------------:|----------------------------:|-------------------------:|------------------------:|
| C2/m      |            0.997 |              0.484 |         0.997 |        0.997 |                       0.875 |                    0.896 |                   0.941 |
| Pna21     |            1     |              0.988 |         1     |        1     |                       1     |                    1     |                   0.993 |
| R-3c      |            0.974 |              0.537 |         0.971 |        0.911 |                       0.733 |                    0.865 |                   0.621 |
| P63/mmc   |            0.957 |              0.704 |         0.98  |        0.997 |                       0.957 |                    0.977 |                   0.997 |
| Ia-3      |            0.994 |              0.918 |         0.992 |        0.992 |                       0.981 |                    0.973 |                   0.977 |
| Fd-3m     |            0.334 |              0.175 |         0.321 |        0.504 |                     nan     |                  nan     |                 nan     |

## bandgap: leave-two-out (Fd-3m removed from training and calibration)

| heldout   |   ('coverage', 'cvplus_lattice') |   ('coverage', 'ncvplus_knn') |   ('width_median', 'cvplus_lattice') |   ('width_median', 'ncvplus_knn') |
|:----------|---------------------------------:|------------------------------:|-------------------------------------:|----------------------------------:|
| C2/m      |                            0.994 |                         0.997 |                                1.144 |                             1.247 |
| Pna21     |                            1     |                         1     |                                1.237 |                             1.146 |
| R-3c      |                            0.954 |                         0.963 |                                1.097 |                             1.316 |
| P63/mmc   |                            0.33  |                         0.35  |                                0.528 |                             0.596 |
| Ia-3      |                            0.987 |                         0.979 |                                1.165 |                             1.013 |

## bandgap: calibration composition, share of residuals >= pooled quantile / share of points

| heldout   |   C2/m |   Pna21 |   R-3c |   P63/mmc |   Ia-3 |   Fd-3m |
|:----------|-------:|--------:|-------:|----------:|-------:|--------:|
| C2/m      | nan    |    0.63 |   0.2  |      0.6  |   0.13 |    3.54 |
| Pna21     |   3.45 |  nan    |   0.09 |      0.23 |   0.02 |    1.56 |
| R-3c      |   0.58 |    0.05 | nan    |      0.6  |   0.36 |    3.54 |
| P63/mmc   |   0.2  |    0.07 |   0.37 |    nan    |   0.75 |    3.57 |
| Ia-3      |   0.06 |    0    |   1.12 |      1    | nan    |    2.81 |
| Fd-3m     |   0.35 |    0.02 |   0.57 |      4.39 |   0.17 |  nan    |

## bandgap: group-level conformal (K = 5 families, illustrative; not in the tables above)

| heldout   |   alpha |   k_needed |   threshold |   coverage_of_heldout | family_event_covered   |
|:----------|--------:|-----------:|------------:|----------------------:|:-----------------------|
| C2/m      |    0.05 |          6 |     inf     |                 1     | True                   |
| C2/m      |    0.1  |          6 |     inf     |                 1     | True                   |
| C2/m      |    0.2  |          5 |       1.263 |                 1     | True                   |
| R-3c      |    0.05 |          6 |     inf     |                 1     | True                   |
| R-3c      |    0.1  |          6 |     inf     |                 1     | True                   |
| R-3c      |    0.2  |          5 |       1.096 |                 0.98  | True                   |
| P63/mmc   |    0.05 |          6 |     inf     |                 1     | True                   |
| P63/mmc   |    0.1  |          6 |     inf     |                 1     | True                   |
| P63/mmc   |    0.2  |          5 |       1.373 |                 0.997 | True                   |
| Ia-3      |    0.05 |          6 |     inf     |                 1     | True                   |
| Ia-3      |    0.1  |          6 |     inf     |                 1     | True                   |
| Ia-3      |    0.2  |          5 |       1.241 |                 0.996 | True                   |
| Fd-3m     |    0.05 |          6 |     inf     |                 1     | True                   |
| Fd-3m     |    0.1  |          6 |     inf     |                 1     | True                   |
| Fd-3m     |    0.2  |          5 |       0.787 |                 0.664 | False                  |
| Pna21     |    0.05 |          6 |     inf     |                 1     | True                   |
| Pna21     |    0.1  |          6 |     inf     |                 1     | True                   |
| Pna21     |    0.2  |          5 |       1.268 |                 1     | True                   |

## formation: LOSO coverage of the held-out lattice

| heldout   |   cvplus_lattice |   cvplus_random_in |   ncvplus_knn |   cqr_cvplus |   cvplus_lattice_noFd3m_cal |   ncvplus_knn_noFd3m_cal |   cqr_cvplus_noFd3m_cal |
|:----------|-----------------:|-------------------:|--------------:|-------------:|----------------------------:|-------------------------:|------------------------:|
| C2/m      |            0.997 |              0.248 |         1     |        0.939 |                       0.758 |                    0.723 |                   0.781 |
| Pna21     |            1     |              0.834 |         1     |        1     |                       1     |                    1     |                   1     |
| R-3c      |            1     |              0.44  |         1     |        1     |                       0.997 |                    0.997 |                   0.977 |
| P63/mmc   |            0.915 |              0.288 |         0.926 |        0.877 |                       0.821 |                    0.846 |                   0.758 |
| Ia-3      |            0.994 |              0.714 |         0.979 |        0.99  |                       0.939 |                    0.91  |                   0.919 |
| Fd-3m     |            0.542 |              0.117 |         0.517 |        0.562 |                     nan     |                  nan     |                 nan     |

## formation: LOSO median width

| heldout   |   cvplus_lattice |   cvplus_random_in |   ncvplus_knn |   cqr_cvplus |   cvplus_lattice_noFd3m_cal |   ncvplus_knn_noFd3m_cal |   cqr_cvplus_noFd3m_cal |
|:----------|-----------------:|-------------------:|--------------:|-------------:|----------------------------:|-------------------------:|------------------------:|
| C2/m      |            0.389 |              0.144 |         0.39  |        0.386 |                       0.304 |                    0.299 |                   0.364 |
| Pna21     |            0.417 |              0.117 |         0.416 |        0.366 |                       0.406 |                    0.391 |                   0.329 |
| R-3c      |            0.357 |              0.11  |         0.43  |        0.337 |                       0.314 |                    0.365 |                   0.307 |
| P63/mmc   |            0.355 |              0.099 |         0.376 |        0.339 |                       0.302 |                    0.32  |                   0.254 |
| Ia-3      |            0.355 |              0.122 |         0.298 |        0.336 |                       0.325 |                    0.265 |                   0.29  |
| Fd-3m     |            0.259 |              0.077 |         0.245 |        0.294 |                     nan     |                  nan     |                 nan     |

## formation: LOSO MAE (fold-model ensembles)

| heldout   |   mae_lattice_ensemble |   mae_random_ensemble |
|:----------|-----------------------:|----------------------:|
| C2/m      |                 0.1004 |                0.1008 |
| Pna21     |                 0.0364 |                0.0345 |
| R-3c      |                 0.057  |                0.0598 |
| P63/mmc   |                 0.096  |                0.0886 |
| Ia-3      |                 0.0422 |                0.0508 |
| Fd-3m     |                 0.1193 |                0.1315 |

## formation: coverage with conflicting records excluded from evaluation

| heldout   |   cvplus_lattice |   cvplus_random_in |   ncvplus_knn |   cqr_cvplus |   cvplus_lattice_noFd3m_cal |   ncvplus_knn_noFd3m_cal |   cqr_cvplus_noFd3m_cal |
|:----------|-----------------:|-------------------:|--------------:|-------------:|----------------------------:|-------------------------:|------------------------:|
| C2/m      |            0.997 |              0.234 |         1     |        0.938 |                       0.754 |                    0.718 |                   0.777 |
| Pna21     |            1     |              0.834 |         1     |        1     |                       1     |                    1     |                   1     |
| R-3c      |            1     |              0.443 |         1     |        1     |                       0.997 |                    0.997 |                   0.977 |
| P63/mmc   |            0.915 |              0.288 |         0.926 |        0.877 |                       0.821 |                    0.846 |                   0.758 |
| Ia-3      |            0.998 |              0.716 |         0.983 |        0.994 |                       0.943 |                    0.914 |                   0.922 |
| Fd-3m     |            0.54  |              0.118 |         0.514 |        0.56  |                     nan     |                  nan     |                 nan     |

## formation: leave-two-out (Fd-3m removed from training and calibration)

| heldout   |   ('coverage', 'cvplus_lattice') |   ('coverage', 'ncvplus_knn') |   ('width_median', 'cvplus_lattice') |   ('width_median', 'ncvplus_knn') |
|:----------|---------------------------------:|------------------------------:|-------------------------------------:|----------------------------------:|
| C2/m      |                            0.997 |                         0.997 |                                0.291 |                             0.325 |
| Pna21     |                            1     |                         1     |                                0.301 |                             0.275 |
| R-3c      |                            1     |                         1     |                                0.288 |                             0.37  |
| P63/mmc   |                            0.299 |                         0.35  |                                0.165 |                             0.171 |
| Ia-3      |                            0.996 |                         0.996 |                                0.364 |                             0.279 |

## formation: calibration composition, share of residuals >= pooled quantile / share of points

| heldout   |   C2/m |   Pna21 |   R-3c |   P63/mmc |   Ia-3 |   Fd-3m |
|:----------|-------:|--------:|-------:|----------:|-------:|--------:|
| C2/m      | nan    |    0.23 |   0    |      1    |   0.79 |    2.98 |
| Pna21     |   3.28 |  nan    |   0    |      0.51 |   0.06 |    1.48 |
| R-3c      |   1.58 |    0    | nan    |      0.63 |   0.38 |    2.68 |
| P63/mmc   |   2.1  |    0    |   0    |    nan    |   0.06 |    3.16 |
| Ia-3      |   1.84 |    0    |   0    |      1.14 | nan    |    2.12 |
| Fd-3m     |   0.17 |    0.02 |   0.11 |      4.13 |   0.81 |  nan    |

## formation: group-level conformal (K = 5 families, illustrative; not in the tables above)

| heldout   |   alpha |   k_needed |   threshold |   coverage_of_heldout | family_event_covered   |
|:----------|--------:|-----------:|------------:|----------------------:|:-----------------------|
| C2/m      |    0.05 |          6 |     inf     |                 1     | True                   |
| C2/m      |    0.1  |          6 |     inf     |                 1     | True                   |
| C2/m      |    0.2  |          5 |       0.217 |                 0.997 | True                   |
| R-3c      |    0.05 |          6 |     inf     |                 1     | True                   |
| R-3c      |    0.1  |          6 |     inf     |                 1     | True                   |
| R-3c      |    0.2  |          5 |       0.225 |                 1     | True                   |
| P63/mmc   |    0.05 |          6 |     inf     |                 1     | True                   |
| P63/mmc   |    0.1  |          6 |     inf     |                 1     | True                   |
| P63/mmc   |    0.2  |          5 |       0.225 |                 0.991 | True                   |
| Ia-3      |    0.05 |          6 |     inf     |                 1     | True                   |
| Ia-3      |    0.1  |          6 |     inf     |                 1     | True                   |
| Ia-3      |    0.2  |          5 |       0.211 |                 0.996 | True                   |
| Fd-3m     |    0.05 |          6 |     inf     |                 1     | True                   |
| Fd-3m     |    0.1  |          6 |     inf     |                 1     | True                   |
| Fd-3m     |    0.2  |          5 |       0.174 |                 0.791 | False                  |
| Pna21     |    0.05 |          6 |     inf     |                 1     | True                   |
| Pna21     |    0.1  |          6 |     inf     |                 1     | True                   |
| Pna21     |    0.2  |          5 |       0.231 |                 1     | True                   |

## Pre-registered decisions

```
{
  "bandgap": {
    "P2.1_fd3m_ratio_by_holdout": {
      "C2/m": 3.537,
      "R-3c": 3.542,
      "P63/mmc": 3.566,
      "Ia-3": 2.805,
      "Pna21": 1.557
    },
    "P2.1_supported": true,
    "H2_ablation": {
      "delta_coverage": {
        "C2/m": -0.12,
        "Pna21": 0.0,
        "R-3c": -0.243,
        "P63/mmc": 0.0,
        "Ia-3": -0.013
      },
      "width_ratio": {
        "C2/m": 0.756,
        "Pna21": 0.981,
        "R-3c": 0.775,
        "P63/mmc": 0.802,
        "Ia-3": 0.85
      },
      "n_support": 2,
      "n_reject": 1,
      "decision": "partial / mixed"
    },
    "H2_leave_two_out": {
      "delta_coverage": {
        "C2/m": -0.003,
        "Pna21": 0.0,
        "R-3c": -0.02,
        "P63/mmc": -0.627,
        "Ia-3": -0.004
      },
      "width_ratio": {
        "C2/m": 0.601,
        "Pna21": 0.543,
        "R-3c": 0.719,
        "P63/mmc": 0.282,
        "Ia-3": 0.634
      },
      "n_support": 1,
      "n_reject": 0,
      "decision": "partial / mixed"
    },
    "G3_fd3m_coverage": {
      "cvplus_lattice": 0.3384223918575064,
      "ncvplus_knn": 0.3256997455470738,
      "cqr_cvplus": 0.5038167938931297
    },
    "G3_fd3m_width_median": {
      "cvplus_lattice": 0.8641113545367909,
      "ncvplus_knn": 0.8462155721959332,
      "cqr_cvplus": 1.6217314448707254
    },
    "G3_positive": false,
    "T2.3_subfamily_coverage": {
      "Fd-3m (gamma~30)": 0.3441295546558704,
      "Fd-3m (gamma~120)": 0.3287671232876712
    },
    "T2.3_one_subfamily_carries_failure": false
  },
  "formation": {
    "P2.1_fd3m_ratio_by_holdout": {
      "C2/m": 2.977,
      "R-3c": 2.676,
      "P63/mmc": 3.158,
      "Ia-3": 2.117,
      "Pna21": 1.48
    },
    "P2.1_supported": true,
    "H2_ablation": {
      "delta_coverage": {
        "C2/m": -0.239,
        "Pna21": 0.0,
        "R-3c": -0.003,
        "P63/mmc": -0.094,
        "Ia-3": -0.054
      },
      "width_ratio": {
        "C2/m": 0.781,
        "Pna21": 0.975,
        "R-3c": 0.879,
        "P63/mmc": 0.85,
        "Ia-3": 0.915
      },
      "n_support": 1,
      "n_reject": 1,
      "decision": "partial / mixed"
    },
    "H2_leave_two_out": {
      "delta_coverage": {
        "C2/m": 0.0,
        "Pna21": 0.0,
        "R-3c": 0.0,
        "P63/mmc": -0.615,
        "Ia-3": 0.002
      },
      "width_ratio": {
        "C2/m": 0.749,
        "Pna21": 0.723,
        "R-3c": 0.808,
        "P63/mmc": 0.464,
        "Ia-3": 1.024
      },
      "n_support": 1,
      "n_reject": 1,
      "decision": "partial / mixed"
    },
    "G3_fd3m_coverage": {
      "cvplus_lattice": 0.5419847328244275,
      "ncvplus_knn": 0.5165394402035624,
      "cqr_cvplus": 0.5623409669211196
    },
    "G3_fd3m_width_median": {
      "cvplus_lattice": 0.2585988264009746,
      "ncvplus_knn": 0.2451484330214237,
      "cqr_cvplus": 0.2942828219672008
    },
    "G3_positive": false,
    "T2.3_subfamily_coverage": {
      "Fd-3m (gamma~30)": 0.6153846153846154,
      "Fd-3m (gamma~120)": 0.4178082191780822
    },
    "T2.3_one_subfamily_carries_failure": false
  }
}
```
