# Phase 1 follow-up checks (exploratory; alpha = 0.1)

Phase 1 files re-analysed: 50. Stored intervals reproduced exactly: yes.
v2 de-duplication removes 57 entries (Phase 1 removed 9); 16 conflicting records flagged.

## bandgap: mean coverage, Phase 1 data vs v2 de-duplicated (models unchanged)

| group   |     cqr |   cvplus_lattice |   cvplus_random |   split_global |   split_mondrian |   split_normalized_knn |
|:--------|--------:|-----------------:|----------------:|---------------:|-----------------:|-----------------------:|
| ALL     | -0.0002 |          -0.0002 |          0.0007 |        -0.0008 |           0.0006 |                -0.0011 |
| C2/m    | -0.0003 |           0      |          0      |         0.0016 |           0.011  |                 0.0008 |
| Pna21   | -0.0015 |           0      |         -0      |        -0.0004 |           0.0002 |                 0.0002 |
| R-3c    |  0.0135 |           0.0021 |          0.0202 |         0.02   |          -0.0011 |                 0.0199 |
| P63/mmc | -0.0023 |          -0.0001 |         -0.0006 |        -0.0017 |          -0.0007 |                -0.0045 |
| Ia-3    | -0.0055 |          -0.0001 |         -0.0003 |        -0.0004 |           0.0004 |                -0.0016 |
| Fd-3m   | -0.0006 |          -0.0008 |         -0.0008 |        -0.0081 |           0      |                -0.0072 |

Coverage on v2 data with conflicting records excluded from evaluation:

| group   |   cqr |   cvplus_lattice |   cvplus_random |   split_global |   split_mondrian |   split_normalized_knn |
|:--------|------:|-----------------:|----------------:|---------------:|-----------------:|-----------------------:|
| ALL     | 0.899 |            0.985 |           0.909 |          0.888 |            0.899 |                  0.892 |
| C2/m    | 0.886 |            1     |           0.985 |          0.968 |            0.908 |                  0.975 |
| Pna21   | 0.93  |            1     |           0.998 |          0.995 |            0.913 |                  0.991 |
| R-3c    | 0.881 |            0.999 |           0.984 |          0.969 |            0.889 |                  0.972 |
| P63/mmc | 0.899 |            0.988 |           0.904 |          0.866 |            0.902 |                  0.856 |
| Ia-3    | 0.932 |            0.999 |           0.985 |          0.976 |            0.888 |                  0.972 |
| Fd-3m   | 0.857 |            0.931 |           0.598 |          0.533 |            0.899 |                  0.569 |

## bandgap: Mondrian, pooled vs one-point-per-group calibration (Dunn et al. 2023)

| lattice   |   n_cal_groups |   pooled_cov_point |   pooled_cov_group |   single_sub_cov_group |   single_sub_frac_infinite |   pooled_width |   single_sub_width_median |
|:----------|---------------:|-------------------:|-------------------:|-----------------------:|---------------------------:|---------------:|--------------------------:|
| C2/m      |          11.76 |              0.901 |              0.874 |                  0.939 |                       0.08 |          0.404 |                   inf     |
| Pna21     |          13    |              0.913 |              0.905 |                  0.926 |                       0    |          0.318 |                     0.366 |
| R-3c      |          13.04 |              0.889 |              0.895 |                  0.926 |                       0    |          0.376 |                     0.666 |
| P63/mmc   |          13.72 |              0.902 |              0.881 |                  0.915 |                       0    |          0.732 |                     0.942 |
| Ia-3      |          12.44 |              0.885 |              0.868 |                  0.911 |                       0    |          0.312 |                     0.501 |
| Fd-3m     |           9.36 |              0.898 |              0.927 |                  0.943 |                       0.36 |          1.372 |                   inf     |

## formation: mean coverage, Phase 1 data vs v2 de-duplicated (models unchanged)

| group   |     cqr |   cvplus_lattice |   cvplus_random |   split_global |   split_mondrian |   split_normalized_knn |
|:--------|--------:|-----------------:|----------------:|---------------:|-----------------:|-----------------------:|
| ALL     | -0.0011 |          -0.0001 |         -0.0004 |         0.0017 |           0.0007 |                 0.0025 |
| C2/m    |  0.0006 |           0      |         -0      |         0.0016 |           0.001  |                 0.0041 |
| Pna21   | -0.0096 |           0      |          0.0003 |         0.0015 |          -0.0005 |                 0.0004 |
| R-3c    |  0.0131 |           0      |          0.0003 |         0.0013 |           0.006  |                -0.0024 |
| P63/mmc | -0.0032 |          -0.0001 |         -0.0006 |         0.002  |           0.0011 |                 0.0027 |
| Ia-3    | -0.0017 |          -0      |          0.0013 |         0.0056 |           0.0002 |                 0.0069 |
| Fd-3m   | -0.0024 |           0      |          0.0044 |         0.0048 |           0      |                 0.0107 |

Coverage on v2 data with conflicting records excluded from evaluation:

| group   |   cqr |   cvplus_lattice |   cvplus_random |   split_global |   split_mondrian |   split_normalized_knn |
|:--------|------:|-----------------:|----------------:|---------------:|-----------------:|-----------------------:|
| ALL     | 0.903 |            0.997 |           0.92  |          0.899 |            0.897 |                  0.896 |
| C2/m    | 0.88  |            1     |           0.972 |          0.944 |            0.888 |                  0.935 |
| Pna21   | 0.917 |            1     |           0.99  |          0.98  |            0.911 |                  0.967 |
| R-3c    | 0.917 |            1     |           0.975 |          0.953 |            0.907 |                  0.984 |
| P63/mmc | 0.895 |            0.993 |           0.941 |          0.926 |            0.903 |                  0.919 |
| Ia-3    | 0.931 |            1     |           0.996 |          0.99  |            0.856 |                  0.975 |
| Fd-3m   | 0.869 |            0.989 |           0.653 |          0.588 |            0.907 |                  0.598 |

## formation: Mondrian, pooled vs one-point-per-group calibration (Dunn et al. 2023)

| lattice   |   n_cal_groups |   pooled_cov_point |   pooled_cov_group |   single_sub_cov_group |   single_sub_frac_infinite |   pooled_width |   single_sub_width_median |
|:----------|---------------:|-------------------:|-------------------:|-----------------------:|---------------------------:|---------------:|--------------------------:|
| C2/m      |          11.76 |              0.882 |              0.855 |                  0.947 |                       0.08 |          0.073 |                   inf     |
| Pna21     |          13    |              0.911 |              0.91  |                  0.928 |                       0    |          0.058 |                     0.07  |
| R-3c      |          13.04 |              0.908 |              0.885 |                  0.932 |                       0    |          0.076 |                     0.103 |
| P63/mmc   |          13.72 |              0.903 |              0.894 |                  0.932 |                       0    |          0.088 |                     0.116 |
| Ia-3      |          12.44 |              0.854 |              0.871 |                  0.91  |                       0    |          0.04  |                     0.047 |
| Fd-3m     |           9.36 |              0.906 |              0.9   |                  0.945 |                       0.36 |          0.187 |                   inf     |

## Phase 1 by Fd-3m sub-family (cell setting)

| target    | method               | group             |   cov_mean |   width_median |   mae |   n_test_mean |
|:----------|:---------------------|:------------------|-----------:|---------------:|------:|--------------:|
| bandgap   | split_global         | Fd-3m             |      0.543 |          0.626 | 0.329 |         79.8  |
| bandgap   | split_global         | Fd-3m (gamma~30)  |      0.542 |          0.626 | 0.325 |         49.56 |
| bandgap   | split_global         | Fd-3m (gamma~120) |      0.552 |          0.626 | 0.333 |         30.24 |
| bandgap   | split_mondrian       | Fd-3m             |      0.898 |          1.361 | 0.329 |         79.8  |
| bandgap   | split_mondrian       | Fd-3m (gamma~30)  |      0.903 |          1.361 | 0.325 |         49.56 |
| bandgap   | split_mondrian       | Fd-3m (gamma~120) |      0.887 |          1.361 | 0.333 |         30.24 |
| bandgap   | split_normalized_knn | Fd-3m             |      0.577 |          0.649 | 0.329 |         79.8  |
| bandgap   | split_normalized_knn | Fd-3m (gamma~30)  |      0.575 |          0.658 | 0.325 |         49.56 |
| bandgap   | split_normalized_knn | Fd-3m (gamma~120) |      0.588 |          0.631 | 0.333 |         30.24 |
| bandgap   | cvplus_random        | Fd-3m             |      0.6   |          0.688 | 0.321 |         79.8  |
| bandgap   | cvplus_random        | Fd-3m (gamma~30)  |      0.607 |          0.681 | 0.318 |         49.56 |
| bandgap   | cvplus_random        | Fd-3m (gamma~120) |      0.596 |          0.685 | 0.327 |         30.24 |
| bandgap   | cvplus_lattice       | Fd-3m             |      0.932 |          1.654 | 0.324 |         79.8  |
| bandgap   | cvplus_lattice       | Fd-3m (gamma~30)  |      0.932 |          1.659 | 0.318 |         49.56 |
| bandgap   | cvplus_lattice       | Fd-3m (gamma~120) |      0.933 |          1.633 | 0.334 |         30.24 |
| bandgap   | cqr                  | Fd-3m             |      0.859 |          1.369 | 0.329 |         79.8  |
| bandgap   | cqr                  | Fd-3m (gamma~30)  |      0.894 |          1.35  | 0.325 |         49.56 |
| bandgap   | cqr                  | Fd-3m (gamma~120) |      0.8   |          1.388 | 0.333 |         30.24 |
| formation | split_global         | Fd-3m             |      0.582 |          0.091 | 0.045 |         79.8  |
| formation | split_global         | Fd-3m (gamma~30)  |      0.557 |          0.091 | 0.045 |         49.56 |
| formation | split_global         | Fd-3m (gamma~120) |      0.616 |          0.091 | 0.046 |         30.24 |
| formation | split_mondrian       | Fd-3m             |      0.906 |          0.176 | 0.045 |         79.8  |
| formation | split_mondrian       | Fd-3m (gamma~30)  |      0.915 |          0.176 | 0.045 |         49.56 |
| formation | split_mondrian       | Fd-3m (gamma~120) |      0.889 |          0.176 | 0.046 |         30.24 |
| formation | split_normalized_knn | Fd-3m             |      0.586 |          0.091 | 0.045 |         79.8  |
| formation | split_normalized_knn | Fd-3m (gamma~30)  |      0.556 |          0.09  | 0.045 |         49.56 |
| formation | split_normalized_knn | Fd-3m (gamma~120) |      0.627 |          0.089 | 0.046 |         30.24 |
| formation | cvplus_random        | Fd-3m             |      0.646 |          0.1   | 0.044 |         79.8  |
| formation | cvplus_random        | Fd-3m (gamma~30)  |      0.618 |          0.099 | 0.044 |         49.56 |
| formation | cvplus_random        | Fd-3m (gamma~120) |      0.688 |          0.1   | 0.044 |         30.24 |
| formation | cvplus_lattice       | Fd-3m             |      0.989 |          0.365 | 0.047 |         79.8  |
| formation | cvplus_lattice       | Fd-3m (gamma~30)  |      1     |          0.355 | 0.046 |         49.56 |
| formation | cvplus_lattice       | Fd-3m (gamma~120) |      0.969 |          0.388 | 0.048 |         30.24 |
| formation | cqr                  | Fd-3m             |      0.87  |          0.233 | 0.045 |         79.8  |
| formation | cqr                  | Fd-3m (gamma~30)  |      0.909 |          0.224 | 0.045 |         49.56 |
| formation | cqr                  | Fd-3m (gamma~120) |      0.796 |          0.25  | 0.046 |         30.24 |
