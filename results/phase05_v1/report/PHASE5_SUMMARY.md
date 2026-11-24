# Phase 5 summary (alpha = 0.1; 500 draws per cell)

## D5.0 code check (random scheme): PASSED (0 of 168 cells outside 4 SE)

## H5 / H5c dispersion vs permutation null (split, m = 20; coverage of new compositions)

| scheme        | heldout   | design   | target    |   mean_real |   mean_null |   nominal |   sd_real |   sd_null |   rho |   rho_lo |   rho_hi |   p_below_real |   p_below_null |
|:--------------|:----------|:---------|:----------|------------:|------------:|----------:|----------:|----------:|------:|---------:|---------:|---------------:|---------------:|
| clustered     | P63/mmc   | l2o      | bandgap   |       0.879 |       0.904 |     0.905 |     0.11  |     0.064 | 1.717 |    1.538 |    1.918 |          0.204 |          0.078 |
| clustered     | P63/mmc   | l2o      | formation |       0.898 |       0.903 |     0.905 |     0.072 |     0.068 | 1.061 |    0.93  |    1.209 |          0.11  |          0.074 |
| clustered     | C2/m      | loso     | bandgap   |       0.841 |       0.907 |     0.905 |     0.108 |     0.068 | 1.604 |    1.43  |    1.792 |          0.306 |          0.082 |
| clustered     | C2/m      | loso     | formation |       0.866 |       0.909 |     0.905 |     0.1   |     0.064 | 1.561 |    1.408 |    1.734 |          0.25  |          0.062 |
| clustered     | Pna21     | loso     | bandgap   |       0.849 |       0.905 |     0.905 |     0.109 |     0.064 | 1.695 |    1.528 |    1.874 |          0.288 |          0.078 |
| clustered     | Pna21     | loso     | formation |       0.883 |       0.901 |     0.905 |     0.112 |     0.066 | 1.692 |    1.513 |    1.881 |          0.174 |          0.074 |
| clustered     | R-3c      | loso     | bandgap   |       0.847 |       0.905 |     0.905 |     0.147 |     0.066 | 2.228 |    1.97  |    2.527 |          0.268 |          0.086 |
| clustered     | R-3c      | loso     | formation |       0.85  |       0.904 |     0.905 |     0.119 |     0.066 | 1.797 |    1.629 |    1.982 |          0.358 |          0.076 |
| clustered     | P63/mmc   | loso     | bandgap   |       0.879 |       0.908 |     0.905 |     0.103 |     0.063 | 1.634 |    1.44  |    1.867 |          0.202 |          0.058 |
| clustered     | P63/mmc   | loso     | formation |       0.897 |       0.904 |     0.905 |     0.074 |     0.07  | 1.05  |    0.919 |    1.197 |          0.096 |          0.086 |
| clustered     | Ia-3      | loso     | bandgap   |       0.809 |       0.905 |     0.905 |     0.176 |     0.065 | 2.72  |    2.381 |    3.07  |          0.362 |          0.066 |
| clustered     | Ia-3      | loso     | formation |       0.761 |       0.901 |     0.905 |     0.216 |     0.072 | 3.001 |    2.713 |    3.343 |          0.462 |          0.092 |
| clustered     | Fd-3m     | loso     | bandgap   |       0.898 |       0.906 |     0.905 |     0.071 |     0.066 | 1.074 |    0.974 |    1.177 |          0.116 |          0.084 |
| clustered     | Fd-3m     | loso     | formation |       0.915 |       0.912 |     0.905 |     0.075 |     0.062 | 1.198 |    1.049 |    1.343 |          0.09  |          0.064 |
| one_per_group | P63/mmc   | l2o      | bandgap   |       0.906 |       0.903 |     0.905 |     0.067 |     0.069 | 0.975 |    0.872 |    1.099 |          0.08  |          0.078 |
| one_per_group | P63/mmc   | l2o      | formation |       0.906 |       0.906 |     0.905 |     0.066 |     0.062 | 1.079 |    0.965 |    1.21  |          0.078 |          0.058 |
| one_per_group | C2/m      | loso     | bandgap   |       0.901 |       0.905 |     0.905 |     0.075 |     0.063 | 1.206 |    1.074 |    1.343 |          0.086 |          0.074 |
| one_per_group | C2/m      | loso     | formation |       0.909 |       0.911 |     0.905 |     0.071 |     0.066 | 1.068 |    0.951 |    1.183 |          0.082 |          0.062 |
| one_per_group | Pna21     | loso     | bandgap   |       0.907 |       0.902 |     0.905 |     0.068 |     0.068 | 1.009 |    0.898 |    1.138 |          0.08  |          0.074 |
| one_per_group | Pna21     | loso     | formation |       0.908 |       0.909 |     0.905 |     0.064 |     0.066 | 0.967 |    0.871 |    1.073 |          0.072 |          0.076 |
| one_per_group | R-3c      | loso     | bandgap   |       0.899 |       0.901 |     0.905 |     0.074 |     0.069 | 1.083 |    0.981 |    1.194 |          0.102 |          0.088 |
| one_per_group | R-3c      | loso     | formation |       0.907 |       0.906 |     0.905 |     0.071 |     0.065 | 1.079 |    0.963 |    1.206 |          0.094 |          0.068 |
| one_per_group | P63/mmc   | loso     | bandgap   |       0.909 |       0.901 |     0.905 |     0.068 |     0.071 | 0.954 |    0.823 |    1.105 |          0.07  |          0.082 |
| one_per_group | P63/mmc   | loso     | formation |       0.901 |       0.905 |     0.905 |     0.069 |     0.066 | 1.048 |    0.913 |    1.195 |          0.068 |          0.08  |
| one_per_group | Ia-3      | loso     | bandgap   |       0.905 |       0.91  |     0.905 |     0.074 |     0.066 | 1.126 |    1     |    1.274 |          0.09  |          0.062 |
| one_per_group | Ia-3      | loso     | formation |       0.906 |       0.904 |     0.905 |     0.071 |     0.066 | 1.074 |    0.974 |    1.195 |          0.072 |          0.08  |
| one_per_group | Fd-3m     | loso     | bandgap   |       0.899 |       0.909 |     0.905 |     0.075 |     0.068 | 1.114 |    0.988 |    1.261 |          0.096 |          0.074 |
| one_per_group | Fd-3m     | loso     | formation |       0.905 |       0.909 |     0.905 |     0.069 |     0.063 | 1.09  |    0.971 |    1.22  |          0.086 |          0.064 |

```
{
  "H5_bandgap": {
    "verdict": "supported",
    "n_families": 7,
    "n_rho_le_support": 7,
    "n_rho_gt_reject": 0,
    "rho": {
      "P63/mmc (l2o)": 0.975,
      "C2/m (loso)": 1.206,
      "Pna21 (loso)": 1.009,
      "R-3c (loso)": 1.083,
      "P63/mmc (loso)": 0.954,
      "Ia-3 (loso)": 1.126,
      "Fd-3m (loso)": 1.114
    }
  },
  "H5_formation": {
    "verdict": "supported",
    "n_families": 7,
    "n_rho_le_support": 7,
    "n_rho_gt_reject": 0,
    "rho": {
      "P63/mmc (l2o)": 1.079,
      "C2/m (loso)": 1.068,
      "Pna21 (loso)": 0.967,
      "R-3c (loso)": 1.079,
      "P63/mmc (loso)": 1.048,
      "Ia-3 (loso)": 1.074,
      "Fd-3m (loso)": 1.09
    }
  },
  "H5c_bandgap": {
    "verdict": "supported",
    "n_families": 7,
    "n_harmed": 6,
    "rho": {
      "P63/mmc (l2o)": 1.717,
      "C2/m (loso)": 1.604,
      "Pna21 (loso)": 1.695,
      "R-3c (loso)": 2.228,
      "P63/mmc (loso)": 1.634,
      "Ia-3 (loso)": 2.72,
      "Fd-3m (loso)": 1.074
    },
    "mean_cov_new": {
      "P63/mmc (l2o)": 0.879,
      "C2/m (loso)": 0.841,
      "Pna21 (loso)": 0.849,
      "R-3c (loso)": 0.847,
      "P63/mmc (loso)": 0.879,
      "Ia-3 (loso)": 0.809,
      "Fd-3m (loso)": 0.898
    }
  },
  "H5c_formation": {
    "verdict": "supported",
    "n_families": 7,
    "n_harmed": 4,
    "rho": {
      "P63/mmc (l2o)": 1.061,
      "C2/m (loso)": 1.561,
      "Pna21 (loso)": 1.692,
      "R-3c (loso)": 1.797,
      "P63/mmc (loso)": 1.05,
      "Ia-3 (loso)": 3.001,
      "Fd-3m (loso)": 1.198
    },
    "mean_cov_new": {
      "P63/mmc (l2o)": 0.898,
      "C2/m (loso)": 0.866,
      "Pna21 (loso)": 0.883,
      "R-3c (loso)": 0.85,
      "P63/mmc (loso)": 0.897,
      "Ia-3 (loso)": 0.761,
      "Fd-3m (loso)": 0.915
    }
  }
}
```

## H5b offset correction (one_per_group, m = 20)

| target    | design   | heldout   |   m |   width_split |   width_offset |   width_ratio |   cov_split |   cov_offset | passes   |
|:----------|:---------|:----------|----:|--------------:|---------------:|--------------:|------------:|-------------:|:---------|
| bandgap   | loso     | Fd-3m     |  20 |         2.655 |          1.513 |         0.57  |       0.899 |        0.898 | True     |
| bandgap   | loso     | C2/m      |  20 |         1.346 |          0.895 |         0.665 |       0.901 |        0.907 | True     |
| bandgap   | loso     | R-3c      |  20 |         1.297 |          0.869 |         0.669 |       0.899 |        0.902 | True     |
| bandgap   | loso     | Ia-3      |  20 |         0.872 |          0.496 |         0.568 |       0.905 |        0.903 | True     |
| bandgap   | l2o      | P63/mmc   |  20 |         1.924 |          1.191 |         0.619 |       0.906 |        0.907 | True     |
| formation | loso     | Fd-3m     |  20 |         0.409 |          0.256 |         0.627 |       0.905 |        0.905 | True     |
| formation | loso     | C2/m      |  20 |         0.32  |          0.244 |         0.765 |       0.909 |        0.909 | True     |
| formation | loso     | R-3c      |  20 |         0.186 |          0.14  |         0.751 |       0.907 |        0.91  | True     |
| formation | loso     | Ia-3      |  20 |         0.21  |          0.159 |         0.76  |       0.906 |        0.905 | True     |
| formation | loso     | Pna21     |  20 |         0.143 |          0.105 |         0.737 |       0.908 |        0.904 | True     |
| formation | l2o      | P63/mmc   |  20 |         0.38  |          0.299 |         0.787 |       0.906 |        0.903 | True     |

```
{
  "bandgap": {
    "verdict": "supported",
    "n_cells": 5,
    "n_pass": 5,
    "needed": 4,
    "n_ratio_gt_reject": 0
  },
  "formation": {
    "verdict": "supported",
    "n_cells": 6,
    "n_pass": 6,
    "needed": 5,
    "n_ratio_gt_reject": 0
  }
}
```

## Practical endpoint m* (mean >= 0.88, P(cov < 0.8) <= 0.1)

| heldout   | design   | target    |   clustered offset_full |   clustered split |   one_per_group offset_full |   one_per_group split |   random offset_full |   random split |
|:----------|:---------|:----------|------------------------:|------------------:|----------------------------:|----------------------:|---------------------:|---------------:|
| C2/m      | loso     | bandgap   |                      -1 |                -1 |                          20 |                    20 |                   20 |             20 |
| C2/m      | loso     | formation |                      -1 |                50 |                          20 |                    20 |                   20 |             20 |
| Fd-3m     | loso     | bandgap   |                      30 |                30 |                          30 |                    20 |                   10 |             20 |
| Fd-3m     | loso     | formation |                      50 |                20 |                          20 |                     9 |                   20 |             20 |
| Ia-3      | loso     | bandgap   |                      -1 |                -1 |                          20 |                    10 |                   20 |             20 |
| Ia-3      | loso     | formation |                      -1 |                -1 |                          20 |                    20 |                   10 |             10 |
| P63/mmc   | l2o      | bandgap   |                      -1 |                -1 |                          20 |                    20 |                   10 |             20 |
| P63/mmc   | l2o      | formation |                      -1 |                30 |                          10 |                    10 |                   20 |             20 |
| P63/mmc   | loso     | bandgap   |                     100 |                -1 |                          20 |                    20 |                   20 |             20 |
| P63/mmc   | loso     | formation |                      -1 |                20 |                          10 |                    20 |                   20 |             10 |
| Pna21     | loso     | bandgap   |                      -1 |                -1 |                          20 |                    20 |                   20 |             20 |
| Pna21     | loso     | formation |                     100 |                50 |                          20 |                    20 |                   10 |             10 |
| R-3c      | loso     | bandgap   |                      -1 |               100 |                          20 |                    30 |                   20 |             10 |
| R-3c      | loso     | formation |                      -1 |                -1 |                          20 |                    20 |                   10 |             10 |
(-1 = not reached within the m grid)

## Widths at m = 20 (one_per_group), in units of the family SD, next to LOSO CV+

| design   | heldout   | target    |   mean_cov_offset_full |   mean_cov_split |   width_sd_offset_full |   width_sd_split |   loso_coverage |   loso_width_sd |
|:---------|:----------|:----------|-----------------------:|-----------------:|-----------------------:|-----------------:|----------------:|----------------:|
| l2o      | P63/mmc   | bandgap   |                  0.907 |            0.906 |                  1.372 |            2.217 |           0.33  |           0.608 |
| l2o      | P63/mmc   | formation |                  0.903 |            0.906 |                  3.11  |            3.955 |           0.299 |           1.714 |
| loso     | C2/m      | bandgap   |                  0.907 |            0.901 |                  1.405 |            2.112 |           0.997 |           2.984 |
| loso     | C2/m      | formation |                  0.909 |            0.909 |                  2.815 |            3.681 |           0.997 |           4.481 |
| loso     | Fd-3m     | bandgap   |                  0.898 |            0.899 |                  1.886 |            3.309 |           0.338 |           1.077 |
| loso     | Fd-3m     | formation |                  0.905 |            0.905 |                  2.623 |            4.185 |           0.542 |           2.646 |
| loso     | Ia-3      | bandgap   |                  0.903 |            0.905 |                  0.503 |            0.885 |           0.992 |           1.862 |
| loso     | Ia-3      | formation |                  0.905 |            0.906 |                  2.398 |            3.154 |           0.994 |           5.347 |
| loso     | P63/mmc   | bandgap   |                  0.908 |            0.909 |                  1.513 |            1.541 |           0.957 |           2.156 |
| loso     | P63/mmc   | formation |                  0.903 |            0.901 |                  3.613 |            3.619 |           0.915 |           3.697 |
| loso     | Pna21     | bandgap   |                  0.904 |            0.907 |                  0.614 |            0.634 |           1     |           2.395 |
| loso     | Pna21     | formation |                  0.904 |            0.908 |                  1.459 |            1.979 |           1     |           5.768 |
| loso     | R-3c      | bandgap   |                  0.902 |            0.899 |                  0.845 |            1.263 |           0.974 |           1.483 |
| loso     | R-3c      | formation |                  0.91  |            0.907 |                  1.884 |            2.509 |           1     |           4.823 |

## MAE before / after subtracting the estimated offset (m = 20, one_per_group)

| design   | heldout   | target    |   mae_full |   mae_raw |   mae_offset |   bias_full |   b_hat |
|:---------|:----------|:----------|-----------:|----------:|-------------:|------------:|--------:|
| l2o      | P63/mmc   | bandgap   |      0.442 |     0.442 |        0.242 |      -0.427 |  -0.43  |
| loso     | C2/m      | bandgap   |      0.378 |     0.378 |        0.161 |       0.367 |   0.353 |
| loso     | R-3c      | bandgap   |      0.394 |     0.396 |        0.19  |       0.384 |   0.336 |
| loso     | P63/mmc   | bandgap   |      0.256 |     0.255 |        0.257 |      -0.127 |  -0.148 |
| loso     | Ia-3      | bandgap   |      0.261 |     0.26  |        0.123 |       0.251 |   0.273 |
| loso     | Fd-3m     | bandgap   |      0.65  |     0.653 |        0.351 |      -0.643 |  -0.593 |
| loso     | Pna21     | bandgap   |      0.118 |     0.118 |        0.124 |      -0     |  -0.027 |
| l2o      | P63/mmc   | formation |      0.115 |     0.115 |        0.065 |       0.099 |   0.097 |
| loso     | C2/m      | formation |      0.1   |     0.101 |        0.049 |      -0.093 |  -0.081 |
| loso     | R-3c      | formation |      0.057 |     0.058 |        0.035 |      -0.053 |  -0.042 |
| loso     | P63/mmc   | formation |      0.096 |     0.096 |        0.071 |       0.064 |   0.062 |
| loso     | Ia-3      | formation |      0.042 |     0.043 |        0.042 |      -0.036 |  -0.031 |
| loso     | Fd-3m     | formation |      0.119 |     0.12  |        0.061 |       0.117 |   0.109 |
| loso     | Pna21     | formation |      0.036 |     0.037 |        0.026 |      -0.031 |  -0.028 |

## Clustered vs one-per-group dispersion (exploratory)

| design   | heldout   | target    |   m |   clustered |   one_per_group |   sd_ratio_clustered_to_one_per_group |
|:---------|:----------|:----------|----:|------------:|----------------:|--------------------------------------:|
| l2o      | P63/mmc   | bandgap   |   9 |       0.155 |           0.09  |                                 1.719 |
| l2o      | P63/mmc   | bandgap   |  10 |       0.129 |           0.084 |                                 1.53  |
| l2o      | P63/mmc   | bandgap   |  20 |       0.11  |           0.067 |                                 1.632 |
| l2o      | P63/mmc   | bandgap   |  30 |       0.095 |           0.059 |                                 1.591 |
| l2o      | P63/mmc   | formation |   9 |       0.108 |           0.087 |                                 1.241 |
| l2o      | P63/mmc   | formation |  10 |       0.101 |           0.079 |                                 1.286 |
| l2o      | P63/mmc   | formation |  20 |       0.072 |           0.066 |                                 1.09  |
| l2o      | P63/mmc   | formation |  30 |       0.058 |           0.062 |                                 0.923 |
| loso     | C2/m      | bandgap   |   9 |       0.161 |           0.093 |                                 1.726 |
| loso     | C2/m      | bandgap   |  10 |       0.137 |           0.085 |                                 1.618 |
| loso     | C2/m      | bandgap   |  20 |       0.108 |           0.075 |                                 1.436 |
| loso     | C2/m      | bandgap   |  30 |       0.098 |           0.064 |                                 1.539 |
| loso     | C2/m      | formation |   9 |       0.125 |           0.102 |                                 1.216 |
| loso     | C2/m      | formation |  10 |       0.115 |           0.083 |                                 1.381 |
| loso     | C2/m      | formation |  20 |       0.1   |           0.071 |                                 1.41  |
| loso     | C2/m      | formation |  30 |       0.08  |           0.067 |                                 1.206 |
| loso     | Fd-3m     | bandgap   |   9 |       0.088 |           0.092 |                                 0.956 |
| loso     | Fd-3m     | bandgap   |  10 |       0.08  |           0.091 |                                 0.879 |
| loso     | Fd-3m     | bandgap   |  20 |       0.071 |           0.075 |                                 0.946 |
| loso     | Fd-3m     | bandgap   |  30 |       0.06  |           0.065 |                                 0.925 |
| loso     | Fd-3m     | formation |   9 |       0.108 |           0.088 |                                 1.218 |
| loso     | Fd-3m     | formation |  10 |       0.105 |           0.093 |                                 1.136 |
| loso     | Fd-3m     | formation |  20 |       0.075 |           0.069 |                                 1.086 |
| loso     | Fd-3m     | formation |  30 |       0.065 |           0.058 |                                 1.119 |
| loso     | Ia-3      | bandgap   |   9 |       0.269 |           0.091 |                                 2.97  |
| loso     | Ia-3      | bandgap   |  10 |       0.253 |           0.082 |                                 3.095 |
| loso     | Ia-3      | bandgap   |  20 |       0.176 |           0.074 |                                 2.375 |
| loso     | Ia-3      | bandgap   |  30 |       0.147 |           0.061 |                                 2.399 |
| loso     | Ia-3      | formation |   9 |       0.212 |           0.104 |                                 2.048 |
| loso     | Ia-3      | formation |  10 |       0.227 |           0.09  |                                 2.52  |
| loso     | Ia-3      | formation |  20 |       0.216 |           0.071 |                                 3.057 |
| loso     | Ia-3      | formation |  30 |       0.185 |           0.071 |                                 2.616 |
| loso     | P63/mmc   | bandgap   |   9 |       0.16  |           0.099 |                                 1.619 |
| loso     | P63/mmc   | bandgap   |  10 |       0.153 |           0.087 |                                 1.748 |
| loso     | P63/mmc   | bandgap   |  20 |       0.103 |           0.068 |                                 1.526 |
| loso     | P63/mmc   | bandgap   |  30 |       0.093 |           0.062 |                                 1.495 |
| loso     | P63/mmc   | formation |   9 |       0.108 |           0.097 |                                 1.107 |
| loso     | P63/mmc   | formation |  10 |       0.099 |           0.085 |                                 1.165 |
| loso     | P63/mmc   | formation |  20 |       0.074 |           0.069 |                                 1.075 |
| loso     | P63/mmc   | formation |  30 |       0.062 |           0.06  |                                 1.043 |
| loso     | Pna21     | bandgap   |   9 |       0.139 |           0.089 |                                 1.561 |
| loso     | Pna21     | bandgap   |  10 |       0.131 |           0.09  |                                 1.453 |
| loso     | Pna21     | bandgap   |  20 |       0.109 |           0.068 |                                 1.602 |
| loso     | Pna21     | bandgap   |  30 |       0.09  |           0.059 |                                 1.535 |
| loso     | Pna21     | formation |   9 |       0.152 |           0.089 |                                 1.706 |
| loso     | Pna21     | formation |  10 |       0.143 |           0.088 |                                 1.616 |
| loso     | Pna21     | formation |  20 |       0.112 |           0.064 |                                 1.755 |
| loso     | Pna21     | formation |  30 |       0.083 |           0.057 |                                 1.457 |
| loso     | R-3c      | bandgap   |   9 |       0.187 |           0.093 |                                 2.015 |
| loso     | R-3c      | bandgap   |  10 |       0.169 |           0.083 |                                 2.033 |
| loso     | R-3c      | bandgap   |  20 |       0.147 |           0.074 |                                 1.986 |
| loso     | R-3c      | bandgap   |  30 |       0.118 |           0.07  |                                 1.691 |
| loso     | R-3c      | formation |   9 |       0.161 |           0.093 |                                 1.725 |
| loso     | R-3c      | formation |  10 |       0.15  |           0.089 |                                 1.673 |
| loso     | R-3c      | formation |  20 |       0.119 |           0.071 |                                 1.69  |
| loso     | R-3c      | formation |  30 |       0.106 |           0.068 |                                 1.566 |

## Fd-3m sub-family transfer (exploratory)

| target    | scheme   |   m | method      |   mean_cov |   q05 |   p_below |   width_median |   b_hat |
|:----------|:---------|----:|:------------|-----------:|------:|----------:|---------------:|--------:|
| bandgap   | sub_a2b  |   9 | split       |      0.886 | 0.678 |     0.124 |          2.67  |  -0.642 |
| bandgap   | sub_a2b  |   9 | offset_full |      0.891 | 0.726 |     0.128 |          1.671 |  -0.642 |
| bandgap   | sub_a2b  |  10 | split       |      0.886 | 0.677 |     0.114 |          2.715 |  -0.626 |
| bandgap   | sub_a2b  |  10 | offset_full |      0.889 | 0.719 |     0.142 |          1.699 |  -0.626 |
| bandgap   | sub_a2b  |  20 | split       |      0.884 | 0.781 |     0.068 |          2.614 |  -0.62  |
| bandgap   | sub_a2b  |  20 | offset_full |      0.883 | 0.781 |     0.088 |          1.502 |  -0.62  |
| bandgap   | sub_a2b  |  30 | split       |      0.884 | 0.801 |     0.046 |          2.539 |  -0.625 |
| bandgap   | sub_a2b  |  30 | offset_full |      0.881 | 0.794 |     0.06  |          1.432 |  -0.625 |
| bandgap   | sub_a2b  |  50 | split       |      0.888 | 0.815 |     0.01  |          2.539 |  -0.629 |
| bandgap   | sub_a2b  |  50 | offset_full |      0.887 | 0.815 |     0.022 |          1.416 |  -0.629 |
| bandgap   | sub_a2b  | 100 | split       |      0.887 | 0.842 |     0     |          2.535 |  -0.627 |
| bandgap   | sub_a2b  | 100 | offset_full |      0.888 | 0.842 |     0     |          1.406 |  -0.627 |
| bandgap   | sub_b2a  |   9 | split       |      0.926 | 0.765 |     0.102 |          3.05  |  -0.673 |
| bandgap   | sub_b2a  |   9 | offset_full |      0.927 | 0.729 |     0.108 |          1.814 |  -0.673 |
| bandgap   | sub_b2a  |  10 | split       |      0.935 | 0.785 |     0.084 |          3.091 |  -0.674 |
| bandgap   | sub_b2a  |  10 | offset_full |      0.935 | 0.769 |     0.086 |          1.83  |  -0.674 |
| bandgap   | sub_b2a  |  20 | split       |      0.926 | 0.806 |     0.038 |          2.671 |  -0.67  |
| bandgap   | sub_b2a  |  20 | offset_full |      0.931 | 0.818 |     0.036 |          1.614 |  -0.67  |
| bandgap   | sub_b2a  |  30 | split       |      0.926 | 0.846 |     0.016 |          2.65  |  -0.673 |
| bandgap   | sub_b2a  |  30 | offset_full |      0.933 | 0.85  |     0.018 |          1.544 |  -0.673 |
| bandgap   | sub_b2a  |  50 | split       |      0.922 | 0.846 |     0.006 |          2.639 |  -0.675 |
| bandgap   | sub_b2a  |  50 | offset_full |      0.929 | 0.87  |     0.002 |          1.494 |  -0.675 |
| bandgap   | sub_b2a  | 100 | split       |      0.917 | 0.883 |     0     |          2.639 |  -0.671 |
| bandgap   | sub_b2a  | 100 | offset_full |      0.927 | 0.891 |     0     |          1.438 |  -0.671 |
| formation | sub_a2b  |   9 | split       |      0.766 | 0.555 |     0.538 |          0.393 |   0.099 |
| formation | sub_a2b  |   9 | offset_full |      0.847 | 0.623 |     0.296 |          0.261 |   0.099 |
| formation | sub_a2b  |  10 | split       |      0.774 | 0.603 |     0.548 |          0.393 |   0.098 |
| formation | sub_a2b  |  10 | offset_full |      0.847 | 0.651 |     0.29  |          0.255 |   0.098 |
| formation | sub_a2b  |  20 | split       |      0.75  | 0.616 |     0.66  |          0.386 |   0.098 |
| formation | sub_a2b  |  20 | offset_full |      0.827 | 0.658 |     0.302 |          0.238 |   0.098 |
| formation | sub_a2b  |  30 | split       |      0.748 | 0.623 |     0.71  |          0.386 |   0.1   |
| formation | sub_a2b  |  30 | offset_full |      0.824 | 0.699 |     0.334 |          0.23  |   0.1   |
| formation | sub_a2b  |  50 | split       |      0.733 | 0.637 |     0.832 |          0.381 |   0.099 |
| formation | sub_a2b  |  50 | offset_full |      0.817 | 0.719 |     0.322 |          0.227 |   0.099 |
| formation | sub_a2b  | 100 | split       |      0.726 | 0.671 |     0.966 |          0.376 |   0.099 |
| formation | sub_a2b  | 100 | offset_full |      0.811 | 0.74  |     0.314 |          0.223 |   0.099 |
| formation | sub_b2a  |   9 | split       |      0.972 | 0.874 |     0.016 |          0.514 |   0.147 |
| formation | sub_b2a  |   9 | offset_full |      0.821 | 0.538 |     0.398 |          0.252 |   0.147 |
| formation | sub_b2a  |  10 | split       |      0.979 | 0.915 |     0.006 |          0.518 |   0.148 |
| formation | sub_b2a  |  10 | offset_full |      0.839 | 0.575 |     0.348 |          0.258 |   0.148 |
| formation | sub_b2a  |  20 | split       |      0.982 | 0.939 |     0     |          0.49  |   0.148 |
| formation | sub_b2a  |  20 | offset_full |      0.812 | 0.632 |     0.462 |          0.236 |   0.148 |
| formation | sub_b2a  |  30 | split       |      0.982 | 0.947 |     0     |          0.477 |   0.148 |
| formation | sub_b2a  |  30 | offset_full |      0.801 | 0.64  |     0.49  |          0.23  |   0.148 |
| formation | sub_b2a  |  50 | split       |      0.983 | 0.955 |     0     |          0.474 |   0.148 |
| formation | sub_b2a  |  50 | offset_full |      0.795 | 0.664 |     0.538 |          0.229 |   0.148 |
| formation | sub_b2a  | 100 | split       |      0.981 | 0.972 |     0     |          0.474 |   0.147 |
| formation | sub_b2a  | 100 | offset_full |      0.801 | 0.737 |     0.5   |          0.227 |   0.147 |

## Training arm (descriptive)

| design   | heldout   | target    |   m | arm       | method      |   n_draws |   mean_cov |   sd_cov |   mean_cov_group |   width_median |   mae |
|:---------|:----------|:----------|----:|:----------|:------------|----------:|-----------:|---------:|-----------------:|---------------:|------:|
| l2o      | P63/mmc   | bandgap   |   0 | baseline  | offset_full |        10 |      0.937 |    0.057 |            0.93  |          1.446 | 0.434 |
| l2o      | P63/mmc   | bandgap   |   0 | baseline  | split       |        10 |      0.925 |    0.073 |            0.907 |          2.309 | 0.434 |
| l2o      | P63/mmc   | bandgap   |   0 | fixed     | offset_full |        10 |      0.94  |    0.062 |            0.94  |          1.467 | 0.445 |
| l2o      | P63/mmc   | bandgap   |   0 | fixed     | split       |        10 |      0.94  |    0.057 |            0.931 |          2.329 | 0.445 |
| l2o      | P63/mmc   | bandgap   |  10 | calibrate | offset_full |        10 |      0.946 |    0.04  |            0.947 |          1.336 | 0.445 |
| l2o      | P63/mmc   | bandgap   |  10 | calibrate | split       |        10 |      0.93  |    0.066 |            0.916 |          1.9   | 0.445 |
| l2o      | P63/mmc   | bandgap   |  10 | train     | offset_full |        10 |      0.927 |    0.046 |            0.917 |          1.106 | 0.268 |
| l2o      | P63/mmc   | bandgap   |  10 | train     | split       |        10 |      0.908 |    0.078 |            0.892 |          1.3   | 0.268 |
| l2o      | P63/mmc   | bandgap   |  30 | calibrate | offset_full |        10 |      0.925 |    0.056 |            0.914 |          1.166 | 0.445 |
| l2o      | P63/mmc   | bandgap   |  30 | calibrate | split       |        10 |      0.927 |    0.056 |            0.904 |          1.875 | 0.445 |
| l2o      | P63/mmc   | bandgap   |  30 | train     | offset_full |        10 |      0.905 |    0.087 |            0.886 |          1.07  | 0.216 |
| l2o      | P63/mmc   | bandgap   |  30 | train     | split       |        10 |      0.906 |    0.083 |            0.888 |          1.148 | 0.216 |
| l2o      | P63/mmc   | bandgap   | 100 | calibrate | offset_full |        10 |      0.902 |    0.065 |            0.883 |          1.084 | 0.445 |
| l2o      | P63/mmc   | bandgap   | 100 | calibrate | split       |        10 |      0.91  |    0.033 |            0.885 |          1.625 | 0.445 |
| l2o      | P63/mmc   | bandgap   | 100 | train     | offset_full |        10 |      0.891 |    0.089 |            0.872 |          0.854 | 0.175 |
| l2o      | P63/mmc   | bandgap   | 100 | train     | split       |        10 |      0.892 |    0.095 |            0.869 |          0.959 | 0.175 |
| l2o      | P63/mmc   | formation |   0 | baseline  | offset_full |        10 |      0.846 |    0.114 |            0.843 |          0.247 | 0.105 |
| l2o      | P63/mmc   | formation |   0 | baseline  | split       |        10 |      0.906 |    0.053 |            0.911 |          0.349 | 0.105 |
| l2o      | P63/mmc   | formation |   0 | fixed     | offset_full |        10 |      0.84  |    0.122 |            0.835 |          0.244 | 0.114 |
| l2o      | P63/mmc   | formation |   0 | fixed     | split       |        10 |      0.905 |    0.053 |            0.91  |          0.38  | 0.114 |
| l2o      | P63/mmc   | formation |  10 | calibrate | offset_full |        10 |      0.856 |    0.127 |            0.854 |          0.268 | 0.114 |
| l2o      | P63/mmc   | formation |  10 | calibrate | split       |        10 |      0.897 |    0.051 |            0.906 |          0.369 | 0.114 |
| l2o      | P63/mmc   | formation |  10 | train     | offset_full |        10 |      0.842 |    0.12  |            0.85  |          0.207 | 0.067 |
| l2o      | P63/mmc   | formation |  10 | train     | split       |        10 |      0.855 |    0.079 |            0.865 |          0.227 | 0.067 |
| l2o      | P63/mmc   | formation |  30 | calibrate | offset_full |        10 |      0.878 |    0.074 |            0.882 |          0.272 | 0.114 |
| l2o      | P63/mmc   | formation |  30 | calibrate | split       |        10 |      0.92  |    0.03  |            0.921 |          0.375 | 0.114 |
| l2o      | P63/mmc   | formation |  30 | train     | offset_full |        10 |      0.878 |    0.082 |            0.875 |          0.173 | 0.046 |
| l2o      | P63/mmc   | formation |  30 | train     | split       |        10 |      0.875 |    0.075 |            0.873 |          0.174 | 0.046 |
| l2o      | P63/mmc   | formation | 100 | calibrate | offset_full |        10 |      0.868 |    0.072 |            0.87  |          0.259 | 0.114 |
| l2o      | P63/mmc   | formation | 100 | calibrate | split       |        10 |      0.907 |    0.037 |            0.912 |          0.367 | 0.114 |
| l2o      | P63/mmc   | formation | 100 | train     | offset_full |        10 |      0.916 |    0.068 |            0.906 |          0.145 | 0.031 |
| l2o      | P63/mmc   | formation | 100 | train     | split       |        10 |      0.913 |    0.048 |            0.898 |          0.141 | 0.031 |
| loso     | Fd-3m     | bandgap   |   0 | baseline  | offset_full |        10 |      0.899 |    0.085 |            0.909 |          1.648 | 0.642 |
| loso     | Fd-3m     | bandgap   |   0 | baseline  | split       |        10 |      0.919 |    0.051 |            0.929 |          2.764 | 0.642 |
| loso     | Fd-3m     | bandgap   |   0 | fixed     | offset_full |        10 |      0.898 |    0.096 |            0.91  |          1.557 | 0.656 |
| loso     | Fd-3m     | bandgap   |   0 | fixed     | split       |        10 |      0.911 |    0.056 |            0.926 |          2.634 | 0.656 |
| loso     | Fd-3m     | bandgap   |  10 | calibrate | offset_full |        10 |      0.914 |    0.05  |            0.931 |          1.511 | 0.656 |
| loso     | Fd-3m     | bandgap   |  10 | calibrate | split       |        10 |      0.915 |    0.053 |            0.931 |          2.634 | 0.656 |
| loso     | Fd-3m     | bandgap   |  10 | train     | offset_full |        10 |      0.886 |    0.08  |            0.899 |          1.644 | 0.497 |
| loso     | Fd-3m     | bandgap   |  10 | train     | split       |        10 |      0.891 |    0.073 |            0.908 |          2.343 | 0.497 |
| loso     | Fd-3m     | bandgap   |  30 | calibrate | offset_full |        10 |      0.885 |    0.047 |            0.901 |          1.414 | 0.656 |
| loso     | Fd-3m     | bandgap   |  30 | calibrate | split       |        10 |      0.893 |    0.047 |            0.911 |          2.535 | 0.656 |
| loso     | Fd-3m     | bandgap   |  30 | train     | offset_full |        10 |      0.885 |    0.086 |            0.91  |          1.72  | 0.416 |
| loso     | Fd-3m     | bandgap   |  30 | train     | split       |        10 |      0.882 |    0.083 |            0.905 |          1.894 | 0.416 |
| loso     | Fd-3m     | bandgap   | 100 | calibrate | offset_full |        10 |      0.872 |    0.045 |            0.89  |          1.358 | 0.656 |
| loso     | Fd-3m     | bandgap   | 100 | calibrate | split       |        10 |      0.873 |    0.055 |            0.892 |          2.5   | 0.656 |
| loso     | Fd-3m     | bandgap   | 100 | train     | offset_full |        10 |      0.888 |    0.094 |            0.918 |          1.626 | 0.374 |
| loso     | Fd-3m     | bandgap   | 100 | train     | split       |        10 |      0.878 |    0.096 |            0.912 |          1.497 | 0.374 |
| loso     | Fd-3m     | formation |   0 | baseline  | offset_full |        10 |      0.902 |    0.105 |            0.904 |          0.279 | 0.122 |
| loso     | Fd-3m     | formation |   0 | baseline  | split       |        10 |      0.889 |    0.073 |            0.921 |          0.424 | 0.122 |
| loso     | Fd-3m     | formation |   0 | fixed     | offset_full |        10 |      0.91  |    0.086 |            0.909 |          0.255 | 0.119 |
| loso     | Fd-3m     | formation |   0 | fixed     | split       |        10 |      0.898 |    0.075 |            0.928 |          0.429 | 0.119 |
| loso     | Fd-3m     | formation |  10 | calibrate | offset_full |        10 |      0.922 |    0.072 |            0.922 |          0.265 | 0.119 |
| loso     | Fd-3m     | formation |  10 | calibrate | split       |        10 |      0.895 |    0.08  |            0.926 |          0.425 | 0.119 |
| loso     | Fd-3m     | formation |  10 | train     | offset_full |        10 |      0.936 |    0.054 |            0.924 |          0.28  | 0.077 |
| loso     | Fd-3m     | formation |  10 | train     | split       |        10 |      0.913 |    0.062 |            0.934 |          0.333 | 0.077 |
| loso     | Fd-3m     | formation |  30 | calibrate | offset_full |        10 |      0.919 |    0.061 |            0.914 |          0.26  | 0.119 |
| loso     | Fd-3m     | formation |  30 | calibrate | split       |        10 |      0.893 |    0.073 |            0.924 |          0.422 | 0.119 |
| loso     | Fd-3m     | formation |  30 | train     | offset_full |        10 |      0.913 |    0.095 |            0.91  |          0.25  | 0.06  |
| loso     | Fd-3m     | formation |  30 | train     | split       |        10 |      0.915 |    0.096 |            0.924 |          0.264 | 0.06  |
| loso     | Fd-3m     | formation | 100 | calibrate | offset_full |        10 |      0.893 |    0.042 |            0.89  |          0.235 | 0.119 |
| loso     | Fd-3m     | formation | 100 | calibrate | split       |        10 |      0.898 |    0.038 |            0.926 |          0.416 | 0.119 |
| loso     | Fd-3m     | formation | 100 | train     | offset_full |        10 |      0.88  |    0.105 |            0.862 |          0.207 | 0.05  |
| loso     | Fd-3m     | formation | 100 | train     | split       |        10 |      0.882 |    0.119 |            0.88  |          0.216 | 0.05  |
