# Phase 6 Task 6.4 report (JARVIS-DFT)

## Units

* formation primary: 10/10 units
* formation regimeA: 10/10 units
* formation randinner: 10/10 units
* formation polyex: 10/10 units
* formation h3: 10/10 units
* formation rep2: 10/10 units
* formation locso: 7/7 units
* bandgap primary: 10/10 units
* bandgap regimeA: 10/10 units
* bandgap randinner: 10/10 units
* bandgap polyex: 10/10 units
* bandgap h3: 10/10 units
* bandgap rep2: 10/10 units
* bandgap locso: 7/7 units
* tbmbj primary: 10/10 units

## formation

Primary (CV+, space-group inner folds, alpha 0.1; mean over repeats): 59 groups; median coverage 0.972; below 0.90: 14%; below 0.80: 5%. **G5** (>= 75% of groups >= 0.85): holds (92%).

Worst 10 groups:

|   group |   coverage |   n_repeats |    n |   mae |   bias |   width_sd |
|--------:|-----------:|------------:|-----:|------:|-------:|-----------:|
|     216 |      0.582 |           2 | 4533 | 0.338 |  0.293 |      0.829 |
|       6 |      0.704 |           2 |  790 | 0.425 |  0.385 |      0.445 |
|     187 |      0.771 |           2 | 3487 | 0.266 |  0.198 |      0.864 |
|      99 |      0.809 |           2 | 3478 | 0.235 |  0.139 |      0.756 |
|     221 |      0.841 |           2 | 4665 | 0.189 |  0.087 |      0.611 |
|     129 |      0.864 |           2 | 1263 | 0.181 | -0.132 |      0.777 |
|      25 |      0.892 |           2 |  499 | 0.168 |  0.093 |      0.696 |
|      47 |      0.896 |           2 |  527 | 0.161 |  0.03  |      0.752 |
|      71 |      0.921 |           2 | 1113 | 0.133 |  0.01  |      0.549 |
|     191 |      0.925 |           2 |  951 | 0.138 |  0.037 |      1.054 |

Per-group summary by design and method (alpha 0.10):

| design    | method     |   median_coverage |   share_below_090 |   share_below_080 |   median_width_sd |   groups |
|:----------|:-----------|------------------:|------------------:|------------------:|------------------:|---------:|
| locso     | cvplus     |             0.967 |             0.286 |             0     |             0.718 |        7 |
| locso     | enssd      |             0.569 |             1     |             1     |             0.175 |        7 |
| locso     | ncvplus    |             0.96  |             0.143 |             0     |             0.785 |        7 |
| polyex    | cvplus     |             0.964 |             0.153 |             0.068 |             0.736 |       59 |
| polyex    | cvplus_gw  |             0.953 |             0.186 |             0.068 |             0.661 |       59 |
| polyex    | enssd      |             0.543 |             1     |             1     |             0.168 |       59 |
| polyex    | ncvplus    |             0.96  |             0.119 |             0.051 |             0.735 |       59 |
| primary   | cqr_cvplus |             0.969 |             0.085 |             0.034 |             0.643 |       59 |
| primary   | cvplus     |             0.971 |             0.136 |             0.051 |             0.744 |       59 |
| primary   | cvplus_gw  |             0.962 |             0.153 |             0.068 |             0.654 |       59 |
| primary   | enssd      |             0.557 |             1     |             1     |             0.16  |       59 |
| primary   | ncvplus    |             0.967 |             0.102 |             0.034 |             0.732 |       59 |
| randinner | cvplus     |             0.908 |             0.458 |             0.136 |             0.44  |       59 |
| rep2      | cvplus     |             0.972 |             0.119 |             0.051 |             0.742 |       59 |
| rep2      | cvplus_gw  |             0.963 |             0.136 |             0.068 |             0.662 |       59 |
| rep2      | enssd      |             0.56  |             1     |             1     |             0.158 |       59 |
| rep2      | ncvplus    |             0.964 |             0.102 |             0.034 |             0.72  |       59 |

**G3-J**: 3 groups below 0.8 under CV+ (repeat 1); normalized or CQR lifts 0% to >= 0.85 -> negative.

Regime A: marginal CV+ random coverage 0.924; groups below 0.85: 7%; Mondrian >= 0.88: 100%. **H1-J**: not supported.

| method        |   coverage |   width |    mae |
|:--------------|-----------:|--------:|-------:|
| cqr_split     |     0.8999 |  0.4716 | 0.0919 |
| cvplus_random |     0.924  |  0.4154 | 0.0919 |
| mondrian      |     0.9018 |  0.3979 | 0.0919 |
| nsplit        |     0.8998 |  0.429  | 0.0919 |
| split         |     0.9    |  0.4122 | 0.0919 |

H3 (valid group-level design): group-weighted coverage {"hcp": 0.943, "single_subsample": 0.905}; share of test groups with own coverage < 0.80: {"hcp": 0.056, "single_subsample": 0.136}. **H3-J**: supported.

| method           |   alpha |   mean |    std |
|:-----------------|--------:|-------:|-------:|
| hcp              |    0.05 | 0.9964 | 0.0044 |
| hcp              |    0.1  | 0.9434 | 0.0283 |
| hcp              |    0.2  | 0.8398 | 0.0415 |
| single_subsample |    0.05 | 0.9527 | 0.0493 |
| single_subsample |    0.1  | 0.9049 | 0.0701 |
| single_subsample |    0.2  | 0.8133 | 0.0894 |

## bandgap

Primary (CV+, space-group inner folds, alpha 0.1; mean over repeats): 31 groups; median coverage 0.927; below 0.90: 23%; below 0.80: 6%. **G5** (>= 75% of groups >= 0.85): holds (84%).

Worst 10 groups:

|   group |   coverage |   n_repeats |    n |   mae |   bias |   width_sd |
|--------:|-----------:|------------:|-----:|------:|-------:|-----------:|
|     129 |      0.784 |           2 |  313 | 0.699 |  0.154 |      1.465 |
|      99 |      0.789 |           2 |  279 | 0.763 | -0.409 |      2.168 |
|     225 |      0.805 |           2 | 1339 | 0.695 |  0.276 |      1.201 |
|     221 |      0.827 |           2 |  346 | 0.672 | -0.171 |      1.197 |
|     123 |      0.846 |           2 |  276 | 0.602 | -0.327 |      1.521 |
|       6 |      0.868 |           2 |  285 | 0.541 | -0.414 |      1.6   |
|     187 |      0.898 |           2 |  215 | 0.508 | -0.209 |      1.505 |
|       8 |      0.904 |           2 |  382 | 0.493 | -0.259 |      1.552 |
|     139 |      0.915 |           2 |  293 | 0.477 | -0.084 |      1.398 |
|       5 |      0.915 |           2 |  317 | 0.433 | -0.034 |      1.261 |

Per-group summary by design and method (alpha 0.10):

| design    | method     |   median_coverage |   share_below_090 |   share_below_080 |   median_width_sd |   groups |
|:----------|:-----------|------------------:|------------------:|------------------:|------------------:|---------:|
| locso     | cvplus     |             0.922 |             0.286 |             0     |             1.381 |        7 |
| locso     | enssd      |             0.363 |             1     |             1     |             0.305 |        7 |
| locso     | ncvplus    |             0.919 |             0.143 |             0     |             1.516 |        7 |
| polyex    | cvplus     |             0.921 |             0.323 |             0.097 |             1.375 |       31 |
| polyex    | cvplus_gw  |             0.947 |             0.161 |             0     |             1.606 |       31 |
| polyex    | enssd      |             0.332 |             1     |             1     |             0.275 |       31 |
| polyex    | ncvplus    |             0.91  |             0.452 |             0     |             1.522 |       31 |
| primary   | cqr_cvplus |             0.924 |             0.258 |             0.032 |             1.502 |       31 |
| primary   | cvplus     |             0.925 |             0.194 |             0.065 |             1.371 |       31 |
| primary   | cvplus_gw  |             0.96  |             0.161 |             0     |             1.62  |       31 |
| primary   | enssd      |             0.344 |             1     |             1     |             0.257 |       31 |
| primary   | ncvplus    |             0.923 |             0.387 |             0     |             1.51  |       31 |
| randinner | cvplus     |             0.886 |             0.71  |             0.161 |             1.099 |       31 |
| rep2      | cvplus     |             0.928 |             0.226 |             0.097 |             1.354 |       31 |
| rep2      | cvplus_gw  |             0.956 |             0.161 |             0     |             1.617 |       31 |
| rep2      | enssd      |             0.36  |             1     |             1     |             0.265 |       31 |
| rep2      | ncvplus    |             0.921 |             0.355 |             0     |             1.509 |       31 |

**G3-J**: 2 groups below 0.8 under CV+ (repeat 1); normalized or CQR lifts 100% to >= 0.85 -> positive.

Regime A: marginal CV+ random coverage 0.917; groups below 0.85: 3%; Mondrian >= 0.88: 97%. **H1-J**: heterogeneity absent.

| method        |   coverage |   width |    mae |
|:--------------|-----------:|--------:|-------:|
| cqr_split     |     0.8993 |  2.3832 | 0.4392 |
| cvplus_random |     0.9169 |  1.9889 | 0.4392 |
| mondrian      |     0.9046 |  2.1139 | 0.4392 |
| nsplit        |     0.9018 |  2.1394 | 0.4392 |
| split         |     0.9017 |  2.0097 | 0.4392 |

H3 (valid group-level design): group-weighted coverage {"hcp": 0.981, "single_subsample": 0.915}; share of test groups with own coverage < 0.80: {"hcp": 0.0, "single_subsample": 0.108}. **H3-J**: supported.

| method           |   alpha |   mean |    std |
|:-----------------|--------:|-------:|-------:|
| hcp              |    0.05 | 1      | 0      |
| hcp              |    0.1  | 0.9807 | 0.0086 |
| hcp              |    0.2  | 0.8705 | 0.0333 |
| single_subsample |    0.05 | 1      | 0      |
| single_subsample |    0.1  | 0.9154 | 0.0801 |
| single_subsample |    0.2  | 0.8309 | 0.1093 |

## tbmbj

Primary (CV+, space-group inner folds, alpha 0.1; mean over repeats): 10 groups; median coverage 0.917; below 0.90: 40%; below 0.80: 0%. **G5** (>= 75% of groups >= 0.85): holds (90%).

Worst 10 groups:

|   group |   coverage |   n_repeats |   n |   mae |   bias |   width_sd |
|--------:|-----------:|------------:|----:|------:|-------:|-----------:|
|     225 |      0.839 |           1 | 639 | 0.908 |  0.032 |      0.991 |
|     221 |      0.871 |           1 | 249 | 0.893 | -0.346 |      1.108 |
|       2 |      0.878 |           1 | 229 | 0.817 | -0.078 |      1.377 |
|      62 |      0.895 |           1 | 258 | 0.786 |  0.227 |      1.282 |
|      14 |      0.912 |           1 | 284 | 0.81  |  0.187 |      1.144 |
|      12 |      0.922 |           1 | 308 | 0.61  | -0.068 |      1.604 |
|     166 |      0.93  |           1 | 330 | 0.761 |  0.372 |      1.824 |
|      15 |      0.942 |           1 | 241 | 0.661 |  0.021 |      1.502 |
|     164 |      0.945 |           1 | 238 | 0.676 |  0.267 |      1.314 |
|     216 |      0.964 |           1 | 307 | 0.545 |  0.166 |      1.727 |

Per-group summary by design and method (alpha 0.10):

| design   | method    |   median_coverage |   share_below_090 |   share_below_080 |   median_width_sd |   groups |
|:---------|:----------|------------------:|------------------:|------------------:|------------------:|---------:|
| primary  | cvplus    |             0.917 |               0.4 |                 0 |             1.345 |       10 |
| primary  | cvplus_gw |             1     |               0   |                 0 |            16.124 |       10 |
| primary  | enssd     |             0.243 |               1   |                 1 |             0.191 |       10 |
| primary  | ncvplus   |             0.894 |               0.6 |                 0 |             1.412 |       10 |
