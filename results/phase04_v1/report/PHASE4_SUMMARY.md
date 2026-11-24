# Phase 4 summary (primary outcome: miscoverage by cvplus_lattice, alpha = 0.1)

## bandgap: miscovered (pooled over the 6 LOSO holdouts; 95 % cluster-bootstrap CI)

| score      |   auroc |   auroc_lo |   auroc_hi |    ap |   prevalence |    n |   n_pos |   within_weighted |   within_lo |   within_hi |
|:-----------|--------:|-----------:|-----------:|------:|-------------:|-----:|--------:|------------------:|------------:|------------:|
| width      |   0.08  |      0.054 |      0.114 | 0.067 |        0.123 | 2343 |     289 |             0.534 |       0.432 |       0.617 |
| knn        |   0.434 |      0.376 |      0.491 | 0.101 |        0.123 | 2343 |     289 |             0.404 |       0.333 |       0.491 |
| kde        |   0.404 |      0.35  |      0.457 | 0.095 |        0.123 | 2343 |     289 |             0.401 |       0.336 |       0.475 |
| committee  |   0.15  |      0.112 |      0.198 | 0.071 |        0.123 | 2343 |     289 |             0.563 |       0.471 |       0.637 |
| subsample  |   0.181 |      0.133 |      0.241 | 0.072 |        0.123 | 2343 |     289 |             0.476 |       0.408 |       0.537 |
| classifier |   0.549 |      0.484 |      0.614 | 0.201 |        0.123 | 2343 |     289 |             0.546 |       0.485 |       0.612 |

Within-lattice AUROC per holdout (eligible if >= 10 structures in each class):

| heldout   |   n |   n_pos | eligible   |   width |     knn |     kde |   committee |   subsample |   classifier |
|:----------|----:|--------:|:-----------|--------:|--------:|--------:|------------:|------------:|-------------:|
| Fd-3m     | 393 |     260 | True       |   0.519 |   0.466 |   0.506 |       0.538 |       0.487 |        0.514 |
| C2/m      | 343 |       1 | False      | nan     | nan     | nan     |     nan     |     nan     |      nan     |
| Pna21     | 427 |       0 | False      | nan     | nan     | nan     |     nan     |     nan     |      nan     |
| R-3c      | 350 |       9 | False      | nan     | nan     | nan     |     nan     |     nan     |      nan     |
| P63/mmc   | 351 |      15 | True       |   0.55  |   0.336 |   0.283 |       0.591 |       0.465 |        0.581 |
| Ia-3      | 479 |       4 | False      | nan     | nan     | nan     |     nan     |     nan     |      nan     |

## bandgap: large_error (pooled over the 6 LOSO holdouts; 95 % cluster-bootstrap CI)

| score      |   auroc |   auroc_lo |   auroc_hi |    ap |   prevalence |    n |   n_pos |   within_weighted |   within_lo |   within_hi |
|:-----------|--------:|-----------:|-----------:|------:|-------------:|-----:|--------:|------------------:|------------:|------------:|
| width      |   0.249 |      0.2   |      0.305 | 0.243 |        0.358 | 2343 |     838 |             0.624 |       0.559 |       0.678 |
| knn        |   0.599 |      0.553 |      0.646 | 0.423 |        0.358 | 2343 |     838 |             0.497 |       0.454 |       0.545 |
| kde        |   0.546 |      0.492 |      0.599 | 0.392 |        0.358 | 2343 |     838 |             0.512 |       0.472 |       0.557 |
| committee  |   0.402 |      0.344 |      0.463 | 0.303 |        0.358 | 2343 |     838 |             0.639 |       0.572 |       0.695 |
| subsample  |   0.342 |      0.291 |      0.394 | 0.273 |        0.358 | 2343 |     838 |             0.636 |       0.592 |       0.678 |
| classifier |   0.505 |      0.461 |      0.549 | 0.381 |        0.358 | 2343 |     838 |             0.476 |       0.432 |       0.527 |

Within-lattice AUROC per holdout (eligible if >= 10 structures in each class):

| heldout   |   n |   n_pos | eligible   |   width |     knn |     kde |   committee |   subsample |   classifier |
|:----------|----:|--------:|:-----------|--------:|--------:|--------:|------------:|------------:|-------------:|
| Fd-3m     | 393 |     338 | True       |   0.521 |   0.482 |   0.535 |       0.556 |       0.482 |        0.493 |
| C2/m      | 343 |     193 | True       |   0.675 |   0.615 |   0.547 |       0.681 |       0.608 |        0.508 |
| Pna21     | 427 |       0 | False      | nan     | nan     | nan     |     nan     |     nan     |      nan     |
| R-3c      | 350 |     151 | True       |   0.768 |   0.512 |   0.557 |       0.8   |       0.77  |        0.552 |
| P63/mmc   | 351 |     110 | True       |   0.402 |   0.592 |   0.565 |       0.413 |       0.536 |        0.575 |
| Ia-3      | 479 |      46 | True       |   0.729 |   0.344 |   0.397 |       0.725 |       0.758 |        0.31  |

## bandgap: coverage of retained structures by abstention fraction

| score          |   0.0 |   0.1 |   0.25 |   0.5 |
|:---------------|------:|------:|-------:|------:|
| classifier     | 0.877 | 0.893 |  0.898 | 0.886 |
| committee      | 0.877 | 0.866 |  0.842 | 0.777 |
| kde            | 0.877 | 0.864 |  0.851 | 0.831 |
| knn            | 0.877 | 0.868 |  0.859 | 0.858 |
| oracle_abs_err | 0.877 | 0.954 |  0.995 | 1     |
| subsample      | 0.877 | 0.863 |  0.842 | 0.787 |
| width          | 0.877 | 0.864 |  0.84  | 0.765 |

## formation: miscovered (pooled over the 6 LOSO holdouts; 95 % cluster-bootstrap CI)

| score      |   auroc |   auroc_lo |   auroc_hi |    ap |   prevalence |    n |   n_pos |   within_weighted |   within_lo |   within_hi |
|:-----------|--------:|-----------:|-----------:|------:|-------------:|-----:|--------:|------------------:|------------:|------------:|
| width      |   0.082 |      0.064 |      0.103 | 0.049 |        0.091 | 2343 |     214 |             0.388 |       0.332 |       0.445 |
| knn        |   0.429 |      0.372 |      0.492 | 0.075 |        0.091 | 2343 |     214 |             0.448 |       0.39  |       0.506 |
| kde        |   0.403 |      0.347 |      0.463 | 0.07  |        0.091 | 2343 |     214 |             0.465 |       0.405 |       0.526 |
| committee  |   0.386 |      0.332 |      0.44  | 0.068 |        0.091 | 2343 |     214 |             0.437 |       0.385 |       0.491 |
| subsample  |   0.445 |      0.401 |      0.496 | 0.079 |        0.091 | 2343 |     214 |             0.581 |       0.525 |       0.635 |
| classifier |   0.62  |      0.541 |      0.694 | 0.196 |        0.091 | 2343 |     214 |             0.594 |       0.516 |       0.672 |

Within-lattice AUROC per holdout (eligible if >= 10 structures in each class):

| heldout   |   n |   n_pos | eligible   |   width |     knn |     kde |   committee |   subsample |   classifier |
|:----------|----:|--------:|:-----------|--------:|--------:|--------:|------------:|------------:|-------------:|
| Fd-3m     | 393 |     180 | True       |   0.472 |   0.393 |   0.431 |       0.504 |       0.495 |        0.605 |
| C2/m      | 343 |       1 | False      | nan     | nan     | nan     |     nan     |     nan     |      nan     |
| Pna21     | 427 |       0 | False      | nan     | nan     | nan     |     nan     |     nan     |      nan     |
| R-3c      | 350 |       0 | False      | nan     | nan     | nan     |     nan     |     nan     |      nan     |
| P63/mmc   | 351 |      30 | True       |   0.293 |   0.51  |   0.504 |       0.363 |       0.677 |        0.582 |
| Ia-3      | 479 |       3 | False      | nan     | nan     | nan     |     nan     |     nan     |      nan     |

## formation: large_error (pooled over the 6 LOSO holdouts; 95 % cluster-bootstrap CI)

| score      |   auroc |   auroc_lo |   auroc_hi |    ap |   prevalence |    n |   n_pos |   within_weighted |   within_lo |   within_hi |
|:-----------|--------:|-----------:|-----------:|------:|-------------:|-----:|--------:|------------------:|------------:|------------:|
| width      |   0.289 |      0.244 |      0.337 | 0.427 |        0.547 | 2343 |    1282 |             0.4   |       0.356 |       0.446 |
| knn        |   0.639 |      0.588 |      0.687 | 0.619 |        0.547 | 2343 |    1282 |             0.566 |       0.529 |       0.6   |
| kde        |   0.623 |      0.567 |      0.673 | 0.625 |        0.547 | 2343 |    1282 |             0.606 |       0.572 |       0.638 |
| committee  |   0.457 |      0.402 |      0.51  | 0.519 |        0.547 | 2343 |    1282 |             0.396 |       0.353 |       0.442 |
| subsample  |   0.505 |      0.459 |      0.549 | 0.565 |        0.547 | 2343 |    1282 |             0.511 |       0.466 |       0.556 |
| classifier |   0.472 |      0.434 |      0.513 | 0.532 |        0.547 | 2343 |    1282 |             0.467 |       0.426 |       0.514 |

Within-lattice AUROC per holdout (eligible if >= 10 structures in each class):

| heldout   |   n |   n_pos | eligible   |   width |   knn |   kde |   committee |   subsample |   classifier |
|:----------|----:|--------:|:-----------|--------:|------:|------:|------------:|------------:|-------------:|
| Fd-3m     | 393 |     344 | True       |   0.486 | 0.344 | 0.365 |       0.483 |       0.437 |        0.561 |
| C2/m      | 343 |     265 | True       |   0.548 | 0.692 | 0.656 |       0.546 |       0.489 |        0.557 |
| Pna21     | 427 |      90 | True       |   0.418 | 0.62  | 0.626 |       0.407 |       0.548 |        0.376 |
| R-3c      | 350 |     196 | True       |   0.343 | 0.472 | 0.589 |       0.378 |       0.579 |        0.44  |
| P63/mmc   | 351 |     266 | True       |   0.286 | 0.402 | 0.423 |       0.337 |       0.499 |        0.502 |
| Ia-3      | 479 |     121 | True       |   0.334 | 0.801 | 0.898 |       0.265 |       0.513 |        0.401 |

## formation: coverage of retained structures by abstention fraction

| score          |   0.0 |   0.1 |   0.25 |   0.5 |
|:---------------|------:|------:|-------:|------:|
| classifier     | 0.909 | 0.925 |  0.937 | 0.933 |
| committee      | 0.909 | 0.899 |  0.883 | 0.881 |
| kde            | 0.909 | 0.901 |  0.889 | 0.879 |
| knn            | 0.909 | 0.903 |  0.894 | 0.89  |
| oracle_abs_err | 0.909 | 0.975 |  1     | 1     |
| subsample      | 0.909 | 0.902 |  0.902 | 0.901 |
| width          | 0.909 | 0.899 |  0.878 | 0.822 |

## Flag (kNN-distance conformal p-values, p <= 0.05)

| design   | heldout   |   n_cal |   n_idtest |   false_flag_rate |   n_test |   flag_rate_heldout |   ffr_ci_lo |   ffr_ci_hi |
|:---------|:----------|--------:|-----------:|------------------:|---------:|--------------------:|------------:|------------:|
| loso     | Fd-3m     |     421 |        346 |             0.084 |      393 |               0.148 |       0.057 |       0.118 |
| loso     | C2/m      |     382 |        435 |             0.023 |      343 |               0.087 |       0.011 |       0.042 |
| loso     | Pna21     |     336 |        395 |             0.046 |      427 |               0.026 |       0.027 |       0.071 |
| loso     | R-3c      |     424 |        429 |             0.065 |      350 |               0.586 |       0.044 |       0.093 |
| loso     | P63/mmc   |     368 |        453 |             0.099 |      351 |               0.393 |       0.073 |       0.131 |
| loso     | Ia-3      |     397 |        328 |             0.037 |      479 |               0     |       0.019 |       0.063 |
| l2o      | C2/m      |     349 |        287 |             0.08  |      343 |               0.35  |       0.051 |       0.118 |
| l2o      | Pna21     |     232 |        305 |             0.033 |      427 |               0.014 |       0.016 |       0.059 |
| l2o      | R-3c      |     359 |        285 |             0.053 |      350 |               0.74  |       0.03  |       0.085 |
| l2o      | P63/mmc   |     287 |        358 |             0.014 |      351 |               0.103 |       0.005 |       0.032 |
| l2o      | Ia-3      |     361 |        290 |             0.055 |      479 |               0.002 |       0.032 |       0.088 |

Miscoverage among flagged / unflagged held-out structures:

| design   | target    | heldout   |   miscov_if_not_flagged |   miscov_if_flagged |
|:---------|:----------|:----------|------------------------:|--------------------:|
| l2o      | bandgap   | C2/m      |                   0.009 |               0     |
| l2o      | bandgap   | Ia-3      |                   0.013 |               0     |
| l2o      | bandgap   | P63/mmc   |                   0.717 |               0.25  |
| l2o      | bandgap   | Pna21     |                   0     |               0     |
| l2o      | bandgap   | R-3c      |                   0.088 |               0.031 |
| l2o      | formation | C2/m      |                   0.004 |               0     |
| l2o      | formation | Ia-3      |                   0.004 |               0     |
| l2o      | formation | P63/mmc   |                   0.724 |               0.5   |
| l2o      | formation | Pna21     |                   0     |               0     |
| l2o      | formation | R-3c      |                   0     |               0     |
| loso     | bandgap   | C2/m      |                   0.003 |               0     |
| loso     | bandgap   | Fd-3m     |                   0.687 |               0.517 |
| loso     | bandgap   | Ia-3      |                   0.008 |             nan     |
| loso     | bandgap   | P63/mmc   |                   0.061 |               0.014 |
| loso     | bandgap   | Pna21     |                   0     |               0     |
| loso     | bandgap   | R-3c      |                   0.062 |               0     |
| loso     | formation | C2/m      |                   0     |               0.033 |
| loso     | formation | Fd-3m     |                   0.499 |               0.224 |
| loso     | formation | Ia-3      |                   0.006 |             nan     |
| loso     | formation | P63/mmc   |                   0.089 |               0.08  |
| loso     | formation | Pna21     |                   0     |               0     |
| loso     | formation | R-3c      |                   0     |               0     |

## Lattice level (descriptive; 11 configurations reuse 6 test lattices)

| design   | target    | heldout   |   n |   miscoverage |   abs_bias |   width_median |   knn_mean |   kde_mean |   classifier_mean |   committee_mean |   flagged |
|:---------|:----------|:----------|----:|--------------:|-----------:|---------------:|-----------:|-----------:|------------------:|-----------------:|----------:|
| loso     | bandgap   | Fd-3m     | 393 |         0.662 |      0.643 |          0.864 |     99.773 |      0.756 |             0.99  |            0.079 |     0.148 |
| loso     | formation | Fd-3m     | 393 |         0.458 |      0.117 |          0.259 |     99.773 |      0.756 |             0.99  |            0.022 |     0.148 |
| loso     | bandgap   | C2/m      | 343 |         0.003 |      0.367 |          1.901 |    103.742 |      0.762 |             0.991 |            0.262 |     0.087 |
| loso     | formation | C2/m      | 343 |         0.003 |      0.093 |          0.389 |    103.742 |      0.762 |             0.991 |            0.042 |     0.087 |
| loso     | bandgap   | Pna21     | 427 |         0     |      0     |          2.277 |     98.449 |      0.674 |             0.987 |            0.268 |     0.026 |
| loso     | formation | Pna21     | 427 |         0     |      0.031 |          0.417 |     98.449 |      0.674 |             0.987 |            0.034 |     0.026 |
| loso     | bandgap   | R-3c      | 350 |         0.026 |      0.384 |          1.524 |    132.477 |      0.992 |             0.992 |            0.151 |     0.586 |
| loso     | formation | R-3c      | 350 |         0     |      0.053 |          0.357 |    132.477 |      0.992 |             0.992 |            0.024 |     0.586 |
| loso     | bandgap   | P63/mmc   | 351 |         0.043 |      0.127 |          1.872 |    114.841 |      0.929 |             0.991 |            0.239 |     0.393 |
| loso     | formation | P63/mmc   | 351 |         0.085 |      0.064 |          0.355 |    114.841 |      0.929 |             0.991 |            0.026 |     0.393 |
| loso     | bandgap   | Ia-3      | 479 |         0.008 |      0.251 |          1.836 |     89.123 |      0.784 |             0.994 |            0.166 |     0     |
| loso     | formation | Ia-3      | 479 |         0.006 |      0.036 |          0.355 |     89.123 |      0.784 |             0.994 |            0.02  |     0     |
| l2o      | bandgap   | C2/m      | 343 |         0.006 |      0.026 |          1.144 |    113.851 |      0.914 |             0.992 |            0.147 |     0.35  |
| l2o      | formation | C2/m      | 343 |         0.003 |      0.03  |          0.291 |    113.851 |      0.914 |             0.992 |            0.01  |     0.35  |
| l2o      | bandgap   | Pna21     | 427 |         0     |      0.051 |          1.237 |     97.014 |      0.651 |             0.991 |            0.14  |     0.014 |
| l2o      | formation | Pna21     | 427 |         0     |      0.009 |          0.301 |     97.014 |      0.651 |             0.991 |            0.022 |     0.014 |
| l2o      | bandgap   | R-3c      | 350 |         0.046 |      0.18  |          1.097 |    138.224 |      0.993 |             0.992 |            0.107 |     0.74  |
| l2o      | formation | R-3c      | 350 |         0     |      0.042 |          0.288 |    138.224 |      0.993 |             0.992 |            0.029 |     0.74  |
| l2o      | bandgap   | P63/mmc   | 351 |         0.67  |      0.427 |          0.528 |    112.567 |      0.927 |             0.991 |            0.063 |     0.103 |
| l2o      | formation | P63/mmc   | 351 |         0.701 |      0.099 |          0.165 |    112.567 |      0.927 |             0.991 |            0.021 |     0.103 |
| l2o      | bandgap   | Ia-3      | 479 |         0.013 |      0.038 |          1.165 |     87.716 |      0.814 |             0.993 |            0.168 |     0.002 |
| l2o      | formation | Ia-3      | 479 |         0.004 |      0.021 |          0.364 |     87.716 |      0.814 |             0.993 |            0.023 |     0.002 |

| target    | configs   | score           | outcome     |   n |   spearman |   p_perm |
|:----------|:----------|:----------------|:------------|----:|-----------:|---------:|
| bandgap   | all 11    | width_median    | miscoverage |  11 |     -0.688 |    0.023 |
| bandgap   | all 11    | width_median    | abs_bias    |  11 |     -0.4   |    0.224 |
| bandgap   | all 11    | knn_mean        | miscoverage |  11 |      0.442 |    0.174 |
| bandgap   | all 11    | knn_mean        | abs_bias    |  11 |      0.245 |    0.469 |
| bandgap   | all 11    | kde_mean        | miscoverage |  11 |      0.624 |    0.045 |
| bandgap   | all 11    | kde_mean        | abs_bias    |  11 |      0.209 |    0.537 |
| bandgap   | all 11    | classifier_mean | miscoverage |  11 |      0.046 |    0.898 |
| bandgap   | all 11    | classifier_mean | abs_bias    |  11 |     -0.018 |    0.967 |
| bandgap   | all 11    | committee_mean  | miscoverage |  11 |     -0.624 |    0.045 |
| bandgap   | all 11    | committee_mean  | abs_bias    |  11 |     -0.518 |    0.107 |
| bandgap   | LOSO 6    | width_median    | miscoverage |   6 |     -0.829 |    0.058 |
| bandgap   | LOSO 6    | width_median    | abs_bias    |   6 |     -0.829 |    0.058 |
| bandgap   | LOSO 6    | knn_mean        | miscoverage |   6 |      0.371 |    0.497 |
| bandgap   | LOSO 6    | knn_mean        | abs_bias    |   6 |      0.314 |    0.564 |
| bandgap   | LOSO 6    | kde_mean        | miscoverage |   6 |      0.371 |    0.497 |
| bandgap   | LOSO 6    | kde_mean        | abs_bias    |   6 |      0.2   |    0.714 |
| bandgap   | LOSO 6    | classifier_mean | miscoverage |   6 |      0.029 |    1     |
| bandgap   | LOSO 6    | classifier_mean | abs_bias    |   6 |      0.257 |    0.658 |
| bandgap   | LOSO 6    | committee_mean  | miscoverage |   6 |     -0.829 |    0.058 |
| bandgap   | LOSO 6    | committee_mean  | abs_bias    |   6 |     -0.829 |    0.058 |
| formation | all 11    | width_median    | miscoverage |  11 |     -0.443 |    0.173 |
| formation | all 11    | width_median    | abs_bias    |  11 |     -0.345 |    0.299 |
| formation | all 11    | knn_mean        | miscoverage |  11 |     -0.168 |    0.618 |
| formation | all 11    | knn_mean        | abs_bias    |  11 |      0.418 |    0.204 |
| formation | all 11    | kde_mean        | miscoverage |  11 |      0.051 |    0.883 |
| formation | all 11    | kde_mean        | abs_bias    |  11 |      0.255 |    0.453 |
| formation | all 11    | classifier_mean | miscoverage |  11 |     -0.019 |    0.961 |
| formation | all 11    | classifier_mean | abs_bias    |  11 |     -0.264 |    0.435 |
| formation | all 11    | committee_mean  | miscoverage |  11 |     -0.462 |    0.154 |
| formation | all 11    | committee_mean  | abs_bias    |  11 |      0.127 |    0.715 |
| formation | LOSO 6    | width_median    | miscoverage |   6 |     -0.899 |    0.028 |
| formation | LOSO 6    | width_median    | abs_bias    |   6 |     -0.6   |    0.242 |
| formation | LOSO 6    | knn_mean        | miscoverage |   6 |     -0.145 |    0.778 |
| formation | LOSO 6    | knn_mean        | abs_bias    |   6 |      0.371 |    0.497 |
| formation | LOSO 6    | kde_mean        | miscoverage |   6 |     -0.058 |    0.933 |
| formation | LOSO 6    | kde_mean        | abs_bias    |   6 |      0.029 |    1     |
| formation | LOSO 6    | classifier_mean | miscoverage |   6 |     -0.087 |    0.861 |
| formation | LOSO 6    | classifier_mean | abs_bias    |   6 |     -0.086 |    0.919 |
| formation | LOSO 6    | committee_mean  | miscoverage |   6 |     -0.435 |    0.389 |
| formation | LOSO 6    | committee_mean  | abs_bias    |   6 |      0.029 |    1     |

## Phase 2 supplement: bias and one-sided misses

| design   | target    | heldout   | method         |   bias |   mae |   mae_without_offset |   coverage |   miss_below |   miss_above |
|:---------|:----------|:----------|:---------------|-------:|------:|---------------------:|-----------:|-------------:|-------------:|
| loso     | bandgap   | Fd-3m     | cvplus_lattice |  0.643 | 0.65  |                0.346 |      0.338 |        0.662 |        0     |
| loso     | formation | Fd-3m     | cvplus_lattice | -0.117 | 0.119 |                0.059 |      0.542 |        0     |        0.458 |
| loso     | bandgap   | C2/m      | cvplus_lattice | -0.367 | 0.378 |                0.156 |      0.997 |        0     |        0.003 |
| loso     | formation | C2/m      | cvplus_lattice |  0.093 | 0.1   |                0.044 |      0.997 |        0.003 |        0     |
| loso     | bandgap   | Pna21     | cvplus_lattice |  0     | 0.118 |                0.118 |      1     |        0     |        0     |
| loso     | formation | Pna21     | cvplus_lattice |  0.031 | 0.036 |                0.025 |      1     |        0     |        0     |
| loso     | bandgap   | R-3c      | cvplus_lattice | -0.384 | 0.394 |                0.183 |      0.974 |        0     |        0.026 |
| loso     | formation | R-3c      | cvplus_lattice |  0.053 | 0.057 |                0.032 |      1     |        0     |        0     |
| loso     | bandgap   | P63/mmc   | cvplus_lattice |  0.127 | 0.256 |                0.246 |      0.957 |        0.043 |        0     |
| loso     | formation | P63/mmc   | cvplus_lattice | -0.064 | 0.096 |                0.07  |      0.915 |        0.003 |        0.083 |
| loso     | bandgap   | Ia-3      | cvplus_lattice | -0.251 | 0.261 |                0.118 |      0.992 |        0.002 |        0.006 |
| loso     | formation | Ia-3      | cvplus_lattice |  0.036 | 0.042 |                0.043 |      0.994 |        0.002 |        0.004 |
| l2o      | bandgap   | C2/m      | cvplus_lattice | -0.026 | 0.173 |                0.169 |      0.994 |        0.003 |        0.003 |
| l2o      | formation | C2/m      | cvplus_lattice |  0.03  | 0.059 |                0.048 |      0.997 |        0     |        0.003 |
| l2o      | bandgap   | Pna21     | cvplus_lattice | -0.051 | 0.108 |                0.101 |      1     |        0     |        0     |
| l2o      | formation | Pna21     | cvplus_lattice |  0.009 | 0.032 |                0.031 |      1     |        0     |        0     |
| l2o      | bandgap   | R-3c      | cvplus_lattice | -0.18  | 0.241 |                0.188 |      0.954 |        0     |        0.046 |
| l2o      | formation | R-3c      | cvplus_lattice |  0.042 | 0.056 |                0.045 |      1     |        0     |        0     |
| l2o      | bandgap   | P63/mmc   | cvplus_lattice |  0.427 | 0.442 |                0.235 |      0.33  |        0.67  |        0     |
| l2o      | formation | P63/mmc   | cvplus_lattice | -0.099 | 0.115 |                0.064 |      0.299 |        0.034 |        0.667 |
| l2o      | bandgap   | Ia-3      | cvplus_lattice | -0.038 | 0.094 |                0.095 |      0.987 |        0.004 |        0.008 |
| l2o      | formation | Ia-3      | cvplus_lattice |  0.021 | 0.041 |                0.044 |      0.996 |        0     |        0.004 |

## Pre-registered decisions

```
{
  "bandgap": {
    "H4_supported": false,
    "H4_detail": {
      "diff": 0.3544,
      "diff_lo": 0.295,
      "diff_hi": 0.4166,
      "best_distance": 0.4342
    },
    "G4_triggered": true,
    "pooled_auroc": {
      "width": 0.08,
      "knn": 0.434,
      "kde": 0.404,
      "committee": 0.15,
      "subsample": 0.181,
      "classifier": 0.549
    },
    "within_weighted_auroc": {
      "width": 0.534,
      "knn": 0.404,
      "kde": 0.401,
      "committee": 0.563,
      "subsample": 0.476,
      "classifier": 0.546
    }
  },
  "formation": {
    "H4_supported": false,
    "H4_detail": {
      "diff": 0.3475,
      "diff_lo": 0.2948,
      "diff_hi": 0.4128,
      "best_distance": 0.4295
    },
    "G4_triggered": true,
    "pooled_auroc": {
      "width": 0.082,
      "knn": 0.429,
      "kde": 0.403,
      "committee": 0.386,
      "subsample": 0.445,
      "classifier": 0.62
    },
    "within_weighted_auroc": {
      "width": 0.388,
      "knn": 0.448,
      "kde": 0.465,
      "committee": 0.437,
      "subsample": 0.581,
      "classifier": 0.594
    }
  }
}
```
