# Phase 6 Task 6.5 post hoc analyses (v8.1; NOT pre-registered)

Designed after the v8 results were seen (2026-11-19). Nothing here changes a pre-registered decision.

## 1. Few-shot intervals against a no-model interval from the same labels

Label sets re-created with the registered fewshot8 seeds and checked draw by draw against the saved draws (model widths and coverages identical). No-model interval = offset-corrected full conformal on the labelled y values (centred on their mean; the posthoc5 definition). alpha = 0.1; 500 draws per cell; one structure per composition. ratio < 1: the model interval is narrower.

| target    | symbol   |   unseen_coverage |   m |   cov_mean_split |   cov_q05_split |   cov_q95_split |   p_below_split |   width_sd_split |   width_sd_offset_full |   width_sd_nomodel |   ratio_offset_to_nomodel |   ratio_split_to_nomodel |   share_draws_model_narrower |   cov_mean_nomodel |
|:----------|:---------|------------------:|----:|-----------------:|----------------:|----------------:|----------------:|-----------------:|-----------------------:|-------------------:|--------------------------:|-------------------------:|-----------------------------:|-------------------:|
| formation | F-43m    |             0.582 |  10 |            0.913 |           0.754 |           0.996 |           0.096 |            2.505 |                  1.942 |              3.647 |                     0.532 |                    0.687 |                        0.954 |              0.912 |
| formation | F-43m    |             0.582 |  20 |            0.899 |           0.776 |           0.978 |           0.082 |            2.327 |                  1.658 |              3.352 |                     0.495 |                    0.694 |                        0.996 |              0.908 |
| formation | F-43m    |             0.582 |  30 |            0.906 |           0.819 |           0.972 |           0.026 |            2.331 |                  1.629 |              3.196 |                     0.510 |                    0.729 |                        1.000 |              0.901 |
| formation | Pm       |             0.704 |  10 |            0.912 |           0.743 |           0.995 |           0.122 |            2.156 |                  1.740 |              3.855 |                     0.451 |                    0.559 |                        0.990 |              0.911 |
| formation | Pm       |             0.704 |  20 |            0.908 |           0.775 |           0.983 |           0.070 |            1.889 |                  1.448 |              3.568 |                     0.406 |                    0.529 |                        1.000 |              0.907 |
| formation | Pm       |             0.704 |  30 |            0.905 |           0.798 |           0.976 |           0.058 |            1.743 |                  1.321 |              3.437 |                     0.384 |                    0.507 |                        1.000 |              0.902 |
| formation | P-6m2    |             0.771 |  10 |            0.913 |           0.762 |           0.995 |           0.084 |            2.065 |                  1.784 |              4.262 |                     0.418 |                    0.485 |                        0.922 |              0.914 |
| formation | P-6m2    |             0.771 |  20 |            0.906 |           0.775 |           0.980 |           0.068 |            1.839 |                  1.547 |              3.534 |                     0.438 |                    0.520 |                        0.976 |              0.905 |
| formation | P-6m2    |             0.771 |  30 |            0.906 |           0.814 |           0.973 |           0.036 |            1.800 |                  1.470 |              3.371 |                     0.436 |                    0.534 |                        0.998 |              0.906 |
| bandgap   | P4/nmm   |             0.784 |  10 |            0.914 |           0.729 |           0.997 |           0.110 |            2.134 |                  2.337 |              3.727 |                     0.627 |                    0.572 |                        0.936 |              0.902 |
| bandgap   | P4/nmm   |             0.784 |  20 |            0.907 |           0.792 |           0.983 |           0.066 |            1.950 |                  2.031 |              3.330 |                     0.610 |                    0.586 |                        0.964 |              0.905 |
| bandgap   | P4/nmm   |             0.784 |  30 |            0.903 |           0.799 |           0.972 |           0.054 |            1.929 |                  1.952 |              3.284 |                     0.594 |                    0.587 |                        0.996 |              0.907 |
| bandgap   | P4mm     |             0.789 |  10 |            0.908 |           0.743 |           0.996 |           0.124 |            3.641 |                  3.613 |              4.056 |                     0.891 |                    0.898 |                        0.554 |              0.916 |
| bandgap   | P4mm     |             0.789 |  20 |            0.904 |           0.780 |           0.981 |           0.084 |            3.201 |                  3.258 |              2.709 |                     1.203 |                    1.182 |                        0.386 |              0.906 |
| bandgap   | P4mm     |             0.789 |  30 |            0.906 |           0.811 |           0.980 |           0.038 |            3.201 |                  3.103 |              2.558 |                     1.213 |                    1.251 |                        0.344 |              0.903 |
| bandgap   | Fm-3m    |             0.805 |  10 |            0.908 |           0.754 |           0.995 |           0.106 |            1.749 |                  1.719 |              3.847 |                     0.447 |                    0.455 |                        0.944 |              0.913 |
| bandgap   | Fm-3m    |             0.805 |  20 |            0.904 |           0.784 |           0.986 |           0.078 |            1.657 |                  1.596 |              3.241 |                     0.492 |                    0.511 |                        0.986 |              0.909 |
| bandgap   | Fm-3m    |             0.805 |  30 |            0.906 |           0.810 |           0.971 |           0.034 |            1.646 |                  1.550 |              2.992 |                     0.518 |                    0.550 |                        0.998 |              0.905 |

## 2. Inner-model SD (MatFold's ensemble SD) at group level

Spearman of repeat-1 group coverage with the group median of each score; percentile bootstrap over groups (2000) and Monte Carlo permutation p (100000).

| target    | score            |    rho |   ci_lo |   ci_hi |   n_boot |   p_perm |   n_groups |
|:----------|:-----------------|-------:|--------:|--------:|---------:|---------:|-----------:|
| formation | committee_median | -0.538 |  -0.721 |  -0.299 |     2000 |    0.000 |         59 |
| formation | width_median     |  0.195 |  -0.070 |   0.436 |     2000 |    0.138 |         59 |
| formation | knn_median       | -0.246 |  -0.522 |   0.028 |     2000 |    0.061 |         59 |
| bandgap   | committee_median | -0.552 |  -0.811 |  -0.176 |     2000 |    0.002 |         31 |
| bandgap   | width_median     |  0.277 |  -0.070 |   0.602 |     2000 |    0.131 |         31 |
| bandgap   | knn_median       | -0.319 |  -0.591 |   0.066 |     2000 |    0.081 |         31 |

Coverage retained inside the three worst groups when the structures with the highest SD (pooled threshold over all held-out structures of the target) are abstained:

| target    | symbol   |   abstain |    n |   share_abstained |   retained_coverage |
|:----------|:---------|----------:|-----:|------------------:|--------------------:|
| formation | F-43m    |     0.000 | 4533 |             0.000 |               0.590 |
| formation | Pm       |     0.000 |  790 |             0.000 |               0.680 |
| formation | P-6m2    |     0.000 | 3487 |             0.000 |               0.771 |
| formation | F-43m    |     0.100 | 4533 |             0.136 |               0.650 |
| formation | Pm       |     0.100 |  790 |             0.400 |               0.926 |
| formation | P-6m2    |     0.100 | 3487 |             0.299 |               0.848 |
| formation | F-43m    |     0.250 | 4533 |             0.283 |               0.715 |
| formation | Pm       |     0.250 |  790 |             0.499 |               0.975 |
| formation | P-6m2    |     0.250 | 3487 |             0.477 |               0.900 |
| formation | F-43m    |     0.500 | 4533 |             0.537 |               0.788 |
| formation | Pm       |     0.500 |  790 |             0.620 |               0.983 |
| formation | P-6m2    |     0.500 | 3487 |             0.652 |               0.946 |
| bandgap   | P4/nmm   |     0.000 |  313 |             0.000 |               0.780 |
| bandgap   | P4mm     |     0.000 |  279 |             0.000 |               0.781 |
| bandgap   | Fm-3m    |     0.000 | 1339 |             0.000 |               0.812 |
| bandgap   | P4/nmm   |     0.100 |  313 |             0.093 |               0.803 |
| bandgap   | P4mm     |     0.100 |  279 |             0.247 |               0.795 |
| bandgap   | Fm-3m    |     0.100 | 1339 |             0.311 |               0.842 |
| bandgap   | P4/nmm   |     0.250 |  313 |             0.259 |               0.823 |
| bandgap   | P4mm     |     0.250 |  279 |             0.505 |               0.848 |
| bandgap   | Fm-3m    |     0.250 | 1339 |             0.571 |               0.859 |
| bandgap   | P4/nmm   |     0.500 |  313 |             0.553 |               0.850 |
| bandgap   | P4mm     |     0.500 |  279 |             0.753 |               0.870 |
| bandgap   | Fm-3m    |     0.500 | 1339 |             0.807 |               0.903 |

## 3. Oxides: ABO3 compositions (composition rule only) and other oxides

Pooled over all eligible held-out groups (mean over repeats): formation: ABO3 0.870 (n 1192), other oxides 0.975 (n 10524), non-oxides 0.909 (n 63481); bandgap: ABO3 0.904 (n 491), other oxides 0.934 (n 5311), non-oxides 0.900 (n 10100); tbmbj: ABO3 0.835 (n 139), other oxides 0.916 (n 679), non-oxides 0.903 (n 2265).

| target    |   group | symbol   |   n_abo3 |   coverage_abo3 |   n_other_oxide |   coverage_other_oxide |   coverage_oxide |
|:----------|--------:|:---------|---------:|----------------:|----------------:|-----------------------:|-----------------:|
| formation |     221 | Pm-3m    |      493 |           0.732 |              43 |                  0.895 |            0.745 |
| formation |      99 | P4mm     |       14 |           1.000 |             257 |                  0.778 |            0.790 |
| formation |     216 | F-43m    |        0 |         nan     |              66 |                  0.841 |            0.841 |
| formation |      44 | Imm2     |        0 |         nan     |              54 |                  0.898 |            0.898 |
| formation |      47 | Pmmm     |        0 |         nan     |              68 |                  0.904 |            0.904 |
| formation |     123 | P4/mmm   |        2 |           0.500 |             284 |                  0.917 |            0.914 |
| formation |     194 | P6_3/mmc |       46 |           0.957 |             127 |                  0.921 |            0.931 |
| formation |     225 | Fm-3m    |        0 |         nan     |             341 |                  0.952 |            0.952 |
| bandgap   |     221 | Pm-3m    |       97 |           0.732 |              19 |                  0.789 |            0.741 |
| bandgap   |      99 | P4mm     |        4 |           1.000 |              50 |                  0.740 |            0.759 |
| bandgap   |     225 | Fm-3m    |        0 |         nan     |             186 |                  0.836 |            0.836 |
| bandgap   |     123 | P4/mmm   |        0 |         nan     |              63 |                  0.873 |            0.873 |
| bandgap   |      74 | Imma     |        9 |           1.000 |              98 |                  0.878 |            0.888 |
| bandgap   |       8 | Cm       |        7 |           0.929 |             144 |                  0.896 |            0.897 |
| bandgap   |      11 | P2_1/m   |       14 |           1.000 |             265 |                  0.906 |            0.910 |
| bandgap   |     148 | R-3      |       47 |           0.936 |             155 |                  0.903 |            0.911 |
| tbmbj     |     221 | Pm-3m    |       86 |           0.826 |              13 |                  0.923 |            0.838 |
| tbmbj     |       2 | P-1      |        3 |           1.000 |              86 |                  0.837 |            0.843 |
| tbmbj     |      62 | Pnma     |       36 |           0.806 |              36 |                  0.889 |            0.847 |
| tbmbj     |     225 | Fm-3m    |        0 |         nan     |              68 |                  0.868 |            0.868 |
| tbmbj     |     166 | R-3m     |        0 |         nan     |              88 |                  0.898 |            0.898 |
| tbmbj     |     164 | P-3m1    |        0 |         nan     |              61 |                  0.934 |            0.934 |
| tbmbj     |      14 | P2_1/c   |        6 |           0.833 |              89 |                  0.955 |            0.947 |
| tbmbj     |      12 | C2/m     |        2 |           1.000 |             115 |                  0.948 |            0.949 |

## 4. Degenerate group-level intervals

atom = 1/(K+1) at +infinity; unbounded if atom > alpha, upper end = largest calibration term if atom = alpha. Single subsampling: rank ceil((K+1)(1-alpha)) vs K.

| target    | method                  |   alpha |   K_min |   K_max |    atom | status       |
|:----------|:------------------------|--------:|--------:|--------:|--------:|:-------------|
| formation | hierarchical jackknife+ |   0.050 |      53 |      54 |   0.019 | finite       |
| formation | group-weighted CV+      |   0.050 |      53 |      54 |   0.019 | finite       |
| formation | hierarchical jackknife+ |   0.100 |      53 |      54 |   0.019 | finite       |
| formation | group-weighted CV+      |   0.100 |      53 |      54 |   0.019 | finite       |
| formation | hierarchical jackknife+ |   0.200 |      53 |      54 |   0.019 | finite       |
| formation | group-weighted CV+      |   0.200 |      53 |      54 |   0.019 | finite       |
| formation | H3 HCP                  |   0.050 |      20 |      20 |   0.048 | finite       |
| formation | H3 single subsampling   |   0.050 |      20 |      20 | nan     | largest term |
| formation | H3 HCP                  |   0.100 |      20 |      20 |   0.048 | finite       |
| formation | H3 single subsampling   |   0.100 |      20 |      20 | nan     | finite       |
| formation | H3 HCP                  |   0.200 |      20 |      20 |   0.048 | finite       |
| formation | H3 single subsampling   |   0.200 |      20 |      20 | nan     | finite       |
| bandgap   | hierarchical jackknife+ |   0.050 |      27 |      28 |   0.036 | finite       |
| bandgap   | group-weighted CV+      |   0.050 |      27 |      28 |   0.036 | finite       |
| bandgap   | hierarchical jackknife+ |   0.100 |      27 |      28 |   0.036 | finite       |
| bandgap   | group-weighted CV+      |   0.100 |      27 |      28 |   0.036 | finite       |
| bandgap   | hierarchical jackknife+ |   0.200 |      27 |      28 |   0.036 | finite       |
| bandgap   | group-weighted CV+      |   0.200 |      27 |      28 |   0.036 | finite       |
| bandgap   | H3 HCP                  |   0.050 |      11 |      11 |   0.083 | unbounded    |
| bandgap   | H3 single subsampling   |   0.050 |      11 |      11 | nan     | unbounded    |
| bandgap   | H3 HCP                  |   0.100 |      11 |      11 |   0.083 | finite       |
| bandgap   | H3 single subsampling   |   0.100 |      11 |      11 | nan     | largest term |
| bandgap   | H3 HCP                  |   0.200 |      11 |      11 |   0.083 | finite       |
| bandgap   | H3 single subsampling   |   0.200 |      11 |      11 | nan     | finite       |
| tbmbj     | hierarchical jackknife+ |   0.050 |       9 |       9 |   0.100 | unbounded    |
| tbmbj     | group-weighted CV+      |   0.050 |       9 |       9 |   0.100 | unbounded    |
| tbmbj     | hierarchical jackknife+ |   0.100 |       9 |       9 |   0.100 | largest term |
| tbmbj     | group-weighted CV+      |   0.100 |       9 |       9 |   0.100 | largest term |
| tbmbj     | hierarchical jackknife+ |   0.200 |       9 |       9 |   0.100 | finite       |
| tbmbj     | group-weighted CV+      |   0.200 |       9 |       9 |   0.100 | finite       |
