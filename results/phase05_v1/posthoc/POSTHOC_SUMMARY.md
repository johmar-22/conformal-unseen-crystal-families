# Phase 5 post hoc analyses (NOT pre-registered)

Designed after the Phase 5 results were seen. m = 20 labels, one per composition, alpha = 0.1; 500 draws per cell, the SAME label sets as Phase 5 (re-created with the registered seeds and checked draw by draw against the saved draws).

* icc_resid: ICC(1) of the residuals by composition group; icc_absdev: ICC(1) of |z - mean z|.
* rho_clustered, cov_clustered_split: pre-registered H5c quantities (split, clustered labels).
* cov_clustered_offset: mean coverage of new compositions, offset-corrected interval, clustered labels.
* width_*_sd: median interval width / family SD. nomodel = full conformal on the labelled y values (no model). ratio = median width offset / median width no-model.

| family                        | target    |   icc_resid |   icc_absdev |   rho_clustered |   cov_clustered_split |   cov_clustered_offset |   width_offset_sd |   width_nomodel_sd |   ratio_offset_to_nomodel |   cov_offset |   cov_nomodel |
|:------------------------------|:----------|------------:|-------------:|----------------:|----------------------:|-----------------------:|------------------:|-------------------:|--------------------------:|-------------:|--------------:|
| Fd-3m (LOSO)                  | bandgap   |       0.033 |        0.036 |           1.074 |                 0.898 |                  0.885 |             1.886 |              2.965 |                     0.636 |        0.898 |         0.904 |
| Fd-3m (LOSO)                  | formation |       0.131 |        0.012 |           1.198 |                 0.915 |                  0.88  |             2.623 |              3.345 |                     0.784 |        0.905 |         0.908 |
| C2/m (LOSO)                   | bandgap   |       0.507 |        0.377 |           1.604 |                 0.841 |                  0.799 |             1.405 |              4.314 |                     0.326 |        0.907 |         0.906 |
| C2/m (LOSO)                   | formation |       0.493 |        0.501 |           1.561 |                 0.866 |                  0.815 |             2.815 |              3.666 |                     0.768 |        0.909 |         0.906 |
| Pna21 (LOSO)                  | bandgap   |       0.554 |        0.372 |           1.695 |                 0.849 |                  0.8   |             0.614 |              3.175 |                     0.193 |        0.904 |         0.902 |
| Pna21 (LOSO)                  | formation |       0.49  |        0.112 |           1.692 |                 0.883 |                  0.832 |             1.459 |              3.442 |                     0.424 |        0.904 |         0.908 |
| R-3c (LOSO)                   | bandgap   |       0.79  |        0.678 |           2.228 |                 0.847 |                  0.71  |             0.845 |              3.16  |                     0.267 |        0.902 |         0.899 |
| R-3c (LOSO)                   | formation |       0.713 |        0.463 |           1.797 |                 0.85  |                  0.737 |             1.884 |              3.216 |                     0.586 |        0.91  |         0.908 |
| P63/mmc (LOSO)                | bandgap   |       0.263 |        0.228 |           1.634 |                 0.879 |                  0.845 |             1.513 |              3.41  |                     0.444 |        0.908 |         0.903 |
| P63/mmc (LOSO)                | formation |       0.445 |        0.424 |           1.05  |                 0.897 |                  0.851 |             3.613 |              3.53  |                     1.023 |        0.903 |         0.902 |
| Ia-3 (LOSO)                   | bandgap   |       0.613 |        0.471 |           2.72  |                 0.809 |                  0.755 |             0.503 |              3.221 |                     0.156 |        0.903 |         0.902 |
| Ia-3 (LOSO)                   | formation |       0.812 |        0.586 |           3.001 |                 0.761 |                  0.629 |             2.398 |              3.151 |                     0.761 |        0.905 |         0.904 |
| P63/mmc (refit without Fd-3m) | bandgap   |       0.252 |        0.277 |           1.717 |                 0.879 |                  0.856 |             1.372 |              3.377 |                     0.406 |        0.907 |         0.897 |
| P63/mmc (refit without Fd-3m) | formation |       0.455 |        0.35  |           1.061 |                 0.898 |                  0.849 |             3.11  |              3.537 |                     0.879 |        0.903 |         0.901 |

Spearman ICC(resid) vs clustered dispersion ratio: rho = 0.73, permutation p = 0.0043 (asymptotic 0.0032), n = 14 family-target pairs.
Spearman ICC(|dev|) vs clustered dispersion ratio: rho = 0.53, permutation p = 0.0565.

Verification: model widths and coverages of every re-created draw equal the Phase 5 draws (max |dwidth| 4.4e-16, max |dcov| 1.1e-16).
