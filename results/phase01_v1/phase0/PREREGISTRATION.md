# Pre-registration: Phase 1 (in-distribution baseline)

Registered (UTC): 2026-11-16T08:15:53Z
Config hash: 0e3daf052c148e71abe7ee9d05813f31ab2ac6b40247023980eb19cbfc226275
Code hash:   fa22ab7a0f419bfc8b0bd1b2406473799891e1b888c35ffd84056d8a91d47463
Config file: /content/cpaudit_v1/config/phase01.yaml

This file is written once and never overwritten. Any later change to the scientific
settings is refused by the runner unless explicitly allowed, and is then logged in
DEVIATIONS.json next to this file.

## Design
* Data: NOMAD 2018 public training set (2,400 entries), exact duplicates removed
  (identical tabular record and identical SOAP vector; lowest id kept).
* Unit of exchangeability: composition group (lattice, Al and Ga fraction rounded to 0.1).
  Outer test set, calibration set and all CV folds are group-disjoint and stratified by lattice.
* 25 independent outer splits (test = 20 %), calibration = 25 % of the training pool.
* Model: LightGBM (L2 objective), hyperparameters tuned with Optuna TPE (25 trials, MAE,
  5-fold composition-grouped CV) on the proper-training set of split 0 (mode = once);
  splits [1, 2] are also re-tuned on their own proper-training set to measure sensitivity.
* Interval methods: cvplus_random (K=5), cvplus_lattice (one fold per lattice, the submitted
  design, reference only), split_global, split_mondrian (per lattice), split_normalized_knn
  (k=25, 50 PCA components), cqr (alpha = 0.1 only).
* alpha in [0.05, 0.1, 0.2]; primary alpha = 0.1.

## Hypotheses (Phase 1 part of H1)
H1a  Under a valid in-distribution calibration (cvplus_random and split_global), per-lattice
     coverage is heterogeneous and Fd-3m is under-covered.
H1b  Mondrian (per-lattice) calibration restores per-lattice coverage to nominal, at a width cost.
H1c  (descriptive) The lattice-fold CV+ design of the submitted manuscript over-covers every
     lattice relative to cvplus_random, on identical data.

## Primary endpoint
Mean over the 25 splits of the band-gap coverage of Fd-3m test points under cvplus_random at
alpha = 0.1. Same endpoint for formation energy is co-primary.

## Decision rules (fixed before any Phase 1 result is seen)
* H1a supported   if the primary endpoint < 0.85.
* H1a rejected    if the primary endpoint >= 0.88  -> decision gate G1: the in-distribution
                  heterogeneity claim is dropped and reported as a negative result.
* 0.85-0.88       inconclusive; reported as such, no claim made.
* H1b supported   if, under split_mondrian at alpha = 0.1, the mean coverage of every lattice
                  is >= 0.88 and the observed fraction of splits below 0.90 is not larger than
                  the Beta-binomial expectation for a valid procedure + 0.15.
* Every endpoint, method, alpha and lattice is reported whichever way it falls.

## Sampling-noise reference
For split-conformal methods the expected fraction of splits with coverage below 1 - alpha under
exact validity is computed per split from Beta-binomial(n_test, k, n_cal + 1 - k) with
k = ceil((1 - alpha)(n_cal + 1)) and stored as p_below_if_valid. Because points within a
composition group are correlated, the true split-to-split spread is larger than this
independent-point reference; it is therefore a lower bound on the noise, and the number of
test and calibration composition groups is reported alongside every coverage value.
