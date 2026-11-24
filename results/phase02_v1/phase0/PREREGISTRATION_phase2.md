# Pre-registration: Phase 2 (LOSO mechanism) and Phase 3 (overlap diagnostics)

Registered (UTC): 2026-11-17T07:14:00Z
Config hash: 8b86e24ced22b6dce03f797c2c396f0a101c7ab5a8a95b6a3162dc8aa91badd1
Code hash:   88eab1e03e62226ca0a8c3cd93befe4d28f16249321070dd21aeb23f8ec1d382

Written once, before any Phase 2 or Phase 3 result exists; never overwritten. A later change
to the scientific settings is refused unless explicitly allowed and is then logged in
DEVIATIONS.json next to this file. Phase 1 results were known when this was written; the
Phase 1 follow-up checks (step post1) are therefore labelled exploratory, not confirmatory.

## Data
* NOMAD 2018 training set. v2 de-duplication: same input structure (pymatgen StructureMatcher,
  ltol = stol = 0.001, angle_tol = 0.05 deg, on the Vegard geometry files) AND both targets within
  0.002 eV(/atom) -> one record (lowest id kept). n = 2343 entries after de-duplication.
  The tolerance was chosen after inspecting the distribution of target differences inside
  identical-input clusters (largest within-tolerance difference 0.0015 eV; every conflicting cluster
  differs by >= 0.0125 eV(/atom) in at least one target),
  before any Phase 2 model was fitted.
* Same-input entries whose targets differ by more (16 entries in 8 clusters) are
  kept and flagged; every Phase 2 coverage is also reported with them excluded from evaluation.

## Design (nested leave-one-lattice-out, LOSO)
* For each held-out lattice g (order ['Fd-3m', 'C2/m', 'Pna21', 'R-3c', 'P63/mmc', 'Ia-3']) and target: hyperparameters tuned with Optuna TPE
  (15 trials, 5 random start-up, MAE, early stopping 50, max 1000 rounds) using
  leave-one-lattice-out inner folds over the five training lattices only (g is never seen).
* CV+ with the same five lattice folds (cvplus_lattice). Nested scores with the same folds
  (Gupta et al. 2022): normalized by a kNN distance sigma fitted on each fold's training lattices
  (ncvplus_knn; k = 25, 50 PCA components) and CQR (cqr_cvplus, alpha = 0.1).
  Contrast: CV+ with 5 composition-grouped random folds on the same training lattices
  (cvplus_random_in): in-distribution residuals only.
* Calibration-only ablation (no refit): the same fold models and residuals with the Fd-3m fold
  removed from the calibration terms (suffix _noFd3m_cal), for every g other than Fd-3m.
* Leave-two-out refit (Task 2.2): for every X other than Fd-3m, train and calibrate on the four
  lattices other than X and Fd-3m (4 lattice folds, hyperparameters of LOSO-X), test on X.
* Group-level conformal (Task 2.4, illustrative): family score = 0.9-quantile of a lattice's
  out-of-fold residuals; threshold = ceil((K+1)(1-alpha))-th smallest of K = 5 family scores
  (Dunn et al. 2023). alpha = 0.1 gives an unbounded interval by construction; alpha = 0.2 the
  largest of the five.
* alpha in [0.05, 0.1, 0.2]; primary alpha = 0.1.

## Phase 3 overlap gate (fixed now; diagnostics are computed before LOSO results)
* Classifiers: held-out lattice vs the other five, logistic regression (C = 1.0) and kNN (k = 25)
  on 50 PCA components of the standardized SOAP vector, 5-fold composition-grouped
  out-of-fold probabilities. AUROC of each; density-ratio weights w = p/(1-p) * n_rest/n_g
  for the five calibration lattices; Kish ESS = (sum w)^2 / sum w^2; test-point weight share
  p_test = w(x) / (sum_i w_i + w(x)) (weighted CP is unbounded when p_test > alpha).
* Weighted CV+ is run for holdout g only if max(AUROC) <= 0.9 AND ESS >= 50 AND the
  median p_test < 0.1. Otherwise weighted CP is reported as inapplicable for that unseen
  lattice, with these numbers.

## Decision rules (alpha = 0.1, each target separately)
* P2.1 (Task 2.1 prediction): for every holdout g other than Fd-3m, Fd-3m's share of the
  calibration residuals at or above the pooled conformal quantile, divided by its share of
  calibration points, is >= 2. Supported if this holds in >= 4 of the 5 holdouts.
* H2 (primary contrast = calibration-only ablation, cvplus_lattice): for X other than Fd-3m let
  dcov = coverage(ablation) - coverage(LOSO) and r = median width(ablation) / median width(LOSO).
  H2 supported if r <= 0.80 and dcov <= -0.05 in >= 4 of 5 X.
  H2 rejected (gate G2) if r >= 0.90 and |dcov| <= 0.03 in >= 4 of 5 X. Otherwise partial.
  The same quantities for the leave-two-out refit are reported as the plan's Task 2.2 result.
* G3: positive if, under LOSO-Fd-3m, ncvplus_knn or cqr_cvplus covers >= 0.85 of Fd-3m.
* Task 2.3 (descriptive): a sub-family 'carries the failure' if its LOSO-Fd-3m cvplus_lattice
  coverage is < 0.75 while the other's is >= 0.85.
* Every lattice, method and alpha is reported whichever way it falls.
