# Pre-registration: Phase 4 (detecting failure before it matters)

Registered (UTC): 2026-11-17T19:13:28Z
Config hash: f654359b5ecc268096513ed8e8217e40cffcf85e37d8a410ea8e777265903ce5
Code hash:   e411f88d5433b979fcd8fe02213244bf0e30a8b00b9d65750d84ffd581c000d9
Phase 2 inputs: phase0/phase2_inputs.json (SHA-256 of all 22 units)

Written once, before any Phase 4 score is computed; never overwritten. A later change to the
scientific settings is refused unless explicitly allowed and is then logged in DEVIATIONS.json.

## What was known when this was written
Phase 1, 2 and 3 results, including: LOSO error is mostly a lattice-wide offset; a domain
classifier separates every lattice from the rest (AUROC 1.000); for LOSO-Fd-3m band gap the
median kNN distance of the held-out structures (97-109) lies inside the range of the out-of-fold
lattices (85-131). No Phase 4 score, AUROC or flag rate had been computed.

## Data and configurations
n = 2,343 (Phase 2 de-duplication). 6 nested-LOSO holdouts and 5 leave-two-out holdouts, using
the Phase 2 models without refitting. Primary analysis: the 6 LOSO holdouts pooled (every
structure appears once).

## Outcomes (per structure, per target)
* Primary: miscovered by the reference interval, cvplus_lattice at alpha = 0.1 (the submitted design).
* Secondary: |error| above the in-distribution 0.9-quantile (random-fold out-of-fold |residuals|
  on the same training lattices).

## Scores (higher = more doubt)
width (cvplus_lattice interval width); knn (mean distance to the 25 nearest training structures,
standardized SOAP, PCA-50); kde (Schultz et al. 2025: Epanechnikov KDE on all standardized
features, bandwidth by the nearest-neighbour rule with quantile 0.3, d = 1 - KDE/max_train KDE);
committee (SD of the lattice-fold models); subsample (SD of the random-fold models; LOSO only);
classifier (out-of-fold probability of belonging to the screened set; logistic, PCA-50).

## Analyses
1. Pooled structure-level AUROC and average precision, 95 % cluster-bootstrap CI (2000 resamples of
   composition groups).
2. Within-lattice AUROC for every holdout with >= 10 structures in each class, and their
   n-weighted mean.
3. Selective prediction: interval coverage and MAE of the retained structures when abstaining
   on the highest-score fraction (grid 0.0-0.9).
4. Lattice level (11 configurations, descriptive): Spearman correlation of the mean score with the
   miscoverage rate and with |bias|, permutation p-value (100000 permutations). The 11
   configurations reuse 6 test lattices, so they are not independent.
5. Flag: kNN-distance conformal p-values (Bates et al. 2023), flag at p <= 0.05; report the
   in-distribution false-flag rate and the flagged fraction of every held-out lattice.

## Decision rules (primary outcome, each target separately)
* H4 supported if max(AUROC_knn, AUROC_kde) - AUROC_width >= 0.05, with the 95 % cluster-
  bootstrap CI of that difference (max taken inside each resample) above 0, and
  max(AUROC_knn, AUROC_kde) >= 0.6. Otherwise H4 is not supported.
* G4 triggered if every score has pooled AUROC < 0.7 AND every score's n-weighted within-lattice
  AUROC is < 0.7. Then the recommendation narrows to: obtain labelled calibration data for every
  family screened, and Phase 5 carries the practical message.
* Expected for the flag (stated now): false-flag rate near 0.05 and close to 100 % of every held-out
  lattice flagged, i.e. the flag detects "new family", not "will fail". The paper reports it as such.
* Every score, outcome, lattice and configuration is reported whichever way it falls.
