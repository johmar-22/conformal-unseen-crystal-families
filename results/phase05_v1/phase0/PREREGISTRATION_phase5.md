# Pre-registration: Phase 5 (few-shot, family-specific recalibration)

Registered (UTC): 2026-11-17T20:12:01Z
Config hash: db2e1156ab5e3c6f00937c6f712d6782e0f2d8b48cd1397185b45e1ab635a2ef
Code hash:   fbebadfd10993cf98195f21368a5e719445ceba660654a4bfa8ac86228453e6f
Phase 2 inputs: phase0/phase2_inputs.json (SHA-256 of every file used)

Written once, before any Phase 5 quantity is computed on the real data; never overwritten. A
later change to the scientific settings is refused unless explicitly allowed, and is then
logged in phase0/DEVIATIONS.json.

## What was known when this was written
All Phase 1-4 results. In particular: under LOSO the error of a held-out lattice is mostly one
offset for the whole lattice (|bias|/MAE >= 0.84 in 9 of the 12 LOSO lattice-target pairs; band gap
P63/mmc 0.49 and Pna21 0.00, formation P63/mmc 0.66; leave-two-out P63/mmc 0.97 and 0.86);
the Fd-3m cell-setting sub-families have offsets +0.625/+0.673 eV (band gap) and -0.099/-0.147
eV/atom (formation); no score detects failure (Phase 4, G4 triggered); composition groups, not
structures, are the exchangeability unit in-distribution (post1). No few-shot coverage, width or
dispersion had been computed on the real data; the code was tested on synthetic data only.

## Question
How many labelled structures of a new family restore valid, useful intervals, and does removing
the family offset make them narrower?

## Families and model
7 families: Fd-3m (LOSO); C2/m (LOSO); Pna21 (LOSO); R-3c (LOSO); P63/mmc (LOSO); Ia-3 (LOSO); P63/mmc (refit without Fd-3m). For each, the fixed Phase 2 model mu = mean of the five lattice-fold
models (the predictor whose CV+ intervals were evaluated in Phases 2 and 4), trained without the
family; nothing is refitted. Residual z = y - mu. Targets: band gap and formation energy.
alpha = 0.1.

## Methods (each draw uses the same m labelled structures)
* split: split conformal with |z| scores, k = ceil((1 - alpha)(m + 1)).
* offset_full: full conformal for y = mu + b + e with b the mean residual of the augmented data
  (exact marginal coverage k/(m + 1) under exchangeability; closed form in cpaudit/fewshot.py).
  For this model the full-conformal set is always one interval and equals the jackknife+
  interval with the leave-one-out mean offset (Barber et al. 2021); the identity is verified in
  the self-test, so the exact full-conformal guarantee applies (not only jackknife+'s 1 - 2 alpha).

## Sampling schemes (500 draws per scheme and m)
* random: m random structures, m in [9, 10, 20, 30, 50, 100]; coverage of all other structures of the family.
  Reference: exact finite-population law (negative hypergeometric; mean k/(m + 1), tends to
  Beta(k, m + 1 - k)).
* clustered: whole composition groups in random order until m structures, m in [9, 10, 20, 30, 50, 100] (as when
  several configurations of one composition are computed).
* one_per_group: m composition groups, one structure each (single subsampling, Dunn et al. 2023),
  m in [9, 10, 20, 30, 40].
  For clustered and one_per_group, at least 10 composition groups must remain untouched.
  Primary coverage = cov_new: mean over the untouched groups of the covered fraction (new
  compositions, the screening target). Also recorded: cov_family, the covered fraction of all
  non-calibration structures.
* Permutation null for clustered and one_per_group: the same procedure after permuting the
  residuals across the family's structures (same group sizes, no group structure).
* sub_a2b / sub_b2a (Fd-3m only, exploratory): calibrate on one cell setting (227a, gamma ~30;
  227b, gamma ~120), test on the other; m in [9, 10, 20, 30, 50, 100].

## Why the plan's H5 is re-specified
With random sampling inside a family the mean coverage is fixed by construction: for random
labels the family coverage of split and offset_full is exactly k/(m + 1), and for one_per_group
labels so is the coverage of new compositions (Dunn et al. 2023). Fd-3m heterogeneity cannot move
these means. What composition structure can change is the spread across draws, i.e. the chance
that a particular set of labels gives a badly calibrated interval (calibration-conditional
coverage; Vovk 2012; Angelopoulos & Bates 2023, Sec. 3.2), and, for clustered labels, the mean as
well. H5 is therefore tested on the spread, against a permutation null with the same group sizes
and no group structure.

## Decision rules (each target separately)
D5.0 Code check (random scheme, split and offset_full): in every family and m,
     |mean coverage - k/(m + 1)| <= 4 Monte Carlo SE (split: exact-law SD / sqrt(draws);
     offset_full: empirical SD / sqrt(draws)). A cell outside this band is re-run once with
     10 x the draws and a new seed; if it is still outside, that is a code error: stop and do
     not interpret Phase 5.
H5   (primary; plan H5 for the recommended protocol) one_per_group, split, m = 20:
     rho = SD(cov_new) / SD(cov_new under the permutation null).
     Supported if rho <= 1.25 in >= 6 of 7 families; rejected if rho > 1.5 in >= 2
     families; otherwise partial. Rejection means the labels needed depend on which compositions
     are labelled, beyond exchangeable theory.
H5c  (clustered labels) clustered, split, m = 20: a family counts as harmed if its rho > 1.25 or
     its mean cov_new < k/(m + 1) - 0.03. Supported (clustering harms) if >= 4 families are
     harmed; otherwise not supported.
H5b  (offset correction) Cells fixed in advance from Phase 2 E1 (|bias|/MAE >= 0.80):
     band gap Fd-3m (LOSO), C2/m (LOSO), R-3c (LOSO), Ia-3 (LOSO), P63/mmc (refit); formation Fd-3m (LOSO), C2/m (LOSO), R-3c (LOSO), Ia-3 (LOSO), Pna21 (LOSO), P63/mmc (refit). Statistic at one_per_group, m = 20: width ratio =
     median width (offset_full) / median width (split), and mean cov_new of offset_full.
     Supported if ratio <= 0.8 and mean coverage >= 0.88 in >= 80% of the cells;
     rejected if ratio > 0.95 in at least half of the cells; otherwise partial.
Practical endpoint (reported for every family, scheme and method; stated expectation, not a
     rule): m* = smallest m with mean coverage >= 0.88, P(coverage < 0.8) <= 0.1 and no
     unbounded interval. Expected under exchangeable theory: m* = 30 (Beta P(coverage < 0.80) =
     0.069 at m = 20 and 0.044 at m = 30), later if the test set is small.
Exploratory, with expectations stated now: (i) clustered calibration is more dispersed than
     one_per_group at equal m; (ii) sub-family transfer for Fd-3m: band gap near nominal (offsets
     +0.625 vs +0.673 eV), formation 227a -> 227b below nominal (offsets -0.099 vs -0.147 eV/atom).
Also reported: interval width in units of the family's SD next to the LOSO CV+ width (Phase 2),
     and MAE before and after subtracting the estimated offset.

## Optional training arm (`train5`, descriptive, no decision rule)
Families Fd-3m (LOSO); P63/mmc (refit without Fd-3m). Per draw, 40% of the composition groups are held out as test; from the
other groups, n_cal = 20 structures are drawn for calibration and m in [10, 30, 100] more for training.
Arms: refit with the m structures added, calibrated on n_cal (train); the fixed Phase 2 model
calibrated on m + n_cal (calibrate); the fixed Phase 2 model calibrated on n_cal (fixed); a refit
without the m structures calibrated on n_cal (baseline; separates the refit from the labels).
10 draws per m, Phase 2 hyperparameters, split and offset_full. Question: are new DFT
labels better spent on the model or on calibration?

Every family, target, scheme, m and method is reported whichever way it falls.
