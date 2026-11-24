"""Difficulty estimate for normalized conformal scores.

sigma(x) = mean Euclidean distance from x to its k nearest proper-training points, in the
space of the first `pca_components` principal components of the standardized SOAP vector,
plus a small additive beta (crepes convention; Bostrom 2022). Following Papadopoulos et al.
(2011) and Lei et al. (2018), sigma is fitted on the proper-training set only, so it is a fixed
function for calibration and test points and the calibration/test scores stay exchangeable.
PCA is used because nearest-neighbour distances in ~24k raw dimensions concentrate and lose
contrast; the explained-variance fraction is reported so the choice can be judged.
"""
from __future__ import annotations

import numpy as np
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


class KNNDifficulty:
    def __init__(self, k: int = 25, n_components: int = 50, beta: float = 0.01, seed: int = 0):
        self.k, self.n_components, self.beta, self.seed = k, n_components, beta, seed

    def fit(self, X_proper: np.ndarray):
        Xs = StandardScaler().fit(X_proper)
        self.scaler_ = Xs
        Z = Xs.transform(X_proper)
        self.pca_ = PCA(n_components=self.n_components, svd_solver="randomized",
                        random_state=self.seed).fit(Z)
        self.explained_ = float(self.pca_.explained_variance_ratio_.sum())
        self.nn_ = NearestNeighbors(n_neighbors=self.k).fit(self.pca_.transform(Z))
        return self

    def sigma(self, X: np.ndarray) -> np.ndarray:
        d, _ = self.nn_.kneighbors(self.pca_.transform(self.scaler_.transform(X)))
        return d.mean(axis=1) + self.beta
