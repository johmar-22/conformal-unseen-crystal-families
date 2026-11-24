"""SOAP descriptors for NOMAD 2018 (the recipe behind soap_features_corrected_v2.npz).

Per structure: DScribe SOAP for every atom (species Al, Ga, In, O; periodic; r_cut 6.0 A,
n_max 9, l_max 8; dense output), then the mean, standard deviation, minimum and maximum over
atoms, concatenated: 4 x 5,994 = 23,976 features. Rows follow train.csv.

Usage:  python make_soap_features.py --train-csv train.csv --geometry-dir train/ --out soap.npz
(geometry files: train/<id>/geometry.xyz from the NOMAD 2018 Kaggle data, FHI-aims format)
"""
import argparse, os
import numpy as np
import pandas as pd
from ase.io import read as ase_read
from dscribe.descriptors import SOAP

ap = argparse.ArgumentParser()
ap.add_argument("--train-csv", required=True)
ap.add_argument("--geometry-dir", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

ids = pd.read_csv(a.train_csv)["id"].values
atoms = [ase_read(os.path.join(a.geometry_dir, str(i), "geometry.xyz"), format="aims") for i in ids]
species = sorted({s for at in atoms for s in at.get_chemical_symbols()})
soap = SOAP(species=species, periodic=True, r_cut=6.0, n_max=9, l_max=8, sparse=False)

def aggregate(x):
    return np.concatenate([x.mean(0), x.std(0), x.min(0), x.max(0)])

X = np.array([aggregate(soap.create(at)) for at in atoms])
np.savez_compressed(a.out, features=X)
print(X.shape)
