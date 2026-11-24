"""Download the third-party inputs into results/ and check their SHA-256.

  python tools/fetch_inputs.py              NOMAD train.csv (needed by most steps and the figures)
  python tools/fetch_inputs.py --jarvis     also the JARVIS-DFT release zip (48 MB)
  python tools/fetch_inputs.py --assets TAG also the two release assets of this repository
"""
import argparse, hashlib, os, shutil, urllib.request

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
NOMAD = ("https://raw.githubusercontent.com/csutton7/nomad_2018_kaggle_dataset/master/train.csv",
         "2be6168b578c4be91b7380e7404a2b58912790d8e03c51b28f4f279f2d8bb872")
JARVIS = ("https://ndownloader.figshare.com/files/64391379", "phase06_v1/inputs/jdft_3d-9-24-2025.json.zip")
REPO = "https://github.com/johmar-22/conformal-unseen-crystal-families/releases/download"

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()

def get(url, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with urllib.request.urlopen(url) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f)
    print("downloaded", dest, sha(dest)[:12])

ap = argparse.ArgumentParser()
ap.add_argument("--jarvis", action="store_true")
ap.add_argument("--assets", metavar="TAG")
a = ap.parse_args()

first = os.path.join(ROOT, "phase01_v1", "inputs", "train.csv")
get(NOMAD[0], first)
if sha(first) != NOMAD[1]:
    raise SystemExit("train.csv SHA-256 differs from the file used in the paper")
for d in ("phase02_v1", "phase04_v1", "phase05_v1"):
    os.makedirs(os.path.join(ROOT, d, "inputs"), exist_ok=True)
    shutil.copy2(first, os.path.join(ROOT, d, "inputs", "train.csv"))
if a.jarvis:
    get(JARVIS[0], os.path.join(ROOT, JARVIS[1]))
if a.assets:
    get(f"{REPO}/{a.assets}/soap_features_corrected_v2.npz",
        os.path.join(ROOT, "soap_features_corrected_v2.npz"))
    get(f"{REPO}/{a.assets}/features_v1.npz",
        os.path.join(ROOT, "phase06_v1", "features", "features_v1.npz"))
print("Compare the hashes with docs/DATA.md.")
