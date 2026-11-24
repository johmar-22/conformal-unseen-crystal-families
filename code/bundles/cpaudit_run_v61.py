"""cpaudit_run_v61.py - one-line launcher for cpaudit v6.1 (Phase 6 Tasks 6.2-6.3 on JARVIS-DFT:
data freeze, features, MatFold leave-one-group-out splits). No target is used and no model is fitted.

v6.1 = v6 + the MatFold review: MatFold's own labels (SpacegroupAnalyzer, symprec 0.1) are now the
primary split labels, the JARVIS labels a sensitivity family, and MatFold itself is run on all
structures to prove that its splits equal ours row by row (new step matfold6).

Put this file and cpaudit_v61.zip in  MyDrive/NOMAD_Uncer_Quan/code_bundles/  then, in Colab:

    from google.colab import drive; drive.mount("/content/drive", force_remount=True)
    %run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v61.py" setup

Steps (run in this order; freeze6 and features6 are unchanged from v6, skip them if v6 already ran)
    setup      self-test (about 3-5 min; installs pymatgen, matminer and MatFold once per session)
    freeze6    rebuild the JARVIS analysis table from the release on Drive; stops unless it equals
               the Task 6.1 audit exactly (about 1 min)
    features6  Magpie 145 + density 3 + CrystalNN 122 features for 81,884 structures; uses all CPUs,
               saves a checkpoint every 2,000 structures, resumes after a disconnect
               (about 15-25 min on 8 CPUs, about 1 h on 2)
    splits6    LOSGO (59) and LOCSO (7) holdouts under MatFold labels (primary) and JARVIS labels
               (sensitivity), with the polymorph-exclusion variant (about 5-8 min; reuses the
               space-group labels if v6 already computed them)
    matfold6   MatFold itself on all structures; its labels and splits must equal ours (about 10 min)
    report6    one summary (seconds)
    all        freeze6, features6, splits6, matfold6, report6 in one go

Existing Drive files are only read. New outputs: NOMAD_Uncer_Quan/phase06_v1/{data,features,splits,
report,splits_matfold}/ (v6 files such as holdouts_v1.npz are left untouched).
The launcher prints the SHA-256 of the zip it runs, so you can check it is the new one.
"""
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(globals().get("__file__", "cpaudit_run_v61.py")))
ZIP = os.environ.get("CPAUDIT_ZIP", os.path.join(HERE, "cpaudit_v61.zip"))
DEST = os.environ.get("CPAUDIT_DEST", "/content")
CODE = os.path.join(DEST, "cpaudit_v61")
STEPS = {"setup": [["selftest"]], "selftest": [["selftest"]], "freeze6": [["freeze6"]],
         "features6": [["features6"]], "splits6": [["splits6"]], "matfold6": [["matfold6"]],
         "report6": [["report6"]], "all": [["freeze6"], ["features6"], ["splits6"], ["matfold6"], ["report6"]],
         "posthoc5": [["posthoc5"]], "audit6": [["audit6"]], "figures": [["figures"]]}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _unpack():
    if not os.path.exists(ZIP):
        raise FileNotFoundError(f"cpaudit_v61.zip not found next to the launcher: {ZIP}")
    stamp, digest = os.path.join(CODE, ".zip_sha256"), _sha(ZIP)
    print(f"[launcher] {os.path.basename(ZIP)} sha256 {digest[:12]}", flush=True)
    if os.path.exists(stamp) and open(stamp).read().strip() == digest:
        return
    shutil.rmtree(CODE, ignore_errors=True)
    with zipfile.ZipFile(ZIP) as z:
        z.extractall(DEST)
    with open(stamp, "w") as f:
        f.write(digest)
    print(f"[launcher] unpacked {os.path.basename(ZIP)} -> {CODE}")


def _install():
    need = [p for m, p in (("yaml", "pyyaml"), ("tabulate", "tabulate"), ("pymatgen", "pymatgen"))
            if importlib.util.find_spec(m) is None]
    if need:
        print(f"[launcher] installing {need} (once per session)", flush=True)
        subprocess.run([sys.executable, "-m", "pip", "install", "-q"] + need, check=True)
    # --no-deps: matminer pins pandas<3 and scikit-learn~=1.3 and would change Colab's packages; its
    # featurizers need only numpy, pandas, scikit-learn, sympy, monty and pymatgen (all present).
    # crepes and MAPIE only cross-check the conformal code.
    nodeps = [p for m, p in (("matminer", "matminer==0.10.1"), ("MatFold", "MatFold==1.3.3"), ("crepes", "crepes"),
                             ("mapie", "mapie==0.8.6")) if importlib.util.find_spec(m) is None]
    if nodeps:
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "--no-deps"] + nodeps, check=True)


def _run(args):
    print(f"\n[launcher] python run.py {' '.join(args)}", flush=True)
    p = subprocess.Popen([sys.executable, "-u", os.path.join(CODE, "run.py")] + args, cwd=CODE,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in p.stdout:
        print(line, end="", flush=True)
    if p.wait() != 0:
        raise SystemExit(f"[launcher] step '{' '.join(args)}' failed (exit {p.returncode}); see the output above.")


def main(argv):
    if not argv or argv[0] not in STEPS:
        print(__doc__)
        return
    step, extra = argv[0], argv[1:]
    _unpack()
    _install()
    for cmd in STEPS[step]:
        _run(cmd + extra if extra else cmd)
    print(f"\n[launcher] '{step}' finished.")


main(sys.argv[1:])
