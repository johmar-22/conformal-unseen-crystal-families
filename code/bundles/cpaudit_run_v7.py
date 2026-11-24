"""cpaudit_run_v7.py - one-line launcher for cpaudit v7 (Phase 6 Task 6.4 on JARVIS-DFT: LightGBM
models and conformal intervals, pre-registered). Earlier steps (freeze6 ... report6) are included.

Put this file and cpaudit_v7.zip in  MyDrive/NOMAD_Uncer_Quan/code_bundles/  then, in Colab:

    from google.colab import drive; drive.mount("/content/drive", force_remount=True)
    %run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v7.py" setup

Steps (in this order; every long step is resumable: if Colab disconnects, run the same cell again)
    setup      self-tests (about 5 min; installs pymatgen, matminer, MatFold, optuna once per session)
    setup7     builds the design, hashes every input, times one fit, writes the PRE-REGISTRATION
               (phase07_v1/phase0/PREREGISTRATION_phase6.md). No model is fitted before this.
    tune7      cross-half nested tuning (15 trials each; resumable per trial)
               options: --target formation bandgap tbmbj
    fit7       all model units, in the pre-registered order; options: --target ..., --block ...
               blocks: primary regimeA randinner polyex h3 rep2 locso
    report7    Task 6.4 tables and the decisions G5, H1-J, G3-J, H3-J

Existing Drive files are only read. New outputs: NOMAD_Uncer_Quan/phase07_v1/.
The launcher prints the SHA-256 of the zip it runs, so you can check it is the new one.
"""
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(globals().get("__file__", "cpaudit_run_v7.py")))
ZIP = os.environ.get("CPAUDIT_ZIP", os.path.join(HERE, "cpaudit_v7.zip"))
DEST = os.environ.get("CPAUDIT_DEST", "/content")
CODE = os.path.join(DEST, "cpaudit_v7")
STEPS = {"setup": [["selftest"]], "selftest": [["selftest"]], "freeze6": [["freeze6"]],
         "features6": [["features6"]], "splits6": [["splits6"]], "matfold6": [["matfold6"]],
         "report6": [["report6"]], "setup7": [["setup7"]], "tune7": [["tune7"]], "fit7": [["fit7"]],
         "report7": [["report7"]], "posthoc5": [["posthoc5"]], "audit6": [["audit6"]], "figures": [["figures"]]}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _unpack():
    if not os.path.exists(ZIP):
        raise FileNotFoundError(f"cpaudit_v7.zip not found next to the launcher: {ZIP}")
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
    need = [p for m, p in (("yaml", "pyyaml"), ("tabulate", "tabulate"), ("pymatgen", "pymatgen"),
                           ("lightgbm", "lightgbm"), ("optuna", "optuna"))
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
