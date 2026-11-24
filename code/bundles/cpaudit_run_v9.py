"""cpaudit_run_v9.py - one-line launcher for cpaudit v9: the final paper figures (Fig. 1-8 and S1-S7),
drawn only from result files already on the Drive. Nothing is fitted or recomputed.

Put this file and cpaudit_v9.zip in  MyDrive/NOMAD_Uncer_Quan/code_bundles/  then, in Colab:

    from google.colab import drive; drive.mount("/content/drive", force_remount=True)
    %run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v9.py" figures9 --only fig1

Steps
    figures9 --only fig1 [fig5 ...]  chosen figures -> paper_figures_v5_drafts/<UTC time>_<names>/ (new folder
                                     every run; for checking one figure at a time)
    figures9                         the full set -> paper_figures_v5/ (refused if that folder exists and is
                                     not empty: nothing is ever overwritten), with FIGURES9_MANIFEST.json
                                     (SHA-256 of every input file, overlap-check result)
    setup                            all self-tests (about 5 min; optional for the figures)

Existing Drive files are only read. The earlier steps (Phases 1-8) stay in the menu but are not needed.
The launcher prints the SHA-256 of the zip it runs, so you can check it is the new one.
"""
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(globals().get("__file__", "cpaudit_run_v9.py")))
ZIP = os.environ.get("CPAUDIT_ZIP", os.path.join(HERE, "cpaudit_v9.zip"))
DEST = os.environ.get("CPAUDIT_DEST", "/content")
CODE = os.path.join(DEST, "cpaudit_v9")
STEPS = {"setup": [["selftest"]], "selftest": [["selftest"]], "freeze6": [["freeze6"]],
         "features6": [["features6"]], "splits6": [["splits6"]], "matfold6": [["matfold6"]],
         "report6": [["report6"]], "setup7": [["setup7"]], "tune7": [["tune7"]], "fit7": [["fit7"]],
         "report7": [["report7"]],
         "setup8": [["setup8"]], "scores8": [["scores8"]], "fewshot8": [["fewshot8"]], "hjplus8": [["hjplus8"]],
         "report8": [["report8"]], "figures8": [["figures8"]],
         "posthoc8": [["posthoc8"]], "figures81": [["figures81"]], "figures9": [["figures9"]], "posthoc5": [["posthoc5"]], "audit6": [["audit6"]], "figures": [["figures"]]}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _unpack():
    if not os.path.exists(ZIP):
        raise FileNotFoundError(f"cpaudit_v9.zip not found next to the launcher: {ZIP}")
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


LIGHT = {"figures9"}                    # drawing only: pandas, numpy, scipy, scikit-learn, matplotlib (all in Colab)


def _install(step):
    if step in LIGHT:
        if importlib.util.find_spec("yaml") is None:
            subprocess.run([sys.executable, "-m", "pip", "install", "-q", "pyyaml"], check=True)
        return
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
    _install(step)
    for cmd in STEPS[step]:
        _run(cmd + extra if extra else cmd)
    print(f"\n[launcher] '{step}' finished.")


main(sys.argv[1:])
