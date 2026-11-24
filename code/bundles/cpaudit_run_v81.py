"""cpaudit_run_v81.py - one-line launcher for cpaudit v8.1: POST HOC analyses and figure fixes after the
Phase 6 Task 6.5 results (NOT pre-registered). It reads the v7 and v8 outputs and writes new folders only.

Put this file and cpaudit_v81.zip in  MyDrive/NOMAD_Uncer_Quan/code_bundles/  then, in Colab:

    from google.colab import drive; drive.mount("/content/drive", force_remount=True)
    %run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v81.py" setup

v8.1 steps
    setup      self-tests (about 5 min), including the new v8.1 checks
    posthoc8   few-shot intervals vs a no-model interval from the same labels (the fewshot8 label sets are
               re-created and checked draw by draw); inner-model SD at group level and abstention inside the
               worst groups; ABO3 oxides; degenerate group-level intervals -> phase08_v1/posthoc/ (15-25 min)
    figures81  Fig. 8-10 with the v8.1 fixes -> paper_figures_v4/ (about 1 min)

The v8 steps below stay in the menu, but setup8 / scores8 / fewshot8 / hjplus8 / report8 now stop with
"code changed" (the v8.1 code differs from the code fixed in the v8 addendum). That is intended:
the registered v8 results are final and are not recomputed.

Earlier v8 steps

Steps (in this order; every long step is resumable: if Colab disconnects, run the same cell again)
    setup      self-tests (about 5 min; installs packages once per session)
    setup8     checks that v7 is complete and unchanged (rules, inputs, units) and writes the Task 6.5
               pre-registration ADDENDUM (phase08_v1/phase0/PREREGISTRATION_phase6_v8.md) before any
               v8 result exists
    scores8    distance and failure-detection scores per outer fold (resumable per fold)
    fewshot8   H5-J few-shot recalibration of the 3 worst-covered groups per target (minutes)
    hjplus8    hierarchical jackknife+ (Lee et al. eq. 20); resumable per model; --target bandgap tbmbj formation
    report8    Task 6.5 tables and decisions (H4-J, G4-J, H5-J; descriptive parts)
    figures8   Fig. 8-10 (JARVIS) into paper_figures_v3/

Existing Drive files are only read. New outputs: phase08_v1/posthoc/ and paper_figures_v4/.
The launcher prints the SHA-256 of the zip it runs, so you can check it is the new one.
"""
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(globals().get("__file__", "cpaudit_run_v81.py")))
ZIP = os.environ.get("CPAUDIT_ZIP", os.path.join(HERE, "cpaudit_v81.zip"))
DEST = os.environ.get("CPAUDIT_DEST", "/content")
CODE = os.path.join(DEST, "cpaudit_v81")
STEPS = {"setup": [["selftest"]], "selftest": [["selftest"]], "freeze6": [["freeze6"]],
         "features6": [["features6"]], "splits6": [["splits6"]], "matfold6": [["matfold6"]],
         "report6": [["report6"]], "setup7": [["setup7"]], "tune7": [["tune7"]], "fit7": [["fit7"]],
         "report7": [["report7"]],
         "setup8": [["setup8"]], "scores8": [["scores8"]], "fewshot8": [["fewshot8"]], "hjplus8": [["hjplus8"]],
         "report8": [["report8"]], "figures8": [["figures8"]],
         "posthoc8": [["posthoc8"]], "figures81": [["figures81"]], "posthoc5": [["posthoc5"]], "audit6": [["audit6"]], "figures": [["figures"]]}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _unpack():
    if not os.path.exists(ZIP):
        raise FileNotFoundError(f"cpaudit_v81.zip not found next to the launcher: {ZIP}")
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
