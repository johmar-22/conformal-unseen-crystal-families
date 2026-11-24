"""cpaudit_run_v5.py - one-line launcher for cpaudit v5 (Phase 5 post hoc tables, Phase 6 Task 6.1
JARVIS-DFT data audit, all paper figures redrawn from Drive with an automatic overlap check).

Put this file and cpaudit_v5.zip in  MyDrive/NOMAD_Uncer_Quan/code_bundles/  then, in Colab:

    from google.colab import drive; drive.mount("/content/drive")
    %run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v5.py" setup

Steps (run in this order; nothing refits a model)
    setup      self-test (about 3-5 min; installs pymatgen once per session)
    posthoc5   Phase 5 post hoc tables, NOT pre-registered: no-model baseline from the same labels and
               residual ICC by composition (about 1 min). Re-creates the Phase 5 label sets with the
               registered seeds and stops if any draw differs from the saved Phase 5 draws.
    audit6     Phase 6 Task 6.1: downloads the pinned JARVIS-DFT release(s) to Drive once (about 90 MB),
               counts structures per space group, applies the scope rule fixed in the config
               (about 2-4 min)
    figures    all paper figures (Fig. 1-7, S1-S6) from the Drive files, each checked for overlapping
               text, legends and labels (about 2 min)
    all        posthoc5, audit6, figures in one go

Existing Drive files are only read. New outputs: NOMAD_Uncer_Quan/phase05_v1/posthoc/,
NOMAD_Uncer_Quan/phase06_v1/ and NOMAD_Uncer_Quan/paper_figures_v2/.
The zip is unpacked to /content/cpaudit_v5 (local disk) and re-unpacked automatically whenever
you upload a new zip.
"""
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(globals().get("__file__", "cpaudit_run_v5.py")))
ZIP = os.environ.get("CPAUDIT_ZIP", os.path.join(HERE, "cpaudit_v5.zip"))
DEST = os.environ.get("CPAUDIT_DEST", "/content")
CODE = os.path.join(DEST, "cpaudit_v5")
STEPS = {"setup": [["selftest"]], "selftest": [["selftest"]], "posthoc5": [["posthoc5"]],
         "audit6": [["audit6"]], "figures": [["figures"]],
         "all": [["posthoc5"], ["audit6"], ["figures"]]}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _unpack():
    if not os.path.exists(ZIP):
        raise FileNotFoundError(f"cpaudit_v5.zip not found next to the launcher: {ZIP}")
    stamp, digest = os.path.join(CODE, ".zip_sha256"), _sha(ZIP)
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
    # crepes and MAPIE only cross-check the conformal code; --no-deps keeps Colab's numpy/sklearn untouched
    nodeps = [p for m, p in (("crepes", "crepes"), ("mapie", "mapie==0.8.6")) if importlib.util.find_spec(m) is None]
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
