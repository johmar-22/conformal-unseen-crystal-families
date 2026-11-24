"""cpaudit_run_v2.py - one-line launcher for the Phase 1 follow-up, Phase 2 and Phase 3 bundle.

Put this file and cpaudit_v2.zip in  MyDrive/NOMAD_Uncer_Quan/code_bundles/  then, in Colab:

    from google.colab import drive; drive.mount("/content/drive")
    %run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v2.py" setup

Steps (run in this order)
    setup              self-test + v2 de-duplication / integrity + Phase 2 pre-registration
    post1              exploratory checks on the saved Phase 1 results (minutes, no refitting)
    overlap            Phase 3 overlap diagnostics (minutes; must precede phase2)
    phase2             nested LOSO + leave-two-out for both targets (long; resumes after a disconnect)
    phase2 bandgap     one target only (phase2 formation for the other)
    report2            tables, pre-registered decisions, figure
    selftest | integrity | prereg   single steps

The zip is unpacked to /content/cpaudit_v2 (local disk, never into Drive) and re-unpacked
automatically whenever you upload a new zip. Missing packages are installed once per session.
All outputs go to MyDrive/NOMAD_Uncer_Quan/phase02_v1/. Phase 1 files are only read.
"""
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(globals().get("__file__", "cpaudit_run_v2.py")))
ZIP = os.environ.get("CPAUDIT_ZIP", os.path.join(HERE, "cpaudit_v2.zip"))
DEST = os.environ.get("CPAUDIT_DEST", "/content")
CODE = os.path.join(DEST, "cpaudit_v2")
STEPS = {"setup": [["selftest"], ["integrity"], ["prereg"]], "selftest": [["selftest"]],
         "integrity": [["integrity"]], "prereg": [["prereg"]], "post1": [["post1"]],
         "overlap": [["overlap"]], "phase2": [["phase2"]], "report2": [["report2"]]}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _unpack():
    if not os.path.exists(ZIP):
        raise FileNotFoundError(f"cpaudit_v2.zip not found next to the launcher: {ZIP}")
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
    need = [p for m, p in (("optuna", "optuna>=3.6"), ("yaml", "pyyaml"), ("tabulate", "tabulate"),
                           ("pymatgen", "pymatgen")) if importlib.util.find_spec(m) is None]
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
        if step == "phase2" and extra and not extra[0].startswith("-"):
            cmd = cmd + ["--targets", extra[0]] + extra[1:]
        elif extra:
            cmd = cmd + extra
        _run(cmd)
    print(f"\n[launcher] '{step}' finished.")


main(sys.argv[1:])
