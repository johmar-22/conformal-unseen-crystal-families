"""cpaudit_run_v4.py - one-line launcher for Phase 5 (few-shot, family-specific recalibration).

Put this file and cpaudit_v4.zip in  MyDrive/NOMAD_Uncer_Quan/code_bundles/  then, in Colab:

    from google.colab import drive; drive.mount("/content/drive")
    %run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v4.py" setup

Steps (run in this order)
    setup      self-test + link the Phase 2 inputs + Phase 5 pre-registration (about 5 min)
    phase5     few-shot recalibration draws for 7 families x 2 targets (about 3-6 min; no refitting;
               resumes per family)
    report5    summaries, decisions D5.0 / H5 / H5c / H5b, practical endpoint (about 1-3 min)
    figures    paper figures Fig. 1-7 and S1-S5 from the saved outputs (1-2 min; read only)
    train5     OPTIONAL training arm: refits LightGBM (about 124 fits, roughly 1-3 h on 8 CPUs;
               resumes per draw). Run report5 and figures again afterwards.
    selftest | prereg   single steps

Phase 1, 2 and 4 files are only read. Outputs go to MyDrive/NOMAD_Uncer_Quan/phase05_v1/ and
(figures) MyDrive/NOMAD_Uncer_Quan/paper_figures_v1/.
The zip is unpacked to /content/cpaudit_v4 (local disk)
and re-unpacked automatically whenever you upload a new zip.
"""
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(globals().get("__file__", "cpaudit_run_v4.py")))
ZIP = os.environ.get("CPAUDIT_ZIP", os.path.join(HERE, "cpaudit_v4.zip"))
DEST = os.environ.get("CPAUDIT_DEST", "/content")
CODE = os.path.join(DEST, "cpaudit_v4")
STEPS = {"setup": [["selftest"], ["prereg"]], "selftest": [["selftest"]], "prereg": [["prereg"]],
         "phase5": [["phase5"]], "train5": [["train5"]], "report5": [["report5"]], "figures": [["figures"]]}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _unpack():
    if not os.path.exists(ZIP):
        raise FileNotFoundError(f"cpaudit_v4.zip not found next to the launcher: {ZIP}")
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
    need = [p for m, p in (("yaml", "pyyaml"), ("tabulate", "tabulate")) if importlib.util.find_spec(m) is None]
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
        if extra:
            cmd = cmd + extra
        _run(cmd)
    print(f"\n[launcher] '{step}' finished.")


main(sys.argv[1:])
