"""cpaudit_run_v3.py - one-line launcher for Phase 4 (detecting failure before it matters).

Put this file and cpaudit_v3.zip in  MyDrive/NOMAD_Uncer_Quan/code_bundles/  then, in Colab:

    from google.colab import drive; drive.mount("/content/drive")
    %run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v3.py" setup

Steps (run in this order)
    setup      self-test + link the Phase 2 inputs + Phase 4 pre-registration (a few minutes)
    phase4     detection scores for all 11 configurations (about 20-40 min; resumes per configuration)
    report4    metrics with cluster-bootstrap CIs, abstention curves, flags, decisions, figure
    figures    paper figures (Fig. 1-6, S1-S3) from the saved Phase 1-4 outputs (1-2 min; read only)
               add figure prefixes to redraw only some, e.g.  ... cpaudit_run_v3.py" figures --only fig3 fig4
    selftest | prereg   single steps

No model is refitted. Phase 1 and 2 files are only read. Outputs go to
MyDrive/NOMAD_Uncer_Quan/phase04_v1/ and (figures) MyDrive/NOMAD_Uncer_Quan/paper_figures_v1/.
The zip is unpacked to /content/cpaudit_v3 (local disk)
and re-unpacked automatically whenever you upload a new zip.
"""
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(globals().get("__file__", "cpaudit_run_v3.py")))
ZIP = os.environ.get("CPAUDIT_ZIP", os.path.join(HERE, "cpaudit_v3.zip"))
DEST = os.environ.get("CPAUDIT_DEST", "/content")
CODE = os.path.join(DEST, "cpaudit_v3")
STEPS = {"setup": [["selftest"], ["prereg"]], "selftest": [["selftest"]], "prereg": [["prereg"]],
         "phase4": [["phase4"]], "report4": [["report4"]], "figures": [["figures"]]}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _unpack():
    if not os.path.exists(ZIP):
        raise FileNotFoundError(f"cpaudit_v3.zip not found next to the launcher: {ZIP}")
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
