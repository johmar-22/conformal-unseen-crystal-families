"""Redraw all paper figures from results/ with the final bundle (code/cpaudit_v9).

  python tools/render_figures.py --out /tmp/figures_check     (the folder must not exist or be empty)
Needs results/phase02_v1/inputs/train.csv (python tools/fetch_inputs.py).
"""
import argparse, os, subprocess, sys, tempfile
import yaml

here = os.path.dirname(os.path.abspath(__file__))
repo = os.path.abspath(os.path.join(here, ".."))
ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
a = ap.parse_args()
code = os.path.join(repo, "code", "cpaudit_v9")
cfg = yaml.safe_load(open(os.path.join(code, "config", "phase06.yaml")))
cfg["paths"]["drive_root"] = os.path.join(repo, "results")
cfg.setdefault("figures9", {})["figures_dir"] = os.path.abspath(a.out)
cfg["figures9"]["drafts_dir"] = os.path.join(tempfile.mkdtemp(), "drafts")
tmp = os.path.join(tempfile.mkdtemp(), "figures_only.yaml")
yaml.safe_dump(cfg, open(tmp, "w"), sort_keys=False)
env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
sys.exit(subprocess.call([sys.executable, "run.py", "figures9", "--config", tmp], cwd=code, env=env))
