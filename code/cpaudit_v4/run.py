#!/usr/bin/env python3
"""cpaudit v4 command line (Phase 5: few-shot, family-specific recalibration).

  python run.py selftest      conformal, few-shot and detection code vs theory / independent code
  python run.py prereg        link the Phase 2 inputs and write the immutable Phase 5 pre-registration
  python run.py phase5 [--allow-deviation]   few-shot recalibration draws (no refitting)
  python run.py train5 [--allow-deviation]   optional training arm (refits LightGBM; slow)
  python run.py report5       summaries, decisions D5.0 / H5 / H5b, practical endpoint
  python run.py figures [--only fig7 figS4]  paper figures from the saved Phase 1-5 outputs
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cpaudit.config import load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "prereg", "phase5", "train5", "report5", "figures"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--allow-deviation", action="store_true")
    ap.add_argument("--only", nargs="*", default=None, help="figure name prefixes, e.g. fig7 figS")
    a = ap.parse_args()
    if a.command == "selftest":
        from cpaudit.selftest import run_selftest
        sys.exit(run_selftest())
    cfg = load_config(a.config)
    if a.command == "prereg":
        from cpaudit.setup5 import preregister
        preregister(cfg)
    elif a.command == "phase5":
        from cpaudit.phase5 import run_phase5
        run_phase5(cfg, allow_deviation=a.allow_deviation)
    elif a.command == "train5":
        from cpaudit.train5 import run_train5
        run_train5(cfg, allow_deviation=a.allow_deviation)
    elif a.command == "report5":
        from cpaudit.report5 import run_report5
        run_report5(cfg)
    elif a.command == "figures":
        from cpaudit.paperfigs import run_figures
        m = run_figures(cfg, a.only)
        print(f"[figures] {len(m['figures'])} drawn, {len(m['skipped'])} skipped; manifest FIGURES_MANIFEST.json")


if __name__ == "__main__":
    main()
