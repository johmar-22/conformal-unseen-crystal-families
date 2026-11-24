#!/usr/bin/env python3
"""cpaudit v3 command line (Phase 4: detecting failure before it matters).

  python run.py selftest      conformal, structure-matching and detection code vs theory / independent code
  python run.py prereg        link the Phase 2 inputs and write the immutable Phase 4 pre-registration
  python run.py phase4 [--allow-deviation]   detection scores for every held-out structure
  python run.py report4       metrics, selective prediction, lattice level, flags, decisions, figure
  python run.py figures [--only fig1 fig4]   paper figures from the saved Phase 1-4 outputs (read only)
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cpaudit.config import load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "prereg", "phase4", "report4", "figures"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--allow-deviation", action="store_true")
    ap.add_argument("--only", nargs="*", default=None, help="figure name prefixes, e.g. fig1 fig4 figS")
    a = ap.parse_args()
    if a.command == "selftest":
        from cpaudit.selftest import run_selftest
        sys.exit(run_selftest())
    cfg = load_config(a.config)
    if a.command == "prereg":
        from cpaudit.setup4 import preregister
        preregister(cfg)
    elif a.command == "phase4":
        from cpaudit.phase4 import run_phase4
        run_phase4(cfg, allow_deviation=a.allow_deviation)
    elif a.command == "report4":
        from cpaudit.report4 import run_report4
        run_report4(cfg)
    elif a.command == "figures":
        from cpaudit.paperfigs import run_figures
        m = run_figures(cfg, a.only)
        print(f"[figures] {len(m['figures'])} drawn, {len(m['skipped'])} skipped; manifest FIGURES_MANIFEST.json")


if __name__ == "__main__":
    main()
