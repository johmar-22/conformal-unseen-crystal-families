#!/usr/bin/env python3
"""cpaudit v2 command line (Phase 1 follow-up, Phase 2, Phase 3).

  python run.py selftest      conformal and structure-matching code vs theory / independent code
  python run.py integrity     environment + v2 de-duplication + data integrity report
  python run.py prereg        immutable Phase 2/3 pre-registration
  python run.py post1         exploratory checks on the saved Phase 1 outputs (no refitting)
  python run.py overlap       Phase 3 overlap diagnostics (before any LOSO result)
  python run.py phase2 [--targets bandgap,formation] [--allow-deviation]
  python run.py report2       Phase 2/3 tables, pre-registered decisions, figure
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cpaudit.config import load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "integrity", "prereg", "post1", "overlap", "phase2", "report2"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--targets", default="bandgap,formation")
    ap.add_argument("--allow-deviation", action="store_true")
    a = ap.parse_args()

    if a.command == "selftest":
        from cpaudit.selftest import run_selftest
        sys.exit(run_selftest())
    cfg = load_config(a.config)
    if a.command == "integrity":
        from cpaudit.setup2 import run_integrity
        run_integrity(cfg)
    elif a.command == "prereg":
        from cpaudit.setup2 import preregister
        preregister(cfg)
    elif a.command == "post1":
        from cpaudit.post1 import run_post1
        run_post1(cfg)
    elif a.command == "overlap":
        from cpaudit.overlap import run_overlap
        run_overlap(cfg, allow_deviation=a.allow_deviation)
    elif a.command == "phase2":
        from cpaudit.phase2 import run_phase2
        run_phase2(cfg, targets=[t for t in a.targets.split(",") if t], allow_deviation=a.allow_deviation)
    elif a.command == "report2":
        from cpaudit.report2 import run_report2
        run_report2(cfg)


if __name__ == "__main__":
    main()
