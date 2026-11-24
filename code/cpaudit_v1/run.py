#!/usr/bin/env python3
"""cpaudit command line.

  python run.py selftest                 conformal code vs theory and vs crepes/MAPIE (no data needed)
  python run.py phase0 [--no-hash-large] freeze manifest, environment, integrity checks
  python run.py prereg                   write the immutable Phase 1 pre-registration
  python run.py phase1 [--repeats 0-24] [--targets bandgap,formation] [--allow-deviation]
  python run.py report                   tables, pre-registered endpoint, figure
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cpaudit.config import load_config  # noqa: E402


def _range(s):
    if s is None:
        return None
    out = []
    for part in s.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "phase0", "prereg", "phase1", "report"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--repeats", default=None, help="e.g. 0-24 or 0,3,5-7")
    ap.add_argument("--targets", default="bandgap,formation")
    ap.add_argument("--allow-deviation", action="store_true")
    ap.add_argument("--no-hash-large", action="store_true", help="skip SHA-256 of files > 50 MB in phase0")
    a = ap.parse_args()

    if a.command == "selftest":
        from cpaudit.selftest import run_selftest
        sys.exit(run_selftest())
    cfg = load_config(a.config)
    if a.command == "phase0":
        from cpaudit.phase0 import run_phase0
        run_phase0(cfg, hash_large=not a.no_hash_large)
    elif a.command == "prereg":
        from cpaudit.phase0 import preregister
        preregister(cfg)
    elif a.command == "phase1":
        from cpaudit.phase1 import run_phase1
        run_phase1(cfg, repeats=_range(a.repeats), targets=a.targets.split(","),
                   allow_deviation=a.allow_deviation)
    elif a.command == "report":
        from cpaudit.report import run_report
        run_report(cfg)


if __name__ == "__main__":
    main()
