#!/usr/bin/env python3
"""cpaudit v5 command line (Phase 5 post hoc tables, Phase 6 Task 6.1 JARVIS audit, paper figures).

  python run.py selftest            conformal, few-shot, post hoc, audit and figure-check code vs theory
  python run.py posthoc5            Phase 5 post hoc tables (NOT pre-registered): no-model baseline, ICC
  python run.py audit6              Phase 6 Task 6.1: JARVIS-DFT data audit (counts only) and scope decision
  python run.py figures [--only fig7 figS]   all paper figures from the Drive files, with the overlap check

Nothing here refits a model. Inputs on Drive are only read; outputs go to phase05_v1/posthoc/,
phase06_v1/ and paper_figures_v2/.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cpaudit.config import load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "posthoc5", "audit6", "figures"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--only", nargs="*", default=None, help="figure name prefixes, e.g. fig7 figS")
    a = ap.parse_args()
    if a.command == "selftest":
        from cpaudit.selftest import run_selftest
        sys.exit(run_selftest())
    cfg = load_config(a.config)
    if a.command == "posthoc5":
        from cpaudit.posthoc5 import run_posthoc5
        run_posthoc5(cfg)
    elif a.command == "audit6":
        from cpaudit.audit6 import run_audit6
        run_audit6(cfg)
    elif a.command == "figures":
        from cpaudit.paperfigs import run_figures
        m = run_figures(cfg, a.only)
        n = m["overlap_check"]["n_issues"]
        print(f"[figures] {len(m['figures'])} drawn, {len(m['skipped'])} skipped; overlap issues: {n}; "
              f"manifest FIGURES_MANIFEST.json")
        if m["skipped"]:
            print(f"[figures] skipped: {m['skipped']}")


if __name__ == "__main__":
    main()
