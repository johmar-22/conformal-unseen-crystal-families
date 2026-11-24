#!/usr/bin/env python3
"""cpaudit v6.1 command line (Phase 6 Tasks 6.2-6.3: JARVIS-DFT freeze, features, splits).

  python run.py selftest     all self-tests (v1-v5 tools plus the v6 data, feature and split tools)
  python run.py freeze6      rebuild the analysis table from the pinned release; must equal the audit
  python run.py features6    Magpie (145) + density (3) + CrystalNN (122) features; parallel, resumable
  python run.py splits6      LOSGO / LOCSO holdouts (MatFold labels primary, JARVIS labels as sensitivity)
                             + polymorph exclusion
  python run.py matfold6     run MatFold itself on all structures; its splits must equal ours row by row
  python run.py report6      one summary of freeze, features and splits

Earlier commands still available: posthoc5, audit6, figures. No target is used and no model is fitted.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cpaudit.config import load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "freeze6", "features6", "splits6", "matfold6", "report6",
                                        "posthoc5", "audit6", "figures"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--only", nargs="*", default=None, help="figure name prefixes, e.g. fig7 figS")
    a = ap.parse_args()
    if a.command == "selftest":
        from cpaudit.selftest import run_selftest
        sys.exit(run_selftest())
    cfg = load_config(a.config)
    if a.command == "freeze6":
        from cpaudit.data6 import run_freeze6
        run_freeze6(cfg)
    elif a.command == "features6":
        from cpaudit.features6 import run_features6
        run_features6(cfg)
    elif a.command == "splits6":
        from cpaudit.splits6 import run_splits6
        run_splits6(cfg)
    elif a.command == "matfold6":
        from cpaudit.matfold6 import run_matfold6
        run_matfold6(cfg)
    elif a.command == "report6":
        from cpaudit.report6 import run_report6
        run_report6(cfg)
    elif a.command == "posthoc5":
        from cpaudit.posthoc5 import run_posthoc5
        run_posthoc5(cfg)
    elif a.command == "audit6":
        from cpaudit.audit6 import run_audit6
        run_audit6(cfg)
    elif a.command == "figures":
        from cpaudit.paperfigs import run_figures
        m = run_figures(cfg, a.only)
        print(f"[figures] {len(m['figures'])} drawn, {len(m['skipped'])} skipped; overlap issues: "
              f"{m['overlap_check']['n_issues']}; manifest FIGURES_MANIFEST.json")


if __name__ == "__main__":
    main()
