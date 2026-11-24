#!/usr/bin/env python3
"""cpaudit v8.1 command line (Phase 6 on JARVIS-DFT: Tasks 6.2-6.3 data, 6.4 models and UQ, 6.5 analyses).

  python run.py setup8       checks v7 and writes the Task 6.5 pre-registration ADDENDUM (before any v8 result)
  python run.py scores8      distance and failure-detection scores per outer fold   [--target T ...]
  python run.py fewshot8     H5-J: few-shot recalibration of the 3 worst-covered groups per target
  python run.py hjplus8      hierarchical jackknife+ (Lee et al. eq. 20), resumable  [--target T ...]
  python run.py report8      Task 6.5 tables and the decisions H4-J / G4-J / H5-J (+ descriptive parts)
  python run.py figures8     Fig. 8-10 (JARVIS) into paper_figures_v3/

  v8.1 (POST HOC, not pre-registered; reads the v7 and v8 outputs, writes new folders only):
  python run.py posthoc8     few-shot vs a no-model interval (paired with the fewshot8 draws), inner-model SD
                             at group level, ABO3 oxides, degenerate group-level intervals -> phase08_v1/posthoc/
  python run.py figures81    Fig. 8-10 with the v8.1 presentation fixes -> paper_figures_v4/

  python run.py setup7       design + input hashes + timing + PRE-REGISTRATION (before any fit)
  python run.py tune7        cross-half nested tuning (resumable per trial)   [--target T ...]
  python run.py fit7         all model units, resumable           [--target T ...] [--block B ...]
  python run.py report7      Task 6.4 tables and the decisions G5, H1-J, G3-J, H3-J

  python run.py selftest     all self-tests (v1-v5 tools plus the v6 data, feature and split tools)
  python run.py freeze6      rebuild the analysis table from the pinned release; must equal the audit
  python run.py features6    Magpie (145) + density (3) + CrystalNN (122) features; parallel, resumable
  python run.py splits6      LOSGO / LOCSO holdouts (MatFold labels primary, JARVIS labels as sensitivity)
                             + polymorph exclusion
  python run.py matfold6     run MatFold itself on all structures; its splits must equal ours row by row
  python run.py report6      one summary of freeze, features and splits

Earlier commands still available: posthoc5, audit6, figures. Steps up to report6 use no target and fit
no model; setup7 must run (and write the pre-registration) before tune7 / fit7.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cpaudit.config import load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "freeze6", "features6", "splits6", "matfold6", "report6",
                                        "setup7", "tune7", "fit7", "report7",
                                        "setup8", "scores8", "fewshot8", "hjplus8", "report8", "figures8",
                                        "posthoc8", "figures81",
                                        "posthoc5", "audit6", "figures"])
    ap.add_argument("--config", default=None)
    ap.add_argument("--only", nargs="*", default=None, help="figure name prefixes, e.g. fig7 figS")
    ap.add_argument("--target", nargs="*", default=None, help="formation bandgap tbmbj (default: all)")
    ap.add_argument("--block", nargs="*", default=None,
                    help="primary regimeA randinner polyex h3 rep2 locso (default: all, in this order)")
    ap.add_argument("--allow-deviation", action="store_true",
                    help="run although config, code or inputs changed since the pre-registration (logged)")
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
    elif a.command == "setup7":
        from cpaudit.setup7 import run_setup7
        run_setup7(cfg)
    elif a.command == "tune7":
        from cpaudit.fit7 import run_tune7
        run_tune7(cfg, a.target, a.allow_deviation)
    elif a.command == "fit7":
        from cpaudit.fit7 import run_fit7
        run_fit7(cfg, a.target, a.block, a.allow_deviation)
    elif a.command == "report7":
        from cpaudit.report7 import run_report7
        run_report7(cfg)
    elif a.command == "setup8":
        from cpaudit.setup8 import run_setup8
        run_setup8(cfg)
    elif a.command == "scores8":
        from cpaudit.run8 import run_scores8
        run_scores8(cfg, a.target, a.allow_deviation)
    elif a.command == "fewshot8":
        from cpaudit.run8 import run_fewshot8
        run_fewshot8(cfg, a.allow_deviation)
    elif a.command == "hjplus8":
        from cpaudit.run8 import run_hjplus8
        run_hjplus8(cfg, a.target, a.allow_deviation)
    elif a.command == "report8":
        from cpaudit.report8 import run_report8
        run_report8(cfg, a.allow_deviation)
    elif a.command == "figures8":
        from cpaudit.figs8 import run_figures8
        run_figures8(cfg)
    elif a.command == "posthoc8":
        from cpaudit.posthoc8 import run_posthoc8
        run_posthoc8(cfg)
    elif a.command == "figures81":
        from cpaudit.figs81 import run_figures81
        run_figures81(cfg)
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
