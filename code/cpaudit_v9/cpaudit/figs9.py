"""Paper figures, final set: 8 main figures and 7 supplementary figures, drawn only from the result files on
the mounted Drive (nothing is refitted or recomputed beyond plotting summaries).

  Fig. 1   study design: split designs (a), family sizes (b), families needed for a bounded group-level
           interval (c), workflow with the points where the rules were registered (d)
  Fig. 2   regime A, NOMAD: coverage by lattice (a-d) and group-level Mondrian (e-g)
  Fig. 3   regime B, NOMAD: coverage of the held-out lattice by method
  Fig. 4   mechanism, NOMAD: offset and transfer (a-d), no covariate overlap (e, f)
  Fig. 5   failure detection, NOMAD (a-c) and JARVIS-DFT (d-f), side by side, never pooled
  Fig. 6   few-shot recalibration, NOMAD (a-f) and the dispersion test (g, h)
  Fig. 7   JARVIS-DFT unseen space groups (a-d) and group-level guarantees (e, f)
  Fig. 8   JARVIS-DFT few-shot recalibration of the three worst-covered groups per target
  Fig. S1-S7  calibration tail, formation transfer, detection AUROC matrix, lattice-level detection,
           sub-family transfer and training arm, JARVIS few-shot per group, hierarchical jackknife+ vs CV+

Composite figures reuse the drawing code of the single figures (paperfigs, figs8, figs81): each part is drawn
into a sub-figure of one page-sized figure, its panel letters are re-assigned, and the whole figure goes
through the automatic overlap check (figcheck) before it is saved. Every file read while drawing (through
pandas.read_csv, numpy.load or json.load) is recorded with its SHA-256 in FIGURES9_MANIFEST.json.

Style: figures4papers house style (Chen Liu, CC BY-NC 4.0), adapted in cpaudit/pubstyle.py.
"""
from __future__ import annotations

import contextlib
import json
import os
import re
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from . import figs8 as F8
from . import figs81 as F81
from . import paperfigs as P
from . import pubstyle as PS
from .an8 import sg_tex
from .config import sha256_file
from .pubstyle import DOUBLE_COL_IN, LATTICE_ORDER, LATTICE_TEX, PALETTE, panel_label

import matplotlib.pyplot as plt  # noqa: E402  (pubstyle selects the Agg backend first)
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Rectangle  # noqa: E402

ALPHA = 0.10
NOMINAL = 1 - ALPHA
PAGE_MAX_IN = 9.2          # tallest figure that fits a journal page with its caption

# one encoding for the detection scores in both datasets (Fig. 5, S3)
SCORES = [("width", "CV+ interval width", PALETTE["red_strong"]),
          ("knn", "kNN distance", PALETTE["blue_main"]),
          ("kde", "KDE dissimilarity", PALETTE["teal"]),
          ("committee", "Inner-model SD", PALETTE["violet"]),
          ("subsample", "Subsample-ensemble SD", PALETTE["magenta"]),
          ("classifier", "Domain classifier", PALETTE["grey_mid"])]
SCORE_COL = {k: c for k, _, c in SCORES}
TMARK = {"bandgap": "o", "formation": "s"}            # marker = target, everywhere in Fig. 5
GROUP_COLORS = (PALETTE["blue_main"], PALETTE["teal"], PALETTE["violet"])
GROUP_MARKERS = ("o", "s", "D")


# ------------------------------------------------------------------------------------------ inputs
class Recorder:
    """SHA-256 of every file read while the figures are drawn."""

    def __init__(self):
        self.used: dict[str, str] = {}

    def add(self, p):
        if isinstance(p, (str, os.PathLike)) and os.path.isfile(p):
            p = os.path.abspath(p)
            if p not in self.used:
                self.used[p] = sha256_file(p)


@contextlib.contextmanager
def recording(rec: Recorder):
    o_csv, o_np, o_js = pd.read_csv, np.load, json.load

    def rc(p, *a, **k):
        rec.add(p)
        return o_csv(p, *a, **k)

    def nl(p, *a, **k):
        rec.add(p)
        return o_np(p, *a, **k)

    def jl(fh, *a, **k):
        rec.add(getattr(fh, "name", None))
        return o_js(fh, *a, **k)

    pd.read_csv, np.load, json.load = rc, nl, jl
    try:
        yield rec
    finally:
        pd.read_csv, np.load, json.load = o_csv, o_np, o_js


class Dirs:
    def __init__(self, root, p1, p2, p4, p5, p6, d7, d8, posthoc):
        self.root, self.d7, self.d8, self.posthoc = root, d7, d8, posthoc
        self.inp = P.Inputs(p1, p2, p4, p5, p6)
        self.counted: dict[str, list] = {}    # files whose number (not content) enters a figure

    def f8(self, *parts):
        return os.path.join(self.d8, *parts)

    def need(self, *paths):
        miss = [p for p in paths if not os.path.exists(p)]
        if miss:
            raise _Skip(f"missing {miss[0]}" + (f" (+{len(miss) - 1} more)" if len(miss) > 1 else ""))


class _Skip(Exception):
    pass


# ------------------------------------------------------------------------------------------ composition
_MODS = (P, F8, F81)
_REF = re.compile(r"\(([a-h](?:\s*[,\-]\s*[a-h])*)\)")
_PREFIX = re.compile(r"^([a-h]): ")


def _remap_refs(s: str, mp: dict) -> str:
    """'(a, b)' -> '(e, f)', '(c-f)' -> '(g-j)', 'a: ...' -> 'e: ...' with the letter map mp."""
    def one(m):
        body = m.group(1)
        if "-" in body:
            a, b = [x.strip() for x in body.split("-")]
            letters = [chr(c) for c in range(ord(a), ord(b) + 1)]
            new = [mp.get(x, x) for x in letters]
            ok = all(ord(new[i + 1]) - ord(new[i]) == 1 for i in range(len(new) - 1))
            return f"({new[0]}-{new[-1]})" if ok else "(" + ", ".join(new) + ")"
        return "(" + ", ".join(mp.get(x.strip(), x.strip()) for x in body.split(",")) + ")"
    s = _REF.sub(one, s)
    return _PREFIX.sub(lambda m: f"{mp.get(m.group(1), m.group(1))}: ", s)


def _texts_of(sub):
    for leg in sub.legends:
        yield from leg.get_texts()
    for ax in sub.axes:
        leg = ax.get_legend()
        if leg is not None:
            yield from leg.get_texts()
        yield from ax.texts
        yield ax.title
        yield ax._left_title


def reletter(sub, mp: dict | None):
    if not mp:
        return
    for ax in sub.axes:
        l_ = getattr(ax, "_cp_panel", None)
        if l_ in mp:
            ax._cp_panel = mp[l_]
    for t in _texts_of(sub):
        s = t.get_text()
        n = _remap_refs(s, mp)
        if n != s:
            t.set_text(n)


def _root(fig):
    while hasattr(fig, "_parent") and fig._parent is not None and fig._parent is not fig:
        fig = fig._parent
    return fig


@contextlib.contextmanager
def hosted(sub):
    """Make the single-figure drawing code of paperfigs / figs8 / figs81 draw into the sub-figure `sub`
    and return without saving."""
    root = _root(sub)
    sub.draw_without_rendering = root.draw_without_rendering
    old = [(m, m.new_figure, m.finalize_figure) for m in _MODS]

    def nf(*a, wspace=None, hspace=None, **k):
        # the spacing a single figure asked for (new_figure(..., wspace=, hspace=)) is applied to the grids of
        # its sub-figure, since the layout engine (and its spacing) belongs to the whole figure
        if wspace is not None or hspace is not None:
            ags, sp = sub.add_gridspec, sub.subplots

            def add_gridspec(*ga, **gk):
                gk.setdefault("wspace", wspace)
                gk.setdefault("hspace", hspace)
                return ags(*ga, **gk)

            def subplots(*sa, gridspec_kw=None, **sk):
                g = dict(gridspec_kw or {})
                g.setdefault("wspace", wspace)
                g.setdefault("hspace", hspace)
                return sp(*sa, gridspec_kw=g, **sk)

            sub.add_gridspec, sub.subplots = add_gridspec, subplots
        return sub

    def ff(fig, base, *a, **k):
        return [str(base)]

    for m in _MODS:
        m.new_figure, m.finalize_figure = nf, ff
    try:
        yield sub
    finally:
        for m, a, b in old:
            m.new_figure, m.finalize_figure = a, b


def old_part(fn, *args, mapping=None, **kw):
    """A part that draws an existing single figure into a sub-figure, then re-letters its panels."""
    def draw(sub):
        with hosted(sub):
            saved, why = fn(*args, **kw)
        if not saved:
            raise _Skip(f"{fn.__name__}: {why}")
        reletter(sub, mapping)
    return draw


def composite(out, name, parts, width=DOUBLE_COL_IN, order=None, hspace=0.02):
    """parts: list of (height in inches, draw(sub)). Parts are stacked top to bottom; `order` gives the
    drawing order (a part with leader-line labels is best drawn last, on the final layout)."""
    H = float(sum(h for h, _ in parts))
    if H > PAGE_MAX_IN + 1e-9:
        raise _Skip(f"{name}: {H:.1f} in is taller than a page ({PAGE_MAX_IN} in)")
    fig = PS.new_figure(width, H)
    subs = np.atleast_1d(fig.subfigures(len(parts), 1, height_ratios=[h for h, _ in parts], hspace=hspace))
    for i in (order or range(len(parts))):
        parts[i][1](subs[i])
    _check_inside(fig)
    return PS.finalize_figure(fig, os.path.join(out, name))


def _check_inside(fig):
    """No text may be cut by the figure edge (the overlap check compares items with each other only)."""
    fig.draw_without_rendering()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    from .figcheck import _texts
    for bb, kind, _, t in _texts(fig, r):
        if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > W + 0.5 or bb.y1 > H + 0.5:
            raise _Skip(f"{kind} '{t.get_text()[:40]}' is cut by the figure edge")


def single(out, name, w, h, draw):
    fig = PS.new_figure(w, h)
    draw(fig)
    _check_inside(fig)
    return PS.finalize_figure(fig, os.path.join(out, name))


# ------------------------------------------------------------------------------------------ Fig. 1
ROLES = {"fit": ("fitting and calibration (cross-fitted)", "#9DB9DD", None),
         "train": ("training only (model fixed)", "#DCE6F2", None),
         "cal": ("calibration labels", PALETTE["green_3"], None),
         "test": ("test (coverage measured)", PALETTE["red_2"], None),
         "out": ("not used", "#F2F2F2", "////")}
DESIGNS = [("Random split", "A", [[("fit", 0.8), ("test", 0.2)]] * 6),
           ("Leave one family out", "B", [[("fit", 1.0)]] * 5 + [[("test", 1.0)]]),
           ("Leave two out", "B", [[("fit", 1.0)]] * 4 + [[("test", 1.0)], [("out", 1.0)]]),
           ("Few labels, new family", "B", [[("train", 1.0)]] * 5 + [[("cal", 0.18), ("test", 0.82)]]),
           ("Group-level split", "B", [[("train", 1.0)]] * 2 + [[("cal", 1.0)]] * 3 + [[("test", 1.0)]])]


def _fig1_designs(ax):
    nf = 6
    for i, (lab, reg, cells) in enumerate(DESIGNS):
        for j, cell in enumerate(cells):
            x = j + 0.08
            for role, frac in cell:
                _, col, hatch = ROLES[role]
                w = 0.84 * frac
                r = Rectangle((x, i + 0.2), w, 0.6, facecolor=col, edgecolor="white", lw=0.8, hatch=hatch, zorder=2)
                r.set_gid("bg")
                ax.add_patch(r)
                x += w
    ax.axhline(0.98, color=PALETTE["grey_mid"], lw=0.7, ls=(0, (2, 2)))
    ax.set_xlim(0, nf)
    ax.set_ylim(len(DESIGNS), -0.05)
    ax.set_yticks(np.arange(len(DESIGNS)) + 0.5)
    ax.set_yticklabels([f"{lab} ({reg})" for lab, reg, _ in DESIGNS], fontsize=6.8)
    ax.set_xticks(np.arange(nf) + 0.5)
    ax.set_xticklabels([str(k + 1) for k in range(nf)], fontsize=6.8)
    ax.set_xlabel("Family (lattice or space group)")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    panel_label(ax, "a", "Split designs")


# Families needed for a bounded hierarchical jackknife+ interval: with N families, K = N - 1 leave-one-out
# models, and the interval is finite only if 1/(K + 1) <= alpha, i.e. N >= ceil(1/alpha).
NEEDED = {a: int(np.ceil(1 / a - 1e-9)) for a in (0.10, 0.05)}


def _nomad_sizes(D: Dirs) -> pd.Series:
    """Structures per lattice in the de-duplicated analysis set (leave-one-lattice-out test sizes)."""
    L = D.inp.csv("phase2", "report", "phase2_loso_all.csv")
    L = L[L.design == "loso"].drop_duplicates("heldout")
    tot = set((L.n_train + L.n_test).astype(int))
    if len(tot) != 1:
        raise _Skip(f"inconsistent NOMAD totals in phase2_loso_all.csv: {sorted(tot)}")
    s = L.set_index("heldout").n_test.astype(int).sort_values(ascending=False)
    if int(s.sum()) != tot.pop():
        raise _Skip("NOMAD lattice sizes do not add up to the analysis set")
    return s


def _fig1_sizes(ax, D: Dirs):
    nom = _nomad_sizes(D)
    meta = D.inp.json("phase6", "audit", "AUDIT_DECISION.json")
    rel, tmin = meta["decision"]["release"], int(meta["audit_settings"]["min_group_size"])
    SP = D.inp.csv("phase6", "audit", "audit_spacegroups.csv")
    G = SP[(SP.release == rel) & (SP.scope == "all_chemistries")].sort_values("n", ascending=False)
    rank = np.arange(1, len(G) + 1)
    ax.step(rank, G.n, where="mid", color=PALETTE["blue_main"], lw=1.3)
    e = G.eligible.to_numpy(bool)
    ax.scatter(rank[e], G.n.to_numpy()[e], s=7, color=PALETTE["blue_main"], lw=0, zorder=3)
    rk = np.arange(1, len(nom) + 1)
    ax.plot(rk, nom.to_numpy(), color=PALETTE["red_strong"], lw=1.0, marker="D", ms=3.6, zorder=4)
    ax.axhline(tmin, color=PALETTE["ink"], lw=0.8, ls=(0, (4, 3)))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Families, ranked by size")
    ax.set_ylabel("Structures per family")
    panel_label(ax, "b", "Family sizes")
    return dict(nomad_sizes=nom, nomad_n=int(nom.sum()), nomad_k=int(len(nom)), jarvis_n=int(G.n.sum()),
                jarvis_k=int(e.sum()), tmin=tmin)


def _fig1_needed(ax, D: Dirs, info):
    TH = D.inp.csv("phase6", "audit", "audit_thresholds.csv")
    rel = D.inp.json("phase6", "audit", "AUDIT_DECISION.json")["decision"]["release"]
    xs = None
    for sc, col, ls in (("all_chemistries", PALETTE["blue_main"], "-"), ("oxides_oxidation_rule", PALETTE["ochre"], (0, (3, 2)))):
        G = TH[(TH.release == rel) & (TH.scope == sc)].sort_values("min_group_size")
        ax.plot(G.min_group_size, G.K, color=col, ls=ls, lw=1.4, marker="o", ms=3)
        xs = G.min_group_size.to_numpy()
    nom = info["nomad_sizes"].to_numpy()
    ax.plot(xs, [int((nom >= t).sum()) for t in xs], color=PALETTE["red_strong"], lw=1.0, marker="D", ms=3.6, zorder=4)
    for k, ls in ((NEEDED[0.10], (0, (1, 1.5))), (NEEDED[0.05], (0, (5, 2)))):
        ax.axhline(k, color=PALETTE["grey_dark"], lw=0.8, ls=ls)
    ax.axvline(info["tmin"], color=PALETTE["ink"], lw=0.8, ls=(0, (4, 3)))
    step = 100 if xs.max() - xs.min() >= 300 else 50
    ax.set_xticks(np.arange(np.ceil(xs.min() / step) * step, xs.max() + 1, step))
    ax.set_xlim(xs.min() - 30, xs.max() + 30)
    ax.set_ylim(0, None)
    ax.set_xlabel("Size threshold\n(structures per family)")
    ax.set_ylabel("Eligible families")
    panel_label(ax, "c", "Families available")


# Fig. 1d: one pipeline row per dataset. Every count in a box is read from a result file (_workflow_facts);
# the method lists are built from the method codes that occur in the result files, so a method is named only
# if results for it exist.
STEPS = ("Data", "Representation", "Model", "Uncertainty", "Evaluation", "Analyses")
METHOD_NAMES = [  # (codes found in the result files, name shown), in display order
    ({"cvplus_random", "cvplus_lattice", "cvplus"}, "CV+"),
    ({"cvplus_gw"}, "group-weighted CV+"),
    ({"split_global"}, "split"),
    ({"split_mondrian"}, "Mondrian"),
    ({"split_normalized_knn", "ncvplus_knn", "ncvplus"}, "normalized"),
    ({"cqr", "cqr_cvplus"}, "CQR"),
    ({"group_level"}, "group-level"),
    ({"hcp"}, "HCP"),
    ({"single_subsample"}, "single subsampling"),
    ({"hjplus"}, "HJ+")]


def _method_list(codes: set) -> list:
    return [name for keys, name in METHOD_NAMES if keys & codes]


def _wrap(items, width=17):
    import textwrap
    return textwrap.fill(", ".join(items), width=width, break_long_words=False, break_on_hyphens=False)


def _workflow_facts(D: Dirs, info) -> dict:
    inp = D.inp
    P1 = inp.csv("phase1", "report", "phase1_all_runs.csv")
    LO = inp.csv("phase2", "report", "phase2_loso_all.csv")
    L2 = inp.csv("phase2", "report", "phase2_l2o_all.csv")
    nom_codes = set(P1.method) | set(LO.method) | set(L2.method)
    feat = inp.json("phase6", "features", "features_v1.json")
    frz = inp.json("phase6", "data", "freeze.json")
    blocks = feat["blocks"]
    if sum(blocks.values()) != int(feat["shape"][1]) or int(feat["shape"][0]) != int(frz["n_rows"]):
        raise _Skip("JARVIS feature blocks do not add up to the feature matrix")
    if int(frz["n_rows"]) != info["jarvis_n"]:
        raise _Skip(f"JARVIS frozen table ({frz['n_rows']}) differs from the audit ({info['jarvis_n']})")
    files = sorted(f for f in os.listdir(os.path.join(D.d7, "unitsB")) if f.endswith(".csv"))
    prim = [f for f in files if "_primary_" in f]
    folds = sorted({f.rsplit("_o", 1)[1][:2] for f in prim})
    jar_codes = set()
    for f in prim:
        jar_codes |= set(pd.read_csv(os.path.join(D.d7, "unitsB", f), usecols=["method"]).method)
    jar_codes.discard("enssd")              # MatFold's inner-model SD interval: a baseline, not a conformal method
    for f in os.listdir(os.path.join(D.d7, "unitsH3")):
        if f.endswith(".csv"):
            jar_codes |= set(pd.read_csv(os.path.join(D.d7, "unitsH3", f), usecols=["method"]).method)
    if all(os.path.exists(D.f8("report", f"{t}_hjplus_summary.csv")) for t in ("formation", "bandgap")):
        if all(len(pd.read_csv(D.f8("report", f"{t}_hjplus_summary.csv"))) for t in ("formation", "bandgap")):
            jar_codes.add("hjplus")
    # regime A units are counted, not read (their content is not used in this figure); listed in the manifest
    dA = os.path.join(D.d7, "unitsA")
    a_units = sorted(f for f in os.listdir(dA) if re.fullmatch(r"(formation|bandgap)_A_s\d\d\.csv", f)) \
        if os.path.isdir(dA) else []
    n_a = {t: len([f for f in a_units if f.startswith(t + "_")]) for t in ("formation", "bandgap")}
    D.counted["phase07 unitsA (regime A splits)"] = a_units
    a_line = f"A: {n_a['formation']} random splits\n" if a_units and n_a["formation"] == n_a["bandgap"] else ""
    k_fe = len(pd.read_csv(D.f8("report", "formation_groups.csv")))
    k_bg = len(pd.read_csv(D.f8("report", "bandgap_groups.csv")))
    return dict(
        nomad=[f"NOMAD 2018\n{info['nomad_n']:,} structures\n{info['nomad_k']} lattices",
               "SOAP\ndescriptors",
               "LightGBM,\nOptuna-tuned",
               _wrap(_method_list(nom_codes)),
               f"A: {P1.repeat.nunique()} random splits\nB: leave one lattice\nout ({LO.heldout.nunique()}), leave\n"
               f"two out ({L2.heldout.nunique()})",
               "mechanism,\nfailure detection,\nfew-shot\nrecalibration"],
        jarvis=[f"JARVIS-DFT\n{int(frz['n_rows']):,} structures\n{k_fe} space groups\n({k_bg} for band gap)",
                f"Magpie {blocks['composition']},\ndensity {blocks['density']},\nCrystalNN {blocks['crystalnn']}",
                "LightGBM,\nOptuna-tuned",
                _wrap(_method_list(jar_codes)),
                a_line + f"B: {len(folds)} folds of\nheld-out space\ngroups",
                "failure detection,\ngroup-level\nguarantees,\nfew-shot\nrecalibration"],
        k_fe=k_fe, k_bg=k_bg)


def _fig1_workflow(ax, facts):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    n = len(STEPS)
    slot = 1.0 / n
    bw, bh = slot * 0.84, 0.36
    rows = (("nomad", 0.66), ("jarvis", 0.22))
    boxes = []
    for j, name in enumerate(STEPS):
        ax.text(slot * (j + 0.5), 0.97, name, ha="center", va="top", fontsize=6.8, fontweight="bold",
                color=PALETTE["ink"])
    for key, yc in rows:
        for j, txt in enumerate(facts[key]):
            xc = slot * (j + 0.5)
            b = FancyBboxPatch((xc - bw / 2, yc - bh / 2), bw, bh, boxstyle="round,pad=0,rounding_size=0.012",
                               facecolor="#DCE6F2" if j else "#F2F2F2", edgecolor=PALETTE["blue_main"] if j else
                               PALETTE["grey_mid"], lw=0.9, zorder=2)
            b.set_gid("bg")
            ax.add_patch(b)
            t_ = ax.text(xc, yc, txt, ha="center", va="center", fontsize=5.9, linespacing=1.12, zorder=4,
                         color=PALETTE["ink"])
            boxes.append((b, t_))
            if j < n - 1:
                a = FancyArrowPatch((xc + bw / 2, yc), (xc + slot - bw / 2, yc), arrowstyle="-|>", mutation_scale=6,
                                    color=PALETTE["grey_dark"], lw=0.8, zorder=3, shrinkA=0.5, shrinkB=0.5)
                a.set_gid("bg")
                ax.add_patch(a)
    panel_label(ax, "d", "Workflow")
    return boxes


def _check_boxes(fig, boxes):
    """Every workflow text must lie inside its box (the overlap check treats the boxes as background)."""
    fig.draw_without_rendering()
    r = fig.canvas.get_renderer()
    bad = [t.get_text().replace("\n", " ") for b, t in boxes
           if not b.get_window_extent(r).expanded(0.93, 0.90).contains(*t.get_window_extent(r).p0)
           or not b.get_window_extent(r).expanded(0.93, 0.90).contains(*t.get_window_extent(r).p1)]
    if bad:
        raise _Skip(f"workflow text wider than its box: {bad[0]!r}")


def fig1_study_design(D: Dirs, out):
    D.need(D.inp.path("phase2", "report", "phase2_loso_all.csv"), D.inp.path("phase2", "report", "phase2_l2o_all.csv"),
           D.inp.path("phase1", "report", "phase1_all_runs.csv"), D.inp.path("phase6", "features", "features_v1.json"),
           D.inp.path("phase6", "data", "freeze.json"), D.f8("report", "formation_groups.csv"),
           D.f8("report", "bandgap_groups.csv"), os.path.join(D.d7, "unitsB"), os.path.join(D.d7, "unitsH3"), D.inp.path("phase6", "audit", "audit_spacegroups.csv"),
           D.inp.path("phase6", "audit", "audit_thresholds.csv"), D.inp.path("phase6", "audit", "AUDIT_DECISION.json"))

    def draw(fig):
        gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.95], width_ratios=[1.45, 1.0, 1.0])
        _fig1_designs(fig.add_subplot(gs[0, 0]))
        info = _fig1_sizes(fig.add_subplot(gs[0, 1]), D)
        _fig1_needed(fig.add_subplot(gs[0, 2]), D, info)
        boxes = _fig1_workflow(fig.add_subplot(gs[1, :]), _workflow_facts(D, info))
        h = [Patch(facecolor=c, edgecolor="white", hatch=hh, label=f"{lab} (a)") for lab, c, hh in ROLES.values()]
        h += [Line2D([], [], color=PALETTE["blue_main"], lw=1.3, marker="o", ms=3, label="JARVIS-DFT, all chemistries (b, c)"),
              Line2D([], [], color=PALETTE["ochre"], lw=1.4, ls=(0, (3, 2)), marker="o", ms=3,
                     label="JARVIS-DFT, oxides only (c)"),
              Line2D([], [], color=PALETTE["red_strong"], lw=1.0, marker="D", ms=3.6,
                     label=f"NOMAD 2018, {info['nomad_k']} lattices (b, c)"),
              Line2D([], [], color=PALETTE["ink"], lw=0.8, ls=(0, (4, 3)), label=f"size threshold {info['tmin']} (b, c)"),
              Line2D([], [], color=PALETTE["grey_dark"], lw=0.8, ls=(0, (1, 1.5)),
                     label=r"bounded group-level interval, $\alpha$ = 0.10: $\geq$ " + f"{NEEDED[0.10]} families (c)"),
              Line2D([], [], color=PALETTE["grey_dark"], lw=0.8, ls=(0, (5, 2)),
                     label=r"$\alpha$ = 0.05: $\geq$ " + f"{NEEDED[0.05]} families (c)")]
        fig.legend(handles=h, loc="outside lower center", ncol=3, fontsize=6.4, columnspacing=1.0, handletextpad=0.5)
        _check_boxes(fig, boxes)
    return single(out, "fig1_study_design", DOUBLE_COL_IN, 5.9, draw)


# ------------------------------------------------------------------------------------------ Fig. 5
def _forest(ax, rows, data, title, letter):
    """rows: score keys; data: {(target, level): DataFrame indexed by score with auroc, lo, hi}."""
    offs = {("bandgap", "pooled"): -0.27, ("bandgap", "within"): -0.09, ("formation", "pooled"): 0.09,
            ("formation", "within"): 0.27}
    for i, s in enumerate(rows):
        if s in ("knn", "kde"):
            b = ax.axhspan(i - 0.48, i + 0.48, color=PALETTE["highlight"], alpha=0.13, lw=0, zorder=0)
            b.set_gid("bg")
        for (t, lev), d in data.items():
            if s not in d.index or not np.isfinite(d.loc[s, "auroc"]):
                continue
            y = i + offs[(t, lev)]
            a, lo, hi = (float(d.loc[s, k]) for k in ("auroc", "lo", "hi"))
            if np.isfinite(lo) and np.isfinite(hi):
                ax.plot([lo, hi], [y, y], color=SCORE_COL[s], lw=1.2, solid_capstyle="round", zorder=2)
            ax.scatter(a, y, s=20, marker=TMARK[t], facecolor=SCORE_COL[s] if lev == "pooled" else "white",
                       edgecolor=SCORE_COL[s], lw=1.1, zorder=3)
    ax.axvline(0.5, color=PALETTE["ink"], lw=0.9)
    ax.axvline(0.7, color=PALETTE["grey_mid"], lw=0.9, ls=(0, (2, 2)))
    lab = {k: l for k, l, _ in SCORES}
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab[s] for s in rows], fontsize=6.6)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("AUROC for a missed interval")
    panel_label(ax, letter, title)


def _selective(ax, curves, title, letter, ylab=True):
    """curves: {target: DataFrame(score, abstained, coverage_retained)}."""
    lo = 1.0
    for t, A in curves.items():
        ls = "-" if t == "bandgap" else (0, (3, 1.6))
        for s in ("width", "knn", "committee"):
            S = A[A.score == s].sort_values("abstained")
            if S.empty:
                continue
            ax.plot(S.abstained, S.coverage_retained, color=SCORE_COL[s], lw=1.3, ls=ls, marker=TMARK[t], ms=2.2,
                    mew=0)
            lo = min(lo, S.coverage_retained.min())
        S = A[A.score == "oracle_abs_err"].sort_values("abstained")
        if len(S):
            ax.plot(S.abstained, S.coverage_retained, color=PALETTE["ink"], lw=0.9, ls=ls, alpha=0.6)
    ax.axhline(NOMINAL, color=PALETTE["ink"], lw=0.8, ls=(0, (4, 3)))
    ax.set_xlabel("Fraction abstained")
    if ylab:
        ax.set_ylabel("Coverage of the retained")
    ax.set_ylim(max(0.0, np.floor(lo * 20) / 20 - 0.03), 1.01)
    panel_label(ax, letter, title)


def _nomad_flag(ax, Fl, letter):
    """One row per held-out lattice: leave one lattice out (circle) above, refit without Fd-3m (triangle) below."""
    for i, lat in enumerate(LATTICE_ORDER):
        for design, dy, mk in (("loso", -0.18, "o"), ("l2o", 0.18, "^")):
            r = Fl[(Fl.design == design) & (Fl.heldout == lat)]
            if r.empty:
                continue
            r = r.iloc[0]
            y = i + dy
            ax.plot([r.ffr_ci_lo, r.ffr_ci_hi], [y, y], color=PALETTE["grey_mid"], lw=1.3, zorder=1)
            ax.scatter(r.false_flag_rate, y, s=12, color=PALETTE["grey_mid"], zorder=2, clip_on=False)
            ax.scatter(r.flag_rate_heldout, y, s=26, color=PALETTE["ink"], marker=mk, edgecolor="white", lw=0.5,
                       zorder=3, clip_on=False)
    ax.axvline(0.05, color=PALETTE["red_strong"], lw=0.9, ls=(0, (4, 3)))
    ax.set_yticks(range(len(LATTICE_ORDER)))
    ax.set_yticklabels([LATTICE_TEX[h] for h in LATTICE_ORDER], fontsize=6.6)
    ax.set_ylim(len(LATTICE_ORDER) - 0.5, -0.5)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Flag rate\n" + r"(conformal $p \leq 0.05$)")
    panel_label(ax, letter, "NOMAD: outlier flag")


def _jarvis_group_sd(fig, ax, D: Dirs, letter):
    C = pd.read_csv(os.path.join(D.posthoc, "committee_groups.csv"))
    Sp = pd.read_csv(os.path.join(D.posthoc, "committee_spearman.csv"))
    W = pd.read_csv(D.f8("fewshot", "worst_groups.csv"))
    X, Y, S, items = [], [], [], []
    for t in ("bandgap", "formation"):
        d = C[C.target == t]
        n = len(d)
        x = d.committee_rank.to_numpy(float)
        ax.scatter(x, d.coverage_r1, s=13, marker=TMARK[t], color=SCORE_COL["committee"], lw=0,
                   alpha=1.0 if t == "bandgap" else 0.75, zorder=3)
        X += list(x)
        Y += list(d.coverage_r1)
        S += [20.0] * n          # placement keeps labels this far from every marker
        for _, w in W[W.target == t].iterrows():
            r = d[d.group == w.group]
            if len(r):
                items.append((float(r.committee_rank.iloc[0]), float(r.coverage_r1.iloc[0]),
                              sg_tex(w.symbol)))
    ax.axhline(NOMINAL, color=PALETTE["ink"], lw=0.8, ls=(0, (4, 3)))
    ax.axhline(0.80, color=PALETTE["red_strong"], lw=0.7, ls=(0, (1, 2)))
    ax.set_xscale("log")
    ax.set_xlim(0.75, max(X) * 1.35)
    ax.set_xticks([1, 2, 5, 10, 20, 50])
    ax.set_xticklabels(["1", "2", "5", "10", "20", "50"])
    ax.minorticks_off()
    ax.set_ylim(min(0.5, min(Y) - 0.05), 1.02)
    ax.set_xlabel("Inner-model SD rank\n(1 = highest)")
    ax.set_ylabel("Group coverage (CV+)")
    rho = {r.target: r for _, r in Sp[Sp.score == "committee_median"].iterrows()}
    sub = ", ".join(f"{'FE' if t == 'formation' else 'BG'} $\\rho$ = {rho[t].rho:.2f}" for t in ("formation", "bandgap")
                    if t in rho)
    panel_label(ax, letter, f"JARVIS groups (post hoc)\n{sub}")
    items.sort(key=lambda it: it[1])        # lowest coverage first: those labels have the least room
    return (ax, np.array(X), np.array(Y), np.array(S), items)


def fig5_detection(D: Dirs, out):
    inp = D.inp
    rep = ("phase4", "report")
    D.need(*[inp.path(*rep, n) for n in ("phase4_detection_metrics.csv", "phase4_abstention.csv", "phase4_flags.csv")],
           *[D.f8("report", f"{t}_{k}.csv") for t in ("formation", "bandgap") for k in ("h4_pooled", "h4_within_summary",
                                                                                         "h4_abstention")],
           os.path.join(D.posthoc, "committee_groups.csv"), os.path.join(D.posthoc, "committee_spearman.csv"))
    Mt = inp.csv(*rep, "phase4_detection_metrics.csv")
    Ab = inp.csv(*rep, "phase4_abstention.csv")
    Fl = inp.csv(*rep, "phase4_flags.csv")

    def draw(fig):
        gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1])
        # NOMAD
        data = {}
        for t in ("bandgap", "formation"):
            for lev, key in (("pooled", "pooled"), ("within", "within (n-weighted)")):
                d = Mt[(Mt.target == t) & (Mt.outcome == "miscovered") & (Mt.level == key)].set_index("score")
                data[(t, lev)] = d.rename(columns={"auroc_lo": "lo", "auroc_hi": "hi"})
        _forest(fig.add_subplot(gs[0, 0]), [k for k, _, _ in SCORES], data, "NOMAD: failure detection", "a")
        _selective(fig.add_subplot(gs[0, 1]), {t: Ab[Ab.target == t] for t in ("bandgap", "formation")},
                   "NOMAD: abstention", "b")
        _nomad_flag(fig.add_subplot(gs[0, 2]), Fl, "c")
        # JARVIS
        data = {}
        for t in ("bandgap", "formation"):
            Pp = pd.read_csv(D.f8("report", f"{t}_h4_pooled.csv")).set_index("score")
            Ww = pd.read_csv(D.f8("report", f"{t}_h4_within_summary.csv")).set_index("score").rename(
                columns={"weighted_auroc": "auroc"})
            data[(t, "pooled")], data[(t, "within")] = Pp, Ww
        _forest(fig.add_subplot(gs[1, 0]), [k for k, _, _ in SCORES if k != "subsample"], data,
                "JARVIS-DFT: failure detection", "d")
        _selective(fig.add_subplot(gs[1, 1]), {t: pd.read_csv(D.f8("report", f"{t}_h4_abstention.csv"))
                                               for t in ("bandgap", "formation")}, "JARVIS-DFT: abstention", "e")
        job = _jarvis_group_sd(fig, fig.add_subplot(gs[1, 2]), D, "f")
        h = [Line2D([], [], color=c, lw=1.6, label=l) for _, l, c in SCORES]
        h += [Patch(facecolor=PALETTE["highlight"], alpha=0.25, edgecolor="none", label="distance scores (a, d)"),
              Line2D([], [], marker="o", color=PALETTE["grey_dark"], lw=1.3, ms=4, label="band gap (solid in b, e)"),
              Line2D([], [], marker="s", color=PALETTE["grey_dark"], lw=1.3, ls=(0, (3, 1.6)), ms=4,
                     label="formation energy (dashed in b, e)"),
              Line2D([], [], marker="o", color=PALETTE["grey_dark"], mfc="white", lw=0, ms=4.5,
                     label="filled: pooled; open: within family (a, d)"),
              Line2D([], [], color=PALETTE["ink"], lw=0.9, alpha=0.6, label="oracle: abstain on true |error| (b, e)"),
              Line2D([], [], color=PALETTE["ink"], lw=0.9, label="chance, 0.5 (a, d)"),
              Line2D([], [], color=PALETTE["grey_mid"], lw=0.9, ls=(0, (2, 2)), label="decision threshold 0.70, fixed in advance (a, d)"),
              Line2D([], [], color=PALETTE["ink"], lw=0.8, ls=(0, (4, 3)), label=r"nominal $1-\alpha$ = 0.90 (b, e, f)"),
              Line2D([], [], color=PALETTE["red_strong"], lw=0.7, ls=(0, (1, 2)), label="coverage 0.80 (f)"),
              Line2D([], [], marker="o", color=PALETTE["ink"], lw=0, ms=4.5, label="flag rate, held-out lattice (c)"),
              Line2D([], [], marker="^", color=PALETTE["ink"], lw=0, ms=4.5,
                     label="flag rate, refit without " + LATTICE_TEX["Fd-3m"] + " (c)"),
              Line2D([], [], marker="o", color=PALETTE["grey_mid"], lw=1.3, ms=3.5,
                     label="false-flag rate in distribution, 95 % CI (c)"),
              Line2D([], [], color=PALETTE["red_strong"], lw=0.9, ls=(0, (4, 3)), label="nominal false-flag rate 0.05 (c)")]
        fig.legend(handles=h, loc="outside lower center", ncol=4, fontsize=6.3, columnspacing=1.0, handletextpad=0.4)
        F8._place_labels(fig, *job, PALETTE["red_strong"], hlines=(NOMINAL, 0.80), near=False)
    return single(out, "fig5_failure_detection", DOUBLE_COL_IN, 6.4, draw)


# ------------------------------------------------------------------------------------------ Fig. 7 e, f
def _jarvis_guarantees(sub, D: Dirs, alpha):
    gs = sub.add_gridspec(1, 2)
    lows, axes = [], []
    for j, t in enumerate(("formation", "bandgap")):
        ax = sub.add_subplot(gs[0, j])
        axes.append(ax)
        files = sorted(f for f in os.listdir(os.path.join(D.d7, "unitsH3")) if f.startswith(f"{t}_H3_") and f.endswith(".csv"))
        if not files:
            raise _Skip(f"no {t} H3 units in {D.d7}/unitsH3")
        H = pd.concat([pd.read_csv(os.path.join(D.d7, "unitsH3", f)) for f in files], ignore_index=True)
        xt, xl, x = [], [], 0
        for mk, lab, col in (("hcp", "HCP", PALETTE["blue_main"]), ("single_subsample", "single\nsubsampling", PALETTE["teal"])):
            gw = H[(H.method == mk) & np.isclose(H.alpha, alpha)].groupby(["partition", "reshuffle"]).coverage.mean()
            q = np.quantile(gw, [0.05, 0.95])
            lows.append(q[0])
            ax.plot([x, x], q, color=col, lw=2.4, solid_capstyle="butt")
            ax.plot(x, gw.mean(), "o", color=col, ms=4.5, mec=PALETTE["ink"], mew=0.4, zorder=4)
            xt.append(x)
            xl.append(lab)
            x += 1
        S = pd.read_csv(D.f8("report", f"{t}_hjplus_summary.csv"))
        s_ = S[np.isclose(S.alpha, alpha)]
        if len(s_):
            ax.plot(x, float(s_.group_weighted.iloc[0]), "D", color=PALETTE["violet"], ms=4.8, mec=PALETTE["ink"],
                    mew=0.4, zorder=4)
            xt.append(x)
            xl.append("hierarchical\njackknife+")
            x += 1
        ax.axhline(1 - alpha, color=PALETTE["ink"], lw=0.8, ls=(0, (4, 3)))
        ax.axhline(1 - 2 * alpha, color=PALETTE["violet"], lw=0.8, ls=(0, (1, 1.5)))
        ax.set_xticks(xt, xl, fontsize=6.6)
        ax.set_xlim(-0.6, x - 0.4)
        ax.set_ylabel("Group-weighted coverage" if j == 0 else "")
        panel_label(ax, "ef"[j], f"{F8.TNAME[t]}: group-level guarantees")
    lo_ = np.floor(min([0.78] + [v - 0.02 for v in lows]) * 20) / 20
    for a_ in axes:
        a_.set_ylim(lo_, 1.02)
    h = [Line2D([], [], color=PALETTE["grey_dark"], lw=2.4, marker="o", ms=4.5, mec=PALETTE["ink"],
                label="mean and 5-95 % over partitions \u00d7 reshuffles (e, f)"),
         Line2D([], [], color=PALETTE["violet"], marker="D", lw=0, ms=4.8, mec=PALETTE["ink"],
                label="hierarchical jackknife+, 10 outer folds (e, f)"),
         Line2D([], [], color=PALETTE["ink"], lw=0.8, ls=(0, (4, 3)), label=r"nominal $1-\alpha$ = 0.90"),
         Line2D([], [], color=PALETTE["violet"], lw=0.8, ls=(0, (1, 1.5)),
                label=r"bound of the hierarchical jackknife+, $1-2\alpha$ = 0.80")]
    sub.legend(handles=h, loc="outside lower center", ncol=2, fontsize=6.5, columnspacing=1.2, handletextpad=0.5)


# ------------------------------------------------------------------------------------------ Fig. 8
def _jarvis_fewshot_condensed(fig, D: Dirs, alpha):
    M = pd.read_csv(os.path.join(D.posthoc, "fewshot_by_m.csv"))
    M = M[M.scheme == "one_per_group"]
    W = pd.read_csv(D.f8("fewshot", "worst_groups.csv"))
    mg = sorted(M.m.unique())
    gs = fig.add_gridspec(2, 2)
    h = []
    for c, t in enumerate(("formation", "bandgap")):
        Wt = W[W.target == t].reset_index(drop=True)
        ax = fig.add_subplot(gs[0, c])
        axr = fig.add_subplot(gs[1, c])
        lo, top = 1.0, 1.3
        for j, w in Wt.iterrows():
            col, mk = GROUP_COLORS[j % 3], GROUP_MARKERS[j % 3]
            d = M[(M.target == t) & (M.group == w.group)].sort_values("m")
            dx = (j - 1) * 1.7
            ax.errorbar(d.m + dx, d.cov_mean_split, yerr=[d.cov_mean_split - d.cov_q05_split, d.cov_q95_split - d.cov_mean_split],
                        color=col, marker=mk, ms=3.6, lw=1.0, capsize=0, ls="none")
            ax.axhline(float(w.primary_coverage), color=col, lw=0.9, ls=(0, (1, 1)))
            lo = min(lo, float(w.primary_coverage), float(d.cov_q05_split.min()))
            axr.plot(d.m, d.ratio_offset_to_nomodel, color=col, marker=mk, ms=3.4, lw=1.3)
            axr.plot(d.m, d.ratio_split_to_nomodel, color=col, marker=mk, ms=3.0, mfc="white", lw=0.9, ls=(0, (3, 1.5)))
            top = max(top, float(np.nanmax(d[["ratio_offset_to_nomodel", "ratio_split_to_nomodel"]].to_numpy())) + 0.12)
            h.append(Line2D([], [], color=col, marker=mk, ms=3.6, lw=1.0,
                            label=f"{sg_tex(w.symbol)} ({'a, c' if c == 0 else 'b, d'})"))
        ax.plot(mg, [F81.nominal_split(m, alpha) for m in mg], color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)))
        ax.axhline(0.80, color=PALETTE["red_strong"], lw=0.7, ls=(0, (1, 2)))
        ax.set_ylim(np.floor(lo * 20) / 20 - 0.03, 1.02)
        ax.set_xticks(mg)
        ax.set_xlim(min(mg) - 4, max(mg) + 4)
        ax.set_xlabel("Labelled structures m (one per composition)")
        ax.set_ylabel("Coverage of new structures" if c == 0 else "")
        panel_label(ax, "ab"[c], f"{F8.TNAME[t]}: coverage over draws")
        axr.axhline(1.0, color=PALETTE["grey_dark"], lw=0.8)
        axr.set_ylim(0, top)
        axr.set_xticks(mg)
        axr.set_xlim(min(mg) - 4, max(mg) + 4)
        axr.set_xlabel("Labelled structures m")
        axr.set_ylabel("Width / no-model width" if c == 0 else "")
        panel_label(axr, "cd"[c], f"{F8.TNAME[t]}: width vs no model (post hoc)")
    h += [Line2D([], [], color=PALETTE["grey_dark"], marker="o", ms=3.6, lw=1.0, ls="none",
                 label="split: mean and 5-95 % over 500 draws (a, b)"),
          Line2D([], [], color=PALETTE["grey_dark"], lw=0.9, ls=(0, (1, 1)), label="coverage without labels (a, b)"),
          Line2D([], [], color=PALETTE["grey_mid"], lw=0.8, ls=(0, (4, 3)), label="exact mean k/(m+1) (a, b)"),
          Line2D([], [], color=PALETTE["red_strong"], lw=0.7, ls=(0, (1, 2)), label="0.80 (a, b)"),
          Line2D([], [], color=PALETTE["grey_dark"], marker="o", ms=3.4, lw=1.3, label="offset-corrected / no model (c, d)"),
          Line2D([], [], color=PALETTE["grey_dark"], marker="o", ms=3.0, mfc="white", lw=0.9, ls=(0, (3, 1.5)),
                 label="split / no model (c, d)")]
    fig.legend(handles=h, loc="outside lower center", ncol=3, fontsize=6.5, columnspacing=1.2, handletextpad=0.5)


def fig8_jarvis_fewshot(D: Dirs, out, alpha):
    D.need(os.path.join(D.posthoc, "fewshot_by_m.csv"), D.f8("fewshot", "worst_groups.csv"))
    return single(out, "fig8_jarvis_fewshot", DOUBLE_COL_IN, 5.0, lambda fig: _jarvis_fewshot_condensed(fig, D, alpha))


# ------------------------------------------------------------------------------------------ SI
def figS3_detection_matrix(D: Dirs, out):
    inp = D.inp
    D.need(inp.path("phase4", "report", "phase4_detection_metrics.csv"))
    Mt = inp.csv("phase4", "report", "phase4_detection_metrics.csv")
    settings = [(t, o, lv) for t in ("bandgap", "formation") for o in ("miscovered", "large_error")
                for lv in ("pooled", "within (n-weighted)")]
    cols = [f"{'BG' if t == 'bandgap' else 'FE'}\n{'missed' if o == 'miscovered' else 'large error'}\n"
            f"{'pooled' if lv == 'pooled' else 'within'}" for t, o, lv in settings]
    V = np.full((len(SCORES), len(settings)), np.nan)
    for i, (s, _, _) in enumerate(SCORES):
        for j, (t, o, lv) in enumerate(settings):
            r = Mt[(Mt.target == t) & (Mt.outcome == o) & (Mt.level == lv) & (Mt.score == s)]
            if len(r):
                V[i, j] = float(r.auroc.iloc[0])
    txt = [["" if np.isnan(v) else f"{v:.2f}" for v in row] for row in V]
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    cmap = LinearSegmentedColormap.from_list("auroc", [PALETTE["red_strong"], PALETTE["red_2"], "#FFFFFF", "#B9CDE5",
                                                       PALETTE["blue_main"]])
    norm = TwoSlopeNorm(vmin=0.0, vcenter=0.5, vmax=1.0)

    def draw(fig):
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.025])
        ax, cax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
        P._heatmap(ax, V, [l for _, l, _ in SCORES], cols, cmap, norm, txt, fontsize=6.8)
        ax.xaxis.tick_top()
        ax.tick_params(axis="x", labelsize=6.6)
        ax.axvline(4, color=PALETTE["ink"], lw=1.2)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label("AUROC (0.5 = chance; < 0.5 points the wrong way)", fontsize=6.8)
        cb.outline.set_visible(False)
        panel_label(ax, "", "NOMAD, leave one lattice out: AUROC of each score")
    return single(out, "figS3_detection_auroc_matrix", DOUBLE_COL_IN, 3.0, draw)


def figS4_lattice_level(D: Dirs, out):
    inp = D.inp
    rep = ("phase4", "report")
    D.need(inp.path(*rep, "phase4_lattice_level.csv"), inp.path(*rep, "phase4_lattice_correlations.csv"))
    Lv = inp.csv(*rep, "phase4_lattice_level.csv")
    Co = inp.csv(*rep, "phase4_lattice_correlations.csv")

    def draw(fig):
        ax = fig.add_subplot(111)
        for target, fill in (("bandgap", True), ("formation", False)):
            for _, row in Lv[Lv.target == target].iterrows():
                col = PS.LATTICE_COLORS[row.heldout]
                ax.scatter(row.kde_mean, row.miscoverage, s=34, marker="o" if row.design == "loso" else "^",
                           facecolor=col if fill else "white", edgecolor=col, lw=1.3, zorder=3)
        lines = []
        for target in ("bandgap", "formation"):
            c = Co[(Co.target == target) & (Co.configs == "all 11") & (Co.score == "kde_mean") & (Co.outcome == "miscoverage")]
            if len(c):
                lines.append(f"{'BG' if target == 'bandgap' else 'FE'} $\\rho$ = {c.spearman.iloc[0]:.2f}, "
                             f"p = {c.p_perm.iloc[0]:.2g}")
        ax.axhline(ALPHA, color=PALETTE["ink"], lw=0.9, ls=(0, (4, 3)))
        ax.set_xlabel("Mean KDE dissimilarity of the held-out lattice")
        ax.set_ylabel("Miscoverage of the lattice")
        panel_label(ax, "", "NOMAD, 11 held-out configurations\n" + "; ".join(lines))
        h = [Line2D([], [], marker="o", color=PS.LATTICE_COLORS[l], lw=0, ms=5, label=LATTICE_TEX[l]) for l in LATTICE_ORDER]
        h += [Line2D([], [], marker="o", color=PALETTE["grey_dark"], lw=0, ms=5, label="band gap"),
              Line2D([], [], marker="o", mfc="white", color=PALETTE["grey_dark"], lw=0, ms=5, label="formation"),
              Line2D([], [], marker="^", color=PALETTE["grey_dark"], lw=0, ms=5, label="refit without " + LATTICE_TEX["Fd-3m"]),
              Line2D([], [], color=PALETTE["ink"], lw=0.9, ls=(0, (4, 3)), label=r"$\alpha$ = 0.10")]
        fig.legend(handles=h, loc="outside lower center", ncol=4, fontsize=6.4, columnspacing=1.0, handletextpad=0.4)
    return single(out, "figS4_lattice_level_detection", 4.6, 4.0, draw)


def figS7_hjplus(D: Dirs, out, alpha):
    D.need(*[D.f8("report", f"{t}_hjplus_by_group.csv") for t in ("formation", "bandgap")])

    def draw(fig):
        gs = fig.add_gridspec(1, 2)
        for j, t in enumerate(("formation", "bandgap")):
            G = pd.read_csv(D.f8("report", f"{t}_hjplus_by_group.csv"))
            ax = fig.add_subplot(gs[0, j])
            lo = float(np.nanmin(G[["cov_hjplus", "cov_cvplus", "cov_cvplus_gw"]].to_numpy())) - 0.04
            lo = np.floor(lo * 20) / 20
            ax.plot([lo, 1.01], [lo, 1.01], color=PALETTE["grey_mid"], lw=0.8)
            ax.scatter(G.cov_hjplus, G.cov_cvplus, s=16, facecolor="white", edgecolor=PALETTE["blue_main"], lw=0.9,
                       zorder=3)
            ax.scatter(G.cov_hjplus, G.cov_cvplus_gw, s=16, color=PALETTE["teal"], lw=0, zorder=4)
            for v, ls in ((1 - alpha, (0, (4, 3))), (1 - 2 * alpha, (0, (1, 1.5)))):
                ax.axhline(v, color=PALETTE["ink"], lw=0.7, ls=ls)
                ax.axvline(v, color=PALETTE["ink"], lw=0.7, ls=ls)
            ax.set_xlim(lo, 1.01)
            ax.set_ylim(lo, 1.01)
            ax.set_aspect("equal")
            ax.set_xlabel("Hierarchical jackknife+ coverage")
            ax.set_ylabel("CV+ or group-weighted CV+ coverage" if j == 0 else "")
            panel_label(ax, "ab"[j], f"{F8.TNAME[t]}: {len(G)} held-out groups")
        h = [Line2D([], [], marker="o", color=PALETTE["teal"], lw=0, ms=4.5, label="group-weighted CV+"),
             Line2D([], [], marker="o", mfc="white", color=PALETTE["blue_main"], lw=0, ms=4.5, label="CV+"),
             Line2D([], [], color=PALETTE["grey_mid"], lw=0.8, label="equal coverage"),
             Line2D([], [], color=PALETTE["ink"], lw=0.7, ls=(0, (4, 3)), label=r"$1-\alpha$ = 0.90"),
             Line2D([], [], color=PALETTE["ink"], lw=0.7, ls=(0, (1, 1.5)), label=r"$1-2\alpha$ = 0.80")]
        fig.legend(handles=h, loc="outside lower center", ncol=5, fontsize=6.5, columnspacing=1.0, handletextpad=0.4)
    return single(out, "figS7_hjplus_vs_cvplus", DOUBLE_COL_IN, 3.9, draw)


# ------------------------------------------------------------------------------------------ the set
def jobs(D: Dirs, out, alpha):
    inp = D.inp
    return [
        ("fig1", lambda: fig1_study_design(D, out)),
        ("fig2", lambda: composite(out, "fig2_regimeA_in_distribution", [
            (5.2, old_part(P.fig1_indistribution, inp, out)),
            (3.0, old_part(P.fig2_group_mondrian, inp, out, mapping=dict(a="e", b="f", c="g")))])),
        ("fig3", lambda: composite(out, "fig3_regimeB_unseen_lattice", [(6.0, old_part(P.fig3_unseen_lattice, inp, out))])),
        ("fig4", lambda: composite(out, "fig4_mechanism", [
            (5.9, old_part(P.fig4_mechanism, inp, out)),
            (2.9, old_part(P.fig5_overlap, inp, out, mapping=dict(a="e", b="f")))])),
        ("fig5", lambda: fig5_detection(D, out)),
        ("fig6", lambda: composite(out, "fig6_fewshot_nomad", [
            (5.3, old_part(P.fig7_fewshot, inp, out)),
            (2.6, old_part(P.figS4_dispersion, inp, out, mapping=dict(a="g", b="h")))])),
        ("fig7", lambda: composite(out, "fig7_jarvis_unseen_groups", [
            (5.4, old_part(F8.fig8_jarvis_groups, D.d7, D.d8, out, alpha, leader_lines=True)),
            (2.9, lambda sub: _jarvis_guarantees(sub, D, alpha))], order=[1, 0])),
        ("fig8", lambda: fig8_jarvis_fewshot(D, out, alpha)),
        ("figS1", lambda: composite(out, "figS1_calibration_tail", [(3.1, old_part(P.figS1_calibration_tail, inp, out))])),
        ("figS2", lambda: composite(out, "figS2_transfer_formation", [(3.1, old_part(P.figS2_transfer_formation, inp, out))],
                                    width=4.6)),
        ("figS3", lambda: figS3_detection_matrix(D, out)),
        ("figS4", lambda: figS4_lattice_level(D, out)),
        ("figS5", lambda: composite(out, "figS5_subfamily_and_training",
                                    [(6.6, old_part(P.figS5_subfamily_and_training, inp, out))])),
        ("figS6", lambda: composite(out, "figS6_jarvis_fewshot_per_group",
                                    [(5.25, old_part(F81.fig10_jarvis_fewshot_v81, D.d8, D.posthoc, out, alpha))])),
        ("figS7", lambda: figS7_hjplus(D, out, alpha)),
    ]


def make_figures9(D: Dirs, out: str, only=None, alpha=ALPHA) -> dict:
    PS.apply_publication_style()
    PS.OVERLAP_LOG.clear()
    rec = Recorder()
    done, skipped = {}, {}
    with recording(rec):
        for key, fn in jobs(D, out, alpha):
            if only and not any(key.startswith(o) for o in only):
                continue
            try:
                saved = fn()
                why = None
            except _Skip as e:
                plt.close("all")
                saved, why = None, str(e)
            except Exception as e:      # report and carry on: one bad input must not block the rest
                plt.close("all")
                saved, why = None, f"error: {type(e).__name__}: {e}"
            if saved:
                stem = os.path.splitext(os.path.basename(str(saved[0])))[0]
                issues = PS.OVERLAP_LOG.get(stem, [])
                done[key] = [os.path.basename(str(p)) for p in saved]
                print(f"[figures9] {key}: {stem} | overlap check: {'clean' if not issues else f'{len(issues)} ISSUE(S)'}",
                      flush=True)
                for d in issues:
                    print(f"           {d['kind']}: {d['a']}  <->  {d['b']}", flush=True)
            else:
                skipped[key] = why
                print(f"[figures9] {key}: skipped ({why})", flush=True)
    for p in inp_used(D):
        rec.add(p)
    n_issues = sum(len(v) for v in PS.OVERLAP_LOG.values())
    outside = [p for p in rec.used if not p.startswith(os.path.abspath(D.root) + os.sep)]
    man = dict(created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
               style="figures4papers house style (Chen Liu, CC BY-NC 4.0), adapted in cpaudit/pubstyle.py",
               figures=done, skipped=skipped, inputs_sha256=rec.used, inputs_outside_root=outside,
               files_counted_not_read=D.counted,
               overlap_check=dict(definition="cpaudit/figcheck.py (T1 text/text, T2 text/other axes, L1 legend/data, "
                                             "A1 annotation/data, K1 label/data)",
                                  issues=dict(PS.OVERLAP_LOG), n_issues=n_issues))
    print(f"[figures9] {len(done)} drawn, {len(skipped)} skipped; overlap issues: {n_issues}; "
          f"{len(rec.used)} input files recorded" + (f"; {len(outside)} OUTSIDE the Drive folder" if outside else ""),
          flush=True)
    return man


def inp_used(D: Dirs):
    return list(D.inp.used)


def _dirs_from_cfg(cfg) -> Dirs:
    MOCK = os.environ.get("CPAUDIT_MOCK7", "") == "1"     # same switch as fit7 / run8 / posthoc8
    p = cfg["paths"]
    root = p["drive_root"]
    if not os.path.isdir(root):
        raise FileNotFoundError(f"Drive folder not found: {root}. Mount Drive first.")
    j = lambda k, d: os.path.join(root, p.get(k, d))   # noqa: E731
    d7 = os.path.join(root, "phase07_mock" if MOCK else cfg["models7"]["out_dir"])
    d8 = os.path.join(root, "phase08_mock" if MOCK else cfg["models8"]["out_dir"])
    posthoc = os.path.join(d8, cfg.get("posthoc8", {}).get("out_subdir", "posthoc"))
    return Dirs(root, j("phase1_dir", "phase01_v1"), os.path.join(root, p["phase2_dir"]), j("phase4_dir", "phase04_v1"),
                j("phase5_dir", "phase05_v1"), j("out_dir", "phase06_v1"), d7, d8, posthoc)


def run_figures9(cfg, only=None) -> dict:
    """Full set -> <figures9.figures_dir>/ (refused if that folder exists and is not empty: nothing is ever
    overwritten). With --only, the chosen figures go to a new time-stamped folder under <drafts_dir>/."""
    MOCK = os.environ.get("CPAUDIT_MOCK7", "") == "1"     # same switch as fit7 / run8 / posthoc8
    D = _dirs_from_cfg(cfg)
    F = dict(cfg.get("figures9", {}))
    alpha = float(cfg["models7"]["alpha_primary"])
    full = os.path.join(D.root, "paper_figures_mock9" if MOCK else F.get("figures_dir", "paper_figures_v5"))
    if only:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
        out = os.path.join(D.root, F.get("drafts_dir", "paper_figures_v5_drafts"), stamp + "_" + "-".join(only))
        if os.path.exists(out):
            raise SystemExit(f"[figures9] {out} exists; wait a second and run again.")
    else:
        out = full
        if os.path.isdir(out) and os.listdir(out):
            raise SystemExit(f"[figures9] {out} already exists and is not empty; nothing was overwritten. Rename "
                             "that folder first if the full set should be drawn again (single figures: --only fig1 ...).")
    os.makedirs(out, exist_ok=True)
    print(f"[figures9] reading from {D.root}; writing to {out}", flush=True)
    man = make_figures9(D, out, only, alpha)
    if man["inputs_outside_root"]:
        raise RuntimeError(f"figures read files outside the Drive folder: {man['inputs_outside_root'][:3]}")
    name = "FIGURES9_MANIFEST.json" if not only and not man["skipped"] else "FIGURES9_PARTIAL.json"
    with open(os.path.join(out, name), "w") as fh:
        json.dump(man, fh, indent=2, default=str)
    return man
