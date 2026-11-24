"""Publication figure style for the paper.

Adapted from the house style of figures4papers by Chen Liu
(https://github.com/ChenLiu-1996/figures4papers, scientific-figure-making/references, commit
0ddfbd3, licence CC BY-NC 4.0): minimalist axes (top/right spines off), frameless legends, a
Helvetica/Arial-type sans-serif font, thick axis lines, the semantic blue-green-red-neutral
palette (blue = key/recommended, green = improvement, red = baseline/failure, neutral =
reference), tight_layout finishing and 300-600 dpi raster plus vector export.

Differences from the repository, made on purpose:
* text.usetex is off (Colab has no LaTeX); mathtext gives the few symbols needed.
* Figure sizes are set for a journal page (7.1 in = 180 mm double column) instead of the very
  wide canvases used for slides, so printed text stays >= 7 pt.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

PALETTE = {
    "blue_main": "#0F4D92", "blue_secondary": "#3775BA",
    "green_1": "#DDF3DE", "green_2": "#AADCA9", "green_3": "#8BCF8B",
    "red_1": "#F6CFCB", "red_2": "#E9A6A1", "red_strong": "#B64342",
    "neutral": "#CFCECE", "grey_mid": "#767676", "grey_dark": "#4D4D4D", "ink": "#272727",
    "highlight": "#FFD700", "teal": "#42949E", "violet": "#9A4D8E", "magenta": "#EA84DD",
    "ochre": "#C9A227",   # added: a sixth lattice colour that is not a grey (greys mean 'reference')
}

# one colour per parent lattice; Fd-3m (the failing family) carries the strong red throughout
LATTICE_COLORS = {
    "C2/m": PALETTE["blue_main"], "Pna21": PALETTE["teal"], "R-3c": PALETTE["green_3"],
    "P63/mmc": PALETTE["violet"], "Ia-3": PALETTE["ochre"], "Fd-3m": PALETTE["red_strong"],
}
LATTICE_TEX = {"C2/m": r"$C2/m$", "Pna21": r"$Pna2_1$", "R-3c": r"$R\bar{3}c$", "P63/mmc": r"$P6_3/mmc$",
               "Ia-3": r"$Ia\bar{3}$", "Fd-3m": r"$Fd\bar{3}m$"}
LATTICE_ORDER = ["C2/m", "Pna21", "R-3c", "P63/mmc", "Ia-3", "Fd-3m"]

DOUBLE_COL_IN = 7.1
SINGLE_COL_IN = 3.5


SANS_PREFERENCE = ["Helvetica", "Arial", "Liberation Sans", "FreeSans", "Nimbus Sans", "DejaVu Sans"]


def resolve_sans() -> str:
    """First installed font of the Helvetica/Arial family (Colab ships Liberation Sans)."""
    from matplotlib import font_manager
    have = {f.name for f in font_manager.fontManager.ttflist}
    for name in SANS_PREFERENCE:
        if name in have:
            return name
    return "DejaVu Sans"


def apply_publication_style(font_size: float = 9.0, axes_linewidth: float = 1.2):
    sans = resolve_sans()
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [sans] + [f for f in SANS_PREFERENCE if f != sans],
        # maths in the same face as the text, so space groups (C2/m, R-3c) match the labels
        "mathtext.fontset": "custom",
        "mathtext.rm": sans,
        "mathtext.it": f"{sans}:italic",
        "mathtext.bf": f"{sans}:bold",
        "mathtext.sf": sans,
        "mathtext.cal": sans,
        "font.size": font_size,
        "axes.labelsize": font_size,
        "axes.titlesize": font_size,
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "legend.fontsize": font_size - 1,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": axes_linewidth,
        "xtick.major.width": axes_linewidth,
        "ytick.major.width": axes_linewidth,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": True,
        "savefig.facecolor": "white",
    })


def new_figure(width: float, height: float, w_pad: float = 0.04, h_pad: float = 0.05,
               wspace: float = 0.04, hspace: float = 0.06):
    """Figure with matplotlib's constrained layout, which keeps tick labels, titles, colour bars
    and an 'outside' figure legend from overlapping (more reliable than tight_layout here)."""
    fig = plt.figure(figsize=(width, height), layout="constrained")
    fig.get_layout_engine().set(w_pad=w_pad, h_pad=h_pad, wspace=wspace, hspace=hspace)
    return fig


def panel_label(ax, letter: str, title: str | None = None, title_color: str = "#4D4D4D"):
    """Register a bold panel letter (placed at the panel's top-left, outside its tick labels, when
    the figure is finalised) and an optional left-aligned panel title."""
    ax._cp_panel = letter
    if title:
        ax.set_title(title, loc="left", fontsize=plt.rcParams["font.size"] - 0.5, color=title_color, pad=4)


def _place_panel_letters(fig):
    fig.draw_without_rendering()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for ax in fig.axes:
        letter = getattr(ax, "_cp_panel", None)
        if not letter:
            continue
        (x0, _), (_, y1) = inv.transform(ax.get_tightbbox(r).get_points())
        (ax_x0, _), _ = inv.transform(ax.get_window_extent(r).get_points())
        x = min(x0, ax_x0 - 0.035)
        fig.text(x, y1, letter, fontsize=plt.rcParams["font.size"] + 2.5, fontweight="bold", ha="left",
                 va="top")


def coverage_cmap():
    """Diverging map centred on the nominal level: red = under-coverage, white = nominal,
    blue = over-coverage (conservative)."""
    return LinearSegmentedColormap.from_list(
        "coverage", [PALETTE["red_strong"], PALETTE["red_2"], "#FFFFFF", "#B9CDE5", PALETTE["blue_main"]])


def coverage_norm(nominal: float = 0.9):
    return TwoSlopeNorm(vmin=0.0, vcenter=nominal, vmax=1.0)


def sequential_red():
    return LinearSegmentedColormap.from_list("reds", ["#FFFFFF", PALETTE["red_1"], PALETTE["red_2"],
                                                      PALETTE["red_strong"], "#6E1F1E"])


def sequential_blue():
    return LinearSegmentedColormap.from_list("blues", ["#FFFFFF", "#DCE6F2", "#9DB9DD", PALETTE["blue_secondary"],
                                                       PALETTE["blue_main"]])


def finalize_figure(fig, out_base: str | Path, formats=("png", "pdf", "svg"), dpi: int = 600):
    out_base = Path(out_base)
    os.makedirs(out_base.parent, exist_ok=True)
    _place_panel_letters(fig)
    saved = []
    for fmt in formats:
        p = out_base.with_suffix("." + fmt)
        fig.savefig(p, dpi=dpi if fmt in ("png", "tif", "tiff", "jpg") else None, bbox_inches="tight",
                    pad_inches=0.04)
        saved.append(p)
    plt.close(fig)
    return saved
