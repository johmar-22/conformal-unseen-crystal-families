"""Automatic layout check for the paper figures: finds overlapping text and text or legends that
cover plotted data, after the final layout and before the figure is saved.

What counts as an overlap (all measured in display pixels on the laid-out figure):
  T1  two text items intersect (tick labels, axis labels, titles, panel letters, annotations,
      legend entries, heatmap cell labels);
  T2  a text item intersects the plotting area of a different axes;
  L1  a legend box covers plotted data of its own axes (line vertices and segments, markers,
      error bars, filled bands, bars), or a figure-level legend intersects any axes or text;
  A1  an annotation inside an axes covers plotted data of that axes;
  K1  a tick label or axis label covers plotted data or a spine of its own axes (catches markers
      drawn outside the axes, and the radial / angular labels of polar plots).
Lines crossing other lines are not counted (that is the data). Artists marked as background
(`artist.set_gid("bg")`, e.g. a shaded band that a label is meant to sit on) and heatmap cell
labels (`text.set_gid("cell")`, checked only against other text) are exempt where stated.
"""
from __future__ import annotations

import numpy as np
from matplotlib.collections import LineCollection, PathCollection, PolyCollection
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.transforms import Bbox

SHRINK = 0.75  # px; ignore touching boxes (anti-aliasing, glyph padding)


def _shrink(bb: Bbox, s: float = SHRINK) -> Bbox:
    return Bbox.from_extents(bb.x0 + s, bb.y0 + s, bb.x1 - s, bb.y1 - s)


def _hit(a: Bbox, b: Bbox) -> bool:
    a, b = _shrink(a), _shrink(b)
    return a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1


def _poly(t: Text, bb: Bbox, r) -> np.ndarray:
    """Corners (4 x 2, display px) of the text's own rectangle. For rotated text the axis-aligned
    window extent is larger than the glyph box, so the rectangle is rebuilt from the unrotated
    size and the rotation angle (the extent is centred on the rotated rectangle)."""
    ang = np.deg2rad(t.get_rotation() % 360)
    if abs(np.sin(ang)) < 1e-6 and abs(np.cos(ang)) > 0.999:
        return np.array([[bb.x0, bb.y0], [bb.x1, bb.y0], [bb.x1, bb.y1], [bb.x0, bb.y1]])
    # line spacing: public getter only in matplotlib >= 3.11; older versions keep it in _linespacing
    ls = t.get_linespacing() if hasattr(t, "get_linespacing") else getattr(t, "_linespacing", 1.2)
    u = Text(text=t.get_text(), fontproperties=t.get_fontproperties(), linespacing=ls)
    u.set_figure(t.get_figure())
    e = u.get_window_extent(r)
    w, h = e.width, e.height
    c = np.array([(bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2])
    R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    corners = np.array([[-w / 2, -h / 2], [w / 2, -h / 2], [w / 2, h / 2], [-w / 2, h / 2]])
    return c + corners @ R.T


def _shrink_poly(P: np.ndarray, s: float = SHRINK) -> np.ndarray:
    c = P.mean(axis=0)
    d = P - c
    n = np.linalg.norm(d, axis=1, keepdims=True)
    return c + d * np.clip((n - s) / np.maximum(n, 1e-9), 0, None)


def _sat(P: np.ndarray, Q: np.ndarray) -> bool:
    """True if two convex polygons overlap (separating axis theorem), after shrinking both."""
    P, Q = _shrink_poly(P), _shrink_poly(Q)
    for poly in (P, Q):
        for i in range(len(poly)):
            e = poly[(i + 1) % len(poly)] - poly[i]
            nrm = np.array([-e[1], e[0]])
            if not np.any(nrm):
                continue
            a, b = P @ nrm, Q @ nrm
            if a.max() <= b.min() or b.max() <= a.min():
                return False
    return True


def _rect(bb: Bbox) -> np.ndarray:
    return np.array([[bb.x0, bb.y0], [bb.x1, bb.y0], [bb.x1, bb.y1], [bb.x0, bb.y1]])


def _visible_text(t: Text) -> bool:
    return t.get_visible() and bool(str(t.get_text()).strip()) and t.get_alpha() != 0


def _label(t: Text) -> str:
    s = str(t.get_text()).replace("\n", " ")
    return s if len(s) <= 40 else s[:37] + "..."


def _texts(fig, r):
    """(bbox, kind, axes or None, Text) for every visible text in the figure."""
    out = []
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        for t in (ax.title, ax._left_title, ax._right_title):
            if _visible_text(t):
                out.append((t.get_window_extent(r), "title", ax, t))
        for axis, kind in ((ax.xaxis, "xlabel"), (ax.yaxis, "ylabel")):
            if axis.get_visible() and _visible_text(axis.label):
                out.append((axis.label.get_window_extent(r), kind, ax, axis.label))
            if not axis.get_visible():
                continue
            for tick in axis._update_ticks():   # only the ticks that are actually drawn
                for lab in (tick.label1, tick.label2):
                    if lab.get_visible() and _visible_text(lab):
                        out.append((lab.get_window_extent(r), "tick", ax, lab))
            ot = axis.get_offset_text()
            if _visible_text(ot):
                out.append((ot.get_window_extent(r), "tick", ax, ot))
        for t in ax.texts:
            if _visible_text(t):  # Text.get_window_extent: the text box only, not an annotation arrow
                out.append((Text.get_window_extent(t, r), "cell" if t.get_gid() == "cell" else "annotation", ax, t))
        leg = ax.get_legend()
        if leg is not None and leg.get_visible():
            for t in leg.get_texts():
                if _visible_text(t):
                    out.append((t.get_window_extent(r), "legend", ax, t))
    for t in fig.texts:
        if _visible_text(t):
            out.append((t.get_window_extent(r), "figtext", None, t))
    for leg in fig.legends:
        for t in leg.get_texts():
            if _visible_text(t):
                out.append((t.get_window_extent(r), "figlegend", None, t))
    return out


def _densify(xy: np.ndarray, step: float = 1.0) -> np.ndarray:
    """Points every `step` px along a polyline given in display coordinates (NaN breaks it)."""
    if len(xy) < 2:
        return xy
    pts = [xy[:1]]
    for a, b in zip(xy[:-1], xy[1:]):
        if not (np.isfinite(a).all() and np.isfinite(b).all()):
            continue
        n = int(np.ceil(np.hypot(*(b - a)) / step))
        if n > 1:
            t = np.linspace(0, 1, n + 1)[1:, None]
            pts.append(a + t * (b - a))
        else:
            pts.append(b[None])
    out = np.vstack(pts)
    return out[np.isfinite(out).all(axis=1)]


class _Adder:
    def __init__(self, fn):
        self.append = fn


def _data_shapes(ax, r):
    """Display-space geometry of plotted data in `ax`: list of (kind, object, clip) where kind is
    'pts' (N x 2 points with a radius), 'path' (filled matplotlib Path) or 'bbox', and clip is the
    axes box for artists clipped to the axes (only the visible part can overlap anything)."""
    raw = []
    ab = ax.get_window_extent(r)

    def add(kind, obj, artist):
        raw.append((kind, obj, ab if artist.get_clip_on() else None))

    shapes = _Adder(add)
    for ln in ax.get_lines():
        if not ln.get_visible() or ln.get_gid() == "bg":
            continue
        xy = ln.get_transform().transform(np.asarray(ln.get_xydata(), float))
        if ln.get_linestyle() not in ("None", "none", "", " ") and ln.get_linewidth() > 0:
            shapes.append("pts", (_densify(xy), max(ln.get_linewidth(), 1.0) / 2 * r.points_to_pixels(1)), ln)
        if ln.get_marker() not in (None, "None", "none", "", " "):
            rad = ln.get_markersize() / 2 * r.points_to_pixels(1)
            shapes.append("pts", (xy[np.isfinite(xy).all(axis=1)], rad), ln)
    for c in ax.collections:
        if not c.get_visible() or c.get_gid() == "bg":
            continue
        if isinstance(c, PathCollection):
            off = c.get_offset_transform().transform(np.asarray(c.get_offsets(), float))
            sz = np.asarray(c.get_sizes(), float)
            rad = (np.sqrt(sz.max()) if sz.size else 6.0) / 2 * r.points_to_pixels(1)
            shapes.append("pts", (off[np.isfinite(off).all(axis=1)], rad), c)
        elif isinstance(c, LineCollection):
            tr = c.get_transform()
            lw = (np.max(c.get_linewidths()) if len(c.get_linewidths()) else 1.0) / 2 * r.points_to_pixels(1)
            for seg in c.get_segments():
                shapes.append("pts", (_densify(tr.transform(np.asarray(seg, float))), max(lw, 0.5)), c)
        elif isinstance(c, PolyCollection):
            tr = c.get_transform()
            for p in c.get_paths():
                shapes.append("path", p.transformed(tr), c)
        else:  # other collections: use their bounding box
            try:
                shapes.append("bbox", c.get_window_extent(r), c)
            except Exception:
                pass
    for p in ax.patches:
        if not p.get_visible() or p.get_gid() == "bg" or p is ax.patch:
            continue
        shapes.append("bbox", p.get_window_extent(r), p)
    return raw


def _covers(bb: Bbox, shapes) -> bool:
    for kind, obj, clip in shapes:
        box = bb
        if clip is not None:
            box = Bbox.intersection(bb, clip)
            if box is None or box.width <= 2 * SHRINK or box.height <= 2 * SHRINK:
                continue
        b = _shrink(box)
        if kind == "pts":
            xy, rad = obj
            if xy.size == 0:
                continue
            inside = ((xy[:, 0] + rad > b.x0) & (xy[:, 0] - rad < b.x1) &
                      (xy[:, 1] + rad > b.y0) & (xy[:, 1] - rad < b.y1))
            if inside.any():
                return True
        elif kind == "path":
            if obj.intersects_bbox(b, filled=True):
                return True
        elif kind == "bbox":
            if _hit(box, obj):
                return True
    return False


def check_overlaps(fig) -> list[dict]:
    """Return a list of overlap issues for a laid-out figure (empty list = clean)."""
    fig.draw_without_rendering()
    r = fig.canvas.get_renderer()
    issues = []
    T = _texts(fig, r)
    polys = [_poly(t, bb, r) for bb, k, a, t in T]
    # T1: text vs text
    for i in range(len(T)):
        for j in range(i + 1, len(T)):
            (b1, k1, a1, t1), (b2, k2, a2, t2) = T[i], T[j]
            if t1 is t2:
                continue
            if k1 in ("legend", "figlegend") and k1 == k2:
                continue  # the legend lays out its own entries
            if _hit(b1, b2) and _sat(polys[i], polys[j]):
                issues.append(dict(kind="T1 text/text", a=f"{k1}: {_label(t1)}", b=f"{k2}: {_label(t2)}"))
    # T2: text vs another axes' plotting area
    areas = [(ax, ax.get_window_extent(r)) for ax in fig.axes if ax.get_visible()]
    for (bb, k, a, t), P in zip(T, polys):
        for ax, ab in areas:
            if ax is a or (a is not None and (ax in getattr(a, "child_axes", []) or a in getattr(ax, "child_axes", []))):
                continue
            if k == "legend" and ax is a:
                continue
            if k == "figtext":   # panel letters must also keep a visible gap (4 pt) from other panels
                g = 4 * r.points_to_pixels(1)
                bbx = Bbox.from_extents(bb.x0 - g, bb.y0 - g, bb.x1 + g, bb.y1 + g)
                if _hit(bbx, ab):
                    issues.append(dict(kind="T2 panel letter closer than 4 pt to another panel", a=f"{k}: {_label(t)}",
                                       b=f"axes '{ax.get_title(loc='left')[:30]}'"))
                continue
            if _hit(bb, ab) and _sat(P, _rect(ab)):
                issues.append(dict(kind="T2 text/other axes", a=f"{k}: {_label(t)}", b=f"axes '{ax.get_title(loc='left')[:30]}'"))
    # L1 / A1: legends and annotations vs data
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        shapes = _data_shapes(ax, r)
        leg = ax.get_legend()
        if leg is not None and leg.get_visible() and _covers(leg.get_window_extent(r), shapes):
            issues.append(dict(kind="L1 legend/data", a="legend", b=f"axes '{ax.get_title(loc='left')[:30]}'"))
        for t in ax.texts:
            if _visible_text(t) and t.get_gid() != "cell" and _covers(Text.get_window_extent(t, r), shapes):
                issues.append(dict(kind="A1 annotation/data", a=_label(t), b=f"axes '{ax.get_title(loc='left')[:30]}'"))
    # K1: tick and axis labels vs the data and spines of their own axes
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        shapes = _data_shapes(ax, r)
        for sp in ax.spines.values():
            if sp.get_visible():
                path = sp.get_path().transformed(sp.get_transform())
                shapes.append(("pts", (_densify(np.asarray(path.vertices, float)), max(sp.get_linewidth(), 1.0) / 2), None))
        for bb, k, a, t in T:
            if a is ax and k in ("tick", "xlabel", "ylabel") and _covers(bb, shapes):
                issues.append(dict(kind="K1 label/data", a=f"{k}: {_label(t)}", b=f"axes '{ax.get_title(loc='left')[:30]}'"))
    for leg in fig.legends:
        lb = leg.get_window_extent(r)
        for ax, ab in areas:
            if _hit(lb, ab):
                issues.append(dict(kind="L1 figure legend/axes", a="figure legend", b=f"axes '{ax.get_title(loc='left')[:30]}'"))
        for bb, k, a, t in T:
            if k != "figlegend" and _hit(lb, bb):
                issues.append(dict(kind="L1 figure legend/text", a="figure legend", b=f"{k}: {_label(t)}"))
    # de-duplicate
    seen, out = set(), []
    for d in issues:
        key = (d["kind"], d["a"], d["b"])
        if key not in seen:
            seen.add(key)
            out.append(d)
    return out
