"""
crop: cut a rectangular window out of an imported KiCad board (CSXCAD.KiCadImport.KiCadBoard), so a small
RF structure of a large board can be simulated on its own.

    from crop import crop
    crop(board, x0, y0, x1, y1)     # board mm, y up (KiCadImport coordinates); edits the board in place

Copper and other shapes are clipped to the window (a clipped zone with holes is cut into hole-free pieces), drills and footprints outside
are dropped, and the board outline becomes outline AND window. Ports should sit inside the window, away
from its cut edges: the cut copper ends there as an open stub.
"""
import numpy as np
from shapely.geometry import Polygon, box, MultiPolygon, GeometryCollection


def _polys(g):
    if g.is_empty:
        return []
    if isinstance(g, Polygon):
        return [g]
    if isinstance(g, (MultiPolygon, GeometryCollection)):
        return [p for x in g.geoms for p in _polys(x)]
    return []


def _no_holes(p, depth=0):
    """ Hole-free pieces of a polygon: ToCSX draws only the outline of a shape (KiCad fills encode their holes
        as keyhole cuts in that outline), so a clipped zone with interiors is cut along x = hole centroid
        until no piece has a hole. """
    from shapely.ops import split
    from shapely.geometry import LineString
    if not p.interiors or depth > 60:
        return [p]
    x = p.interiors[0].centroid.x
    (_, y0, _, y1) = p.bounds
    parts = _polys(split(p, LineString([(x, y0 - 1), (x, y1 + 1)])))
    if len(parts) < 2:
        return [p]
    return [q for r in parts for q in _no_holes(r, depth + 1)]


def crop(board, x0, y0, x1, y1):
    from CSXCAD.KiCadImport import Shape
    win = box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
    shapes = []
    for s in board.shapes:
        if len(s.pts) < 3:
            continue
        lo, hi = s.pts.min(0), s.pts.max(0)
        if hi[0] < win.bounds[0] or lo[0] > win.bounds[2] or hi[1] < win.bounds[1] or lo[1] > win.bounds[3]:
            continue
        if lo[0] >= win.bounds[0] and hi[0] <= win.bounds[2] and lo[1] >= win.bounds[1] and hi[1] <= win.bounds[3]:
            shapes.append(s)                        # fully inside: unchanged
            continue
        try:
            g = Polygon(s.pts, [h for h in s.holes if len(h) >= 3]).buffer(0).intersection(win)
        except Exception:
            continue
        for p in (q for r in _polys(g) for q in _no_holes(r)):
            if p.area < 1e-6:
                continue
            ext = np.asarray(p.exterior.coords)[:-1]
            shapes.append(Shape(s.layer, ext, s.net, s.kind, s.ref, [], s.ends, s.pad))
    board.shapes = shapes
    board.drills = [d for d in board.drills if win.contains(box(d.x - d.d / 2, d.y - d.d / 2, d.x + d.d / 2, d.y + d.d / 2))]
    if isinstance(board.footprints, dict):
        board.footprints = {k: f for k, f in board.footprints.items() if win.contains(box(f.x - 1e-3, f.y - 1e-3, f.x + 1e-3, f.y + 1e-3))}
    out = board.outline[0] if board.outline else None
    edge = Polygon(out).buffer(0).intersection(win) if out is not None else win
    ext = max(_polys(edge), key=lambda p: p.area) if _polys(edge) else win
    inner = []
    for h in (board.outline[1:] if board.outline else []):
        hp = Polygon(h)
        if hp.within(ext):
            inner.append(np.asarray(h))
    board.outline = [np.asarray(ext.exterior.coords)[:-1]] + inner
    board.cropped = tuple(win.bounds)
    return board


def window_around(board, refs, margin=3.0):
    """ Window (x0, y0, x1, y1) around the pads of the given footprint references, plus a margin (mm). """
    pts = [s.pts for s in board.shapes if s.kind == 'pad' and s.ref in refs]
    if not pts:
        raise ValueError('no pads of %s' % (refs,))
    p = np.vstack(pts)
    lo, hi = p.min(0) - margin, p.max(0) + margin
    return lo[0], lo[1], hi[0], hi[1]
