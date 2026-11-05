#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2026 Zoltan Doczi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""
kicad_sim: S-parameters of a KiCad board with openEMS, in one step.

    python kicad_sim.py board.kicad_pcb --port J1 --port J2 --fmax 6e9 [--fab JLCPCB:JLC04161H-7628]

Ports are lumped 50 Ohm ports at the outer edge of a pad (given by footprint reference, the
pad on an outer copper layer that is not on a ground net). The port runs vertically from the
sheet of the nearest copper layer under the pad that has a ground plane there to the pad
sheet. The mesh is built from the stackup (every dielectric at least four cells), the copper
of the port nets (thirds rule on axis-aligned track edges) and the port planes. Every port is
excited once, so the result is the full N-port S matrix; it is written as Touchstone.
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

C0 = 299792458.0
HERE = os.path.dirname(os.path.abspath(__file__))
def _find_engine():
    """ The openEMS executable: Z2_OPENEMS, else the engine bundled with a Z2 EMS release package
        (<package>/bin/openEMS[.exe] next to twin/), else openEMS on PATH. """
    name = 'openEMS.exe' if os.name == 'nt' else 'openEMS'
    for c in (os.environ.get('Z2_OPENEMS'), os.path.join(HERE, '..', 'bin', name), shutil.which('openEMS')):
        if c and os.path.isfile(c):
            return os.path.abspath(c)
    return name


ENGINE = _find_engine()
ENGINE_DIR = os.path.dirname(ENGINE) or '.'
_CUDA_EXE = os.environ.get('Z2_OPENEMS_CUDA') or ENGINE
VARIANTS = {
    # CPU engine (multithreaded, SSE); the executable comes from _find_engine()
    'cpu': dict(exe=ENGINE, path=[ENGINE_DIR]),
    # GPU engine (--engine=cuda) of a -DWITH_CUDA=ON build: the release packages ship one executable with both
    # engines; Z2_OPENEMS_CUDA=<executable> and Z2_OPENEMS_CUDA_LIBS=<library folders, os.pathsep separated>
    # point to a separate GPU build
    'cuda': dict(exe=_CUDA_EXE, path=[os.path.dirname(_CUDA_EXE) or '.'],
                 libs=[d for d in os.environ.get('Z2_OPENEMS_CUDA_LIBS', '').split(os.pathsep) if d],
                 engine=['--engine=cuda']),
}
VARIANTS['local'] = VARIANTS['cpu']          # older name of the CPU variant
DEFAULT_VARIANT = 'cpu'
# material options passed to KiCadBoard.ToCSX (see CSXCAD.PCBMaterials)
MATERIALS = dict(dispersion=None, conductor='pec', roughness_um=None, components=False)
# PML: Mur does not absorb the low-frequency/static part of the field, the energy then never reaches the
# end criterion and every run goes to MaxTime (coupon_msl50: 87438 vs 15860 timesteps, same S within 0.005 dB)
BOUNDARY = ['PML_8'] * 6


def pml_pad(lines, lo, hi, n=10):
    """ At least n mesh cells of air outside [lo, hi] on both sides (the PML takes 8 of them): the domain
        grows outward with the outermost cell size. """
    L = sorted(float(v) for v in lines)
    for side in (0, 1):
        out = [v for v in L if (v < lo - 1e-9 if side == 0 else v > hi + 1e-9)]
        if len(out) >= n:
            continue
        step = (L[1] - L[0]) if side == 0 else (L[-1] - L[-2])
        for _ in range(n - len(out)):
            L = [L[0] - step] + L if side == 0 else L + [L[-1] + step]
    return np.array(L)
GROUND = re.compile(r'^(GND|AGND|DGND|PGND|GNDA|GNDD|VSS|0V|EARTH|CHASSIS)', re.I)


def _setup_python():
    if os.name == 'nt':
        os.environ.setdefault('CSXCAD_INSTALL_PATH', ENGINE_DIR)
        os.environ.setdefault('OPENEMS_INSTALL_PATH', ENGINE_DIR)
        os.environ['PATH'] = ENGINE_DIR + os.pathsep + os.environ.get('PATH', '')
    import CSXCAD
    fork = os.environ.get('Z2_CSXCAD_PY', os.path.join(HERE, '..', 'CSXCAD', 'python', 'CSXCAD'))
    if os.path.isdir(fork) and fork not in CSXCAD.__path__:
        CSXCAD.__path__.insert(0, fork)       # KiCadImport / FabStackups of the fork


# ---------------------------------------------------------------- ports
def find_pad(board, ref):
    pads = [s for s in board.shapes if s.kind == 'pad' and s.ref == ref and s.layer in (board.copper_layers[0], board.copper_layers[-1])]
    if not pads:
        raise ValueError('no outer-layer pad of %s' % ref)
    sig = [p for p in pads if not GROUND.match(p.net or '')]
    return (sig or pads)[0]


def _inside(board, layer, x, y, exclude_net):
    from CSXCAD.KiCadImport import PointsInPolygon
    for s in board.ShapesOn(layer):
        if s.net == exclude_net:
            continue
        if PointsInPolygon([(x, y)], s.pts)[0]:
            return s.net
    return None


def make_port(board, pad, number, prefer='edge', side=None):
    # Z2_PORT_FORCE=<json file>: use the port dicts stored there (list, by number) verbatim -- the same port
    # definition for two models of one board (e.g. PADS-import vs gerber path)
    pf = os.environ.get('Z2_PORT_FORCE')
    if pf:
        import json as _j
        for q in _j.load(open(pf)):
            if q.get('number') == number:
                q = dict(q); q.pop('moved_from_pad', None)
                return q
    """ Lumped port at the pad edge that faces the nearest board edge (prefer='edge'), or at the opposite
        edge, the end that faces the line (prefer='inner': edge connectors whose planes are pulled back
        under the board-edge end of the pad, or pads that touch the board outline). """
    # Z2_PORT_TRACK_EXIT=1: prefer the port across the track that leaves the pad (a THT pin pad has copper on
    # both outer layers: a port sheet at the pad edge would short against the pad's own back-side copper)
    if os.environ.get('Z2_PORT_TRACK_EXIT') == '1' and getattr(pad, 'pad', None) != 'TRK':
        moved = port_at_track_exit(board, pad, number)
        if moved is not None:
            return moved
    (bx0, bx1), (by0, by1) = board.BoundingBox()[:2]
    P = pad.pts
    px0, py0 = P.min(0)
    px1, py1 = P.max(0)
    cx, cy = 0.5 * (px0 + px1), 0.5 * (py0 + py1)
    d = {'-x': cx - bx0, '+x': bx1 - cx, '-y': cy - by0, '+y': by1 - cy}
    # the port plane sits slightly inside the pad: a plane exactly on the pad and plane edges can miss
    # both conductors (edge points are not reliably inside a polygon) and leave the port floating
    inset = min(0.1, 0.25 * min(px1 - px0, py1 - py0))
    qx0, qx1, qy0, qy1 = px0 + inset, px1 - inset, py0 + inset, py1 - inset
    layers = board.copper_layers
    i = layers.index(pad.layer)
    order = layers[i + 1:] if i == 0 else list(reversed(layers[:i]))
    # preferred: the pad edge facing the nearest board edge; an edge connector often has the planes pulled
    # back (or cut out) under that end, then the other pad edges are tried, nearest board edge first
    ref_layer = None
    sides = sorted(d, key=d.get)
    if side in d:                                  # forced (e.g. the side where the track leaves the pad)
        opp = {'-x': '+x', '+x': '-x', '-y': '+y', '+y': '-y'}[side]
        sides = [side, opp] + [x for x in sides if x not in (side, opp)]
    elif prefer == 'inner':
        opp = {'-x': '+x', '+x': '-x', '-y': '+y', '+y': '-y'}[sides[0]]
        sides = [opp] + [x for x in sides if x != opp]
    else:
        opp = {'-x': '+x', '+x': '-x', '-y': '+y', '+y': '-y'}[sides[0]]
        sides = [sides[0], opp] + [x for x in sides[1:] if x != opp]   # then the same axis, so the plane stays across the line
    for side in sides:
        if side[1] == 'x':
            xp = qx0 if side == '-x' else qx1
            a, b = (xp, qy0 - inset), (xp, qy1 + inset)
            probe = (xp, cy)
            axis = 0
        else:
            yp = qy0 if side == '-y' else qy1
            a, b = (qx0 - inset, yp), (qx1 + inset, yp)
            probe = (cx, yp)
            axis = 1
        for l in order:
            n = _inside(board, l, probe[0], probe[1], pad.net)
            if n is not None:
                ref_layer, ref_net = l, n
                break
        if ref_layer is not None:
            break
    if ref_layer is None:
        if getattr(pad, 'pad', None) != 'TRK':       # not already a synthetic track pad
            moved = port_at_track_exit(board, pad, number)
            if moved is not None:
                return moved
        raise ValueError('port %d (%s): no reference copper under the pad' % (number, pad.ref))
    zs, zr = board.SheetZ(pad.layer), board.SheetZ(ref_layer)
    return dict(number=number, ref=pad.ref, net=pad.net, layer=pad.layer, ref_layer=ref_layer, ref_net=ref_net,
                start=[a[0], a[1], zr], stop=[b[0], b[1], zs], side=side, normal_axis=axis,
                width=float(np.hypot(b[0] - a[0], b[1] - a[1])))


def port_at_track_exit(board, pad, number, step=0.05, margin=0.15):
    """ Port for a pad with no reference copper under it (through-hole centre pin of a vertical
        SMA: every plane has an antipad there). Walk out along the track that leaves the pad until
        a reference plane is under the track, and put the port there, across the track: that is the
        start of the line over its reference, the plane a VNA de-embeds to. Returns None when no
        track of the pad's net leaves it or no reference copper turns up along the first segment. """
    from shapely.geometry import Polygon, Point
    from CSXCAD.KiCadImport import Shape
    P = Polygon(pad.pts).buffer(0)
    layers = board.copper_layers
    i = layers.index(pad.layer)
    order = layers[i + 1:] if i == 0 else list(reversed(layers[:i]))
    best = None
    for s in board.shapes:
        if s.net != pad.net or s.layer != pad.layer or s.kind not in ('track', 'arc') or s.ends is None:
            continue
        a, b = np.array(s.ends[0], float), np.array(s.ends[1], float)
        for p_in, p_out in ((a, b), (b, a)):
            if not P.buffer(0.05).contains(Point(*p_in)):
                continue
            L = float(np.hypot(*(p_out - p_in)))
            if L < 1e-6:
                continue
            u = (p_out - p_in) / L
            width = float(Polygon(s.pts).buffer(0).area / max(L, 1e-6)) if len(s.pts) > 2 else 0.2
            t = 0.0
            while t <= L:
                q = p_in + u * t
                if not P.contains(Point(*q)):
                    found = next((l for l in order if _inside(board, l, q[0] + u[0] * margin, q[1] + u[1] * margin, pad.net)), None)
                    if found:
                        q = q + u * margin
                        if best is None or t < best[0]:
                            best = (t, q, u, max(width, 0.1))
                        break
                t += step
    if best is None:
        return None
    _, q, u, w = best
    # axis-aligned stand-in pad across the track (lumped ports are axis-aligned boxes)
    horiz = abs(u[0]) >= abs(u[1])
    h, ln = w / 2.0, 0.1
    if horiz:
        pts = np.array([[q[0] - ln, q[1] - h], [q[0] + ln, q[1] - h], [q[0] + ln, q[1] + h], [q[0] - ln, q[1] + h]])
        side = '+x' if u[0] > 0 else '-x'
    else:
        pts = np.array([[q[0] - h, q[1] - ln], [q[0] + h, q[1] - ln], [q[0] + h, q[1] + ln], [q[0] - h, q[1] + ln]])
        side = '+y' if u[1] > 0 else '-y'
    tp = Shape(pad.layer, pts, pad.net, 'pad', pad.ref, [], None, 'TRK')
    port = make_port(board, tp, number, prefer='inner', side=side)
    port['moved_from_pad'] = round(float(best[0]), 3)
    return port


# ---------------------------------------------------------------- mesh
def smooth(lines, res, ratio=1.4, merge=None):
    from CSXCAD.SmoothMeshLines import SmoothMeshLines
    v = np.unique(np.round(np.asarray(lines, float), 6))
    merge = res / 20 if merge is None else merge
    keep = [v[0]]
    for x in v[1:]:
        if x - keep[-1] >= merge:
            keep.append(x)
    if keep[-1] != v[-1]:
        keep[-1] = v[-1]
    return SmoothMeshLines(np.array(keep), res, ratio)


def grade_ratio():
    """ the neighbour-cell ratio limit of grade(): Z2_GRADE_RATIO (default 1.5); Z2_GRADE=0 turns grading off """
    if os.environ.get('Z2_GRADE', '1') in ('0', 'off', 'no'):
        return None
    return float(os.environ.get('Z2_GRADE_RATIO', '1.5'))


def _sig7(v):
    """ v as Write2XML stores it (7 significant digits): mesh lines and the primitives that must lie on them """
    return float('%.6e' % float(v))


def _force_line(lines, v):
    """ Put v on the mesh: move the nearest line onto it when it is closer than a third of the local cell,
        else insert it. Edits the list in place. """
    a = np.asarray(lines, float)
    i = int(np.argmin(np.abs(a - v)))
    if abs(a[i] - v) < 1e-6:
        return
    nb = [abs(a[j] - a[i]) for j in (i - 1, i + 1) if 0 <= j < len(a)]
    if nb and abs(a[i] - v) < min(nb) / 3:
        lines[i] = v
    else:
        lines.append(v)
    lines.sort()


def grade(lines, ratio=1.5, max_iter=20000):
    """ Limit the size ratio of neighbouring cells: SmoothMeshLines only fills gaps larger than the
        resolution, so a 0.05 mm cell next to a 0.27 mm one survives it (ratio 5). Such jumps cost
        accuracy; here the larger cell of the worst pair is split, growing geometrically from the
        smaller one (or halved when the remainder would be smaller than the small cell), until every
        ratio is at most `ratio`. The lines given are kept.

        No cell is made smaller than the smallest cell of the input: halving below it cascaded outward
        on the 8L array (25 um -> 13 um cells, 564 MCells, half the timestep). With that floor the array
        keeps 25 um (512 MCells, timestep x2) and the coupon_msl50 check stays: Z0 50.28 vs 50.25 ohm,
        eeff 3.360 vs 3.358, 283 s vs 382 s. (Leaving every pair below 2x unsplit was tried and is wrong:
        the strip's transverse lines thinned out, Z0 48.1 ohm, eeff 4.40.) """
    L = list(np.unique(np.asarray(lines, float)))
    floor = float(np.diff(L).min()) * 0.999 if len(L) > 1 else 0.0
    for _ in range(max_iter):
        d = np.diff(L)
        if len(d) < 2:
            break
        r = np.maximum(d[1:] / d[:-1], d[:-1] / d[1:])
        big, small = np.maximum(d[1:], d[:-1]), np.minimum(d[1:], d[:-1])
        stuck = (big < small * (1 + ratio)) & (big / 2 < floor)     # only a halving below the floor would help
        r[stuck] = np.minimum(r[stuck], ratio)
        i = int(np.argmax(r))
        if r[i] <= ratio * 1.001:
            break
        if d[i + 1] > d[i]:                  # the big cell is on the right of line i+1: grow into it
            s, b, x0, sign = d[i], d[i + 1], L[i + 1], 1.0
        else:                                # the big cell is on the left of line i+1
            s, b, x0, sign = d[i + 1], d[i], L[i + 1], -1.0
        step = s * ratio if b - s * ratio >= s else b / 2
        L.insert(i + 1 + (1 if sign > 0 else 0), x0 + sign * step)
    return np.array(L)


def mesh_lines(board, ports, fmax, res_div=20, margin=None, coarse=1.0, copper='sheet', refine=1.0):
    er = max([l.epsilon_r or 1.0 for l in board.stackup if l.type in ('prepreg', 'core')] + [1.0])
    res = C0 / (fmax * math.sqrt(er)) / 1e-3 / res_div * coarse
    # Z2_AIR_MARGIN=<mm>: air around the board for radiating structures (antennas need >= lambda/4 at fmin)
    if margin is None and os.environ.get('Z2_AIR_MARGIN'):
        margin = float(os.environ['Z2_AIR_MARGIN'])
    air = margin if margin is not None else max(4.0, C0 / fmax / 1e-3 / 8)
    (bx0, bx1), (by0, by1) = board.BoundingBox()[:2]
    zt = board.SheetZ(board.copper_layers[0])
    zb = board.SheetZ(board.copper_layers[-1])
    X = [bx0 - air, bx0, bx1, bx1 + air]
    Y = [by0 - air, by0, by1, by1 + air]
    Z = [zb - air, zt + air]
    nets = {p['net'] for p in ports}
    pnets = set(nets)                                   # the port nets (signal vias)
    # small boards (coupons, calibration cards): every net gets the conductor rules, not only the port nets --
    # a gap-coupled resonator is on its own net (the ring of the calibration card was meshed only by the coarse
    # grid and broke into 36 islands at fmax 6 GHz, 2026-10-21). Z2_MESH_NETS=ports|all overrides.
    mode = os.environ.get('Z2_MESH_NETS', 'auto')
    if mode == 'all' or (mode == 'auto' and len(board.shapes) <= 3000):
        refs_ = {p.get('ref_net') for p in ports} | {'GND'}          # not the ground net: its pads/fence would cost 3x cells
        nets = ({s.net for s in board.shapes if s.kind != 'zone'} - refs_) | nets
    fine = res

    # refine < 1 (the convergence check) also refines the conductor rule: the cells at the strip edges,
    # not only the largest cell, decide the line impedance and delay
    cf = min(1.0, refine)
    # lines across a strip / pad: 4 -> 8 moved the 50 ohm coupon's Z0 from 50.9 to 49.4 ohm (reference 49.0) with
    # conformal PEC, at 2.4x the cells of a small board; Z2_N_ACROSS for mesh studies
    n_base = int(os.environ.get('Z2_N_ACROSS', '8'))
    n_across = max(4, int(round(n_base / cf)))

    def thirds(lo, hi, r):
        # thirds rule: no line on a metal edge, one r/3 inside and one 2r/3 outside, plus a few
        # lines across the conductor; a line exactly on the edge widens the strip by half a cell
        r = r * cf
        return [lo + r / 3, lo - 2 * r / 3, hi - r / 3, hi + 2 * r / 3] + list(np.linspace(lo + r / 3, hi - r / 3, n_across))

    for s in board.shapes:
        if s.net not in nets or s.kind == 'zone':
            continue
        x0, y0 = s.pts.min(0)
        x1, y1 = s.pts.max(0)
        w, h = x1 - x0, y1 - y0
        from CSXCAD.KiCadImport import PolygonArea
        rect = abs(PolygonArea(s.pts)) > 0.97 * w * h        # axis-aligned rectangle (pad)
        if s.kind == 'track' and s.ends is not None and (abs(s.ends[0][1] - s.ends[1][1]) < 1e-6 or abs(s.ends[0][0] - s.ends[1][0]) < 1e-6):
            horiz = abs(s.ends[0][1] - s.ends[1][1]) < 1e-6
            r = (h if horiz else w) / 3
            fine = min(fine, r * cf)
            if horiz:
                Y += thirds(y0, y1, r)
                X += [min(s.ends[0][0], s.ends[1][0]), max(s.ends[0][0], s.ends[1][0])]
            else:
                X += thirds(x0, x1, r)
                Y += [min(s.ends[0][1], s.ends[1][1]), max(s.ends[0][1], s.ends[1][1])]
        elif rect:
            r = max(min(w, h) / 3, 0.025)                    # floor: sliver pads (fmcw3 coupler, 34 um) gave 11 um cells
            fine = min(fine, r * cf)
            X += thirds(x0, x1, min(r, w / 3))
            Y += thirds(y0, y1, min(r, h / 3))
        elif s.kind in ('track', 'arc') and s.ends is not None:
            # diagonal track / arc segment: the staircase keeps a thin strip connected only when the cells across it
            # are <= w/3 in BOTH axes (a 0.37 mm 10-degree line meshed at 0.7 mm was an open circuit, 2026-10-21),
            # so the strip's bounding box gets a uniform fine grid
            L_ = math.hypot(s.ends[1][0] - s.ends[0][0], s.ends[1][1] - s.ends[0][1])
            wtr = abs(PolygonArea(s.pts)) / L_ if L_ > 1e-6 else min(w, h)
            r = max(min(wtr / 3, res), 0.03)
            fine = min(fine, r * cf)
            X += list(np.arange(x0 - r, x1 + r + 1e-9, r))
            Y += list(np.arange(y0 - r, y1 + r + 1e-9, r))
        else:
            X += [x0, x1]
            Y += [y0, y1]
    # copper drawn in the port footprints themselves (printed couplers, filters: net-less graphics): lines by
    # the thirds rule at their long axis-aligned edges, so narrow gaps between such conductors stay open
    prefs = {p.get('ref') for p in ports} - {None}
    rg = 0.1
    for s in board.shapes:
        if s.kind != 'graphic' or not s.layer.endswith('.Cu') or len(s.pts) < 3:
            continue
        # also net-level copper polygons of the meshed nets (printed antennas drawn as board polygons)
        if s.ref not in prefs and not (mode != 'ports' and s.net in nets):
            continue
        from shapely.geometry import Polygon as _P, Point as _Pt
        P = _P(s.pts).buffer(0)
        q = np.vstack([s.pts, s.pts[:1]])
        for a, b in zip(q[:-1], q[1:]):
            for ax in (0, 1):
                o = 1 - ax
                if abs(a[o] - b[o]) < 1e-6 and abs(a[ax] - b[ax]) > 0.5:      # edge along ax at coordinate a[o]
                    c = a[o]
                    m = 0.5 * (a[ax] + b[ax])
                    probe = [0.0, 0.0]; probe[ax] = m; probe[o] = c + 1e-3
                    sgn = 1 if P.contains(_Pt(*probe)) else -1        # side of the edge the copper is on
                    (Y if ax == 0 else X).extend([c + sgn * rg / 3, c - sgn * 2 * rg / 3])
                    fine = min(fine, rg / 3)
        # curved / slanted thin copper (printed IFA ring arm, meander drawn as one polygon): the staircase keeps it
        # connected only with cells <= w/3 in BOTH axes; the axis-edge lines above do nothing for an arc, the
        # WES0115 ring arm (1 mm wide) broke into islands at 0.3-2.4 mm cells (2026-11-04) -> uniform grid over
        # the bbox at w_eff/3, w_eff = 2*area/perimeter
        if P.area > 0 and P.length > 0:
            w_eff = 2 * P.area / P.length
            bx0, by0, bx1, by1 = P.bounds
            if P.area < 0.6 * (bx1 - bx0) * (by1 - by0) and 0.15 < w_eff < 3.0:
                rr = max(min(w_eff / 3, res), 0.05)
                fine = min(fine, rr * cf)
                X.extend(np.arange(bx0 - rr, bx1 + rr + 1e-9, rr))
                Y.extend(np.arange(by0 - rr, by1 + rr + 1e-9, rr))
                print('  thin polygon %s w_eff %.2f mm -> uniform %.3f mm grid' % (s.net, w_eff, rr))
    # explicit refinement boxes: Z2_REFINE_BOX="x0,y0,x1,y1,r;..." (board mm) -> uniform lines every r. For copper
    # the net rules do not see: the WES0115 IFA ring arm is one polygon with the GND disc (shorted by the IFA stub)
    # and broke into islands on the coarse ground-net grid (2026-11-04)
    for bx in filter(None, os.environ.get('Z2_REFINE_BOX', '').split(';')):
        bx0, by0, bx1, by1, rr = [float(v) for v in bx.split(',')]
        fine = min(fine, rr * cf)
        X.extend(np.arange(bx0, bx1 + 1e-9, rr))
        Y.extend(np.arange(by0, by1 + 1e-9, rr))
        print('  refine box %s -> %.3f mm grid' % (bx, rr))
    # a mesh line through every plated hole of a reference net, so the stitching vias are not lost
    # between coarse lines (capped: on a large board the vias near the ports matter most)
    refs = {p.get('ref_net') for p in ports} - {None}
    holes = [d for d in board.drills if d.plated and d.net in refs]
    if len(holes) > 400:
        c = np.array([p['start'][:2] for p in ports])
        holes.sort(key=lambda d: np.min(np.hypot(c[:, 0] - d.x, c[:, 1] - d.y)))
        holes = holes[:400]
    X += [d.x for d in holes]
    Y += [d.y for d in holes]
    # signal vias: the barrel is a cylinder of the drill diameter; in the staircase mesh it stays connected to the
    # pads only with cells <= d/3 across it (a 0.3 mm via meshed at 0.7 mm was an open circuit, S21 -116 dB,
    # 2026-10-21), so every plated hole of a port net gets a uniform fine grid over the barrel plus a centre line
    for d in board.drills:
        if not d.plated or d.net not in pnets:
            continue
        # ... and over the pad + antipad (clearance hole in the planes): a staircase edge from a pad node to a
        # plane node shorts the via to GND when no line falls inside the clearance (0.3 mm via, 0.3 mm clearance
        # meshed at 0.4 mm: S11 0 dB, S21 -97 dB, 2026-10-21) -> lines every <= 0.1 mm out to d/2 + 0.8 mm
        r = max(min(d.d / 3, 0.1, res), 0.03)
        fine = min(fine, r * cf)
        R_ = d.d / 2 + 0.8
        X += [d.x] + list(np.arange(d.x - R_, d.x + R_ + 1e-9, r))
        Y += [d.y] + list(np.arange(d.y - R_, d.y + R_ + 1e-9, r))
    for p in ports:        # only the port plane gets a line; its transverse extent snaps to the strip lines
        (X if p['normal_axis'] == 0 else Y).append(p['start'][p['normal_axis']])
    for c in board.copper_layers:
        Z.append(board.SheetZ(c))
        if copper == 'solid':          # both copper faces, at least two cells across the outer copper
            z0, z1 = board.CopperZ(c)
            Z += [z0, z1] + ([0.5 * (z0 + z1)] if c in (board.copper_layers[0], board.copper_layers[-1]) else [])
    for l, z0, z1 in board.DielectricZ(copper):
        Z += list(np.linspace(z0, z1, 5))
    zres = min(res, 4 * (zt - board.SheetZ(board.copper_layers[1])) if len(board.copper_layers) > 1 else res)
    # merge distance: half the thirds-rule cell. Lines of neighbouring conductors closer than that only
    # shrink the timestep (fine-pitch pads gave 13 um cells, 8x more timesteps, same result)
    m = min(res / 20, fine / 2)
    return smooth(X, res, merge=m), smooth(Y, res, merge=m), smooth(Z, max(zres, 0.02), merge=0.005), res


def mesh_quality(lines, fmax, er=1.0, res=None):
    """ mesh_quality([X, Y, Z], fmax, er, res) -> dict

    The mesh figures a user reads before trusting a result: cell count, the smallest and largest
    cell per axis (mm), the largest neighbour ratio (grading), the largest cell aspect ratio of
    the tensor grid (max cell of one axis over min cell of another), and cells per wavelength
    at fmax in the dielectric (largest cell). ``grade`` is a one-word verdict: 'ok' when the
    grading is <= 1.6, the aspect <= 50 and the largest cell gives >= 10 cells/wavelength;
    'coarse' when the wavelength rule fails, 'rough' when the grading or the aspect does.
    """
    d = [np.diff(np.unique(np.asarray(v, float))) for v in lines]
    d = [x[x > 0] for x in d]
    if any(len(x) == 0 for x in d):
        return dict(grade='empty')
    mn = [float(x.min()) for x in d]
    mx = [float(x.max()) for x in d]
    ratio = max(float((x[1:] / x[:-1]).max()) if len(x) > 1 else 1.0 for x in d)
    ratio = max(ratio, max(float((x[:-1] / x[1:]).max()) if len(x) > 1 else 1.0 for x in d))
    aspect = max(mx[i] / mn[j] for i in range(3) for j in range(3) if i != j)
    lam = C0 / (fmax * math.sqrt(er)) / 1e-3
    cpw = lam / max(mx)
    q = dict(cells=int(np.prod([len(x) + 1 for x in d])), min_mm=[round(v, 4) for v in mn], max_mm=[round(v, 4) for v in mx],
             grading=round(ratio, 3), aspect=round(aspect, 1), cells_per_wavelength=round(cpw, 1))
    if res:
        q['max_over_res'] = round(max(mx) / res, 2)
    q['grade'] = 'ok' if (ratio <= 1.6 and aspect <= 50 and cpw >= 10) else ('coarse' if cpw < 10 else 'rough')
    return q


# ---------------------------------------------------------------- simulation
def final_mesh(board, ports, fmax, coarse=1.0, copper='sheet', refine=1.0):
    """ final_mesh(board, ports, fmax, ...) -> X, Y, Z, res

    The mesh exactly as build() hands it to the engine: mesh_lines() plus the PML padding, the
    lumped-component hints, the forced port planes and the grading. Shared with the twin server's
    /api/mesh (mesh preview in the KiCad plugin), so what the user sees is what gets simulated. """
    X, Y, Z, res = mesh_lines(board, ports, fmax, coarse=coarse, copper=copper, refine=refine)
    if any(str(b).startswith('PML') for b in BOUNDARY):
        (bx0, bx1), (by0, by1) = board.BoundingBox()[:2]
        zs = [board.SheetZ(c) for c in board.copper_layers]
        X, Y, Z = pml_pad(X, bx0, bx1), pml_pad(Y, by0, by1), pml_pad(Z, min(zs), max(zs))
    hints = getattr(board, 'mesh_hints', None) or {}
    if any(hints.get(a) for a in 'xyz'):          # lumped components: lines at the element ends
        X = smooth(list(X) + hints.get('x', []), res, merge=res / 40)
        Y = smooth(list(Y) + hints.get('y', []), res, merge=res / 40)
    # every port plane must lie exactly on a mesh line: the merge step of smooth() can fold it into a
    # neighbouring line (33 um away on gh_oresat_amp_thru), and a zero-thickness port box between lines
    # is empty (no excitation, probes read 0, S-parameters NaN)
    X, Y = list(X), list(Y)
    for p in ports:
        ax = p['normal_axis']
        _force_line(X if ax == 0 else Y, _sig7(p['start'][ax]))
    gr = grade_ratio()
    if gr:                                        # M4: no cell more than `gr` times its neighbour (coupon: +22 % cells, same dt)
        X, Y, Z = (grade(v, gr) for v in (X, Y, Z))
    return X, Y, Z, res


def build(board, ports, fmin, fmax, excite, copper='sheet', with_mask=False, coarse=1.0, nrts=200000, end=1e-4, materials=None,
          fields=None, refine=1.0):
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS
    # Lossless PEC loops (plane - via - plane) keep a static magnetic field from the DC part of the
    # excitation forever, so the energy end criterion may never be met. The ports' RF spectrum is
    # complete long before: cap the simulated time by the board's transit time.
    (bx0, bx1), (by0, by1) = board.BoundingBox()[:2]
    er = max([l.epsilon_r or 1.0 for l in board.stackup if l.type in ('prepreg', 'core')] + [1.0])
    t_max = max(5e-9, 25 * math.hypot(bx1 - bx0, by1 - by0) * 1e-3 * math.sqrt(er) / C0)
    FDTD = openEMS(EndCriteria=end, NrTS=nrts, MaxTime=t_max)
    if fmin > 0.2 * fmax:
        f0 = 0.5 * (fmin + fmax)
        FDTD.SetGaussExcite(f0, fmax - f0)
    else:
        # no DC in the excitation: a DC part charges floating copper and keeps a static field, so the energy
        # never falls below about -20 dB and every run went to MaxTime. The spectrum at DC is still -24 dB,
        # enough for the S-parameters down to fmax/100 (ratios of port waves).
        d = 0.025 * fmax
        FDTD.SetGaussExcite(0.5 * fmax + d, 0.5 * fmax - d)
    FDTD.SetBoundaryCond(BOUNDARY)
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mat = dict(MATERIALS, **(materials or {}))
    if mat.get('dispersion') == 'debye':
        mat.setdefault('band', (max(fmin, fmax / 50.0), fmax))
    board.mesh_hints = None      # a cached board keeps the hints of an earlier run with components
    board.ToCSX(CSX, frequency=0.5 * (fmin + fmax), copper=copper, with_mask=with_mask, **mat)
    X, Y, Z, res = final_mesh(board, ports, fmax, coarse=coarse, copper=copper, refine=refine)
    # Z2_EXTRA_PY=<file defining extra(CSX, X, Y, Z, res) -> (X, Y, Z)>: 3D parts the board file does not hold
    # (TH SMA body + coax, cable, instrument case); the port is then usually moved with Z2_PORT_FORCE
    xp = os.environ.get('Z2_EXTRA_PY')
    if xp:
        ns = {'__file__': xp}
        exec(compile(open(xp).read(), xp, 'exec'), ns)
        X, Y, Z = ns['extra'](CSX, list(X), list(Y), list(Z), res)
    g = CSX.GetGrid()
    g.SetDeltaUnit(1e-3)
    for ax, v in zip('xyz', (X, Y, Z)):
        # Write2XML stores primitive coordinates with 7 significant digits, mesh lines at full
        # precision; a sheet or port that is not exactly on its mesh line is dropped ("unused primitive").
        # Lines and port coordinates are rounded to the same 7 significant digits, so they match exactly.
        # (5 decimals did not: an assumed stackup puts In1 at 1.0783333 mm, written as 1.078333 for the
        # sheet but meshed at 1.07833, and both inner planes vanished on gh_oresat_amp_thru.)
        g.SetLines(ax, np.unique([_sig7(x) for x in v]))
    pobj = []
    for p in ports:
        st = [_sig7(v) for v in p['start']]
        sp = [_sig7(v) for v in p['stop']]
        pobj.append(FDTD.AddLumpedPort(p['number'], p.get('z0', 50.0), st, sp, 'z', excite=1 if p['number'] == excite else 0))
    nl = [len(g.GetLines(a)) for a in 'xyz']
    dumps = []
    if fields and excite == 1:                    # current maps of the port-1 excitation
        import fields as fmod
        dumps = fmod.add_dumps(CSX, board, fields['freqs'], fields.get('layers'))
    er = max([l.epsilon_r or 1.0 for l in board.stackup if l.type in ('prepreg', 'core')] + [1.0])
    quality = mesh_quality([g.GetLines(a) for a in 'xyz'], fmax, er, res)
    return FDTD, CSX, pobj, dict(mesh=nl, mcells=float(np.prod(nl) / 1e6), res_mm=res, dumps=dumps, quality=quality)


_CUDA_OK = None


def cuda_available():
    """ True when the engine has the CUDA engine (its --help lists it) and an NVIDIA driver is installed. """
    global _CUDA_OK
    if _CUDA_OK is None:
        _CUDA_OK = False
        exe = VARIANTS['cuda']['exe']
        try:
            env = dict(os.environ, PATH=os.pathsep.join(VARIANTS['cuda']['path']) + os.pathsep + os.environ.get('PATH', ''))
            h = subprocess.run([exe, '--help'], env=env, capture_output=True, text=True, timeout=60)
            if 'cuda' in (h.stdout + h.stderr).lower():
                import ctypes
                if os.name == 'nt':
                    ctypes.WinDLL('nvcuda.dll')
                else:
                    ctypes.CDLL('libcuda.so.1')
                _CUDA_OK = True
        except (OSError, subprocess.TimeoutExpired):
            pass
    return _CUDA_OK


_ENGINE_OPTS = {}


def engine_default_args(exe, env):
    """ Default engine options the twin adds to every run when this engine build knows them:
        conformal PEC (contour path + free face area, a=0.5), which brings a 50 ohm coupon with the thirds
        rule from 52.8 to 50.9 ohm (reference 49.0) at 4 lines across the strip and converges 2-3x faster in
        the lines across. Z2_CONFORMAL=0 turns it off; Z2_OPENEMS_ARGS still adds / overrides. """
    if os.environ.get('Z2_CONFORMAL', '1') in ('0', 'off', 'no'):
        return []
    if exe not in _ENGINE_OPTS:
        try:
            h = subprocess.run([exe, '--help'], env=env, capture_output=True, text=True, timeout=60).stdout
        except (OSError, subprocess.TimeoutExpired):
            h = ''
        _ENGINE_OPTS[exe] = ['--conformal-pec=1', '--conformal-min-area=0.5'] if 'conformal-pec' in h else []
    return list(_ENGINE_OPTS[exe])


class Watchdog(object):
    """ Early stop of a run that is going nowhere, from the engine's progress lines.

        unstable: energy NaN/inf, or more than grow_db above its lowest value after the excitation
                  -> the engine is killed, the run is marked failed (nothing worth evaluating)
        stagnant: after min_exc excitation lengths the energy has not dropped by stall_db over a window of
                  max(stall_frac of the elapsed timesteps, stall_exc excitation lengths), still above end_db + margin_db
                  -> an "ABORT" file stops openEMS gracefully: the port data are complete up to there
                  (a static field from floating copper or a lossless loop; the RF part has decayed)
    """
    RX = re.compile(r'Timestep:\s*(\d+).*Energy:\s*~?\s*([-+\w.]+)\s*\(\s*([-+]?\s*[\w.]+)\s*dB\)')    # "(- 0.00dB)" has a space

    def __init__(self, sim, end_db=-40.0, min_exc=3.0, stall_frac=0.25, stall_db=0.5, margin_db=3.0, grow_db=10.0,
                 min_points=8, stall_exc=3.0):
        self.sim, self.end_db, self.min_exc = sim, end_db, min_exc
        self.stall_frac, self.stall_db, self.margin_db, self.grow_db = stall_frac, stall_db, margin_db, grow_db
        self.min_points, self.stall_exc = min_points, stall_exc
        self.exc = None
        self.hist = []
        self.verdict = None
        self.aborted = False

    def __call__(self, line):
        m = re.search(r'Excitation signal length is:\s*(\d+)', line)
        if m:
            self.exc = int(m.group(1))
            return
        m = self.RX.search(line)
        if not m:
            return
        ts, db = int(m.group(1)), m.group(3).replace(' ', '')
        try:
            db = float(db)
        except ValueError:
            db = float('nan')
        if not np.isfinite(db) or 'nan' in m.group(2).lower() or 'inf' in m.group(2).lower():
            self.verdict = 'unstable: energy %s at timestep %d' % (m.group(2), ts)
            raise RuntimeError(self.verdict)
        # an explosion during the excitation hides behind the relative dB (the running maximum resets):
        # the absolute energy must not jump by a million between two reports
        try:
            e_abs = float(m.group(2))
        except ValueError:
            e_abs = float('nan')
        if np.isfinite(e_abs) and self.hist and self.hist[-1][2] > 0 and e_abs > 1e6 * self.hist[-1][2]:
            self.verdict = 'unstable: energy %s at timestep %d (x%.0e in one report)' % (m.group(2), ts, e_abs / self.hist[-1][2])
            raise RuntimeError(self.verdict)
        self.hist.append((ts, db, e_abs))
        exc = self.exc or 0
        if ts < self.min_exc * max(exc, 1):
            return
        after = [d for t, d, _ in self.hist if t >= exc]
        if after and db > min(after) + self.grow_db:
            self.verdict = 'unstable: energy grew %.1f dB above its minimum (timestep %d)' % (db - min(after), ts)
            raise RuntimeError(self.verdict)
        if self.aborted or db <= self.end_db + self.margin_db:
            return
        # a plateau, not a slow decay: high-Q arrays lose ~1 dB per few thousand steps and still reach the end
        # criterion; stopping those truncates the ringing and shifts the resonance (2026-10-13, 24 GHz arrays)
        span = max(ts * self.stall_frac, self.stall_exc * exc)
        win = [(t, d) for t, d, _ in self.hist if t >= ts - span]
        if len(win) >= self.min_points and win[-1][0] - win[0][0] >= 0.9 * span and win[0][1] - db < self.stall_db:
            self.verdict = 'stagnant: %.1f dB, %.2f dB drop over the last %d timesteps; stopped at %d' % (
                db, win[0][1] - db, ts - win[0][0], ts)
            open(os.path.join(self.sim, 'ABORT'), 'w').close()      # openEMS stops and writes its data
            self.aborted = True


def engine_figures(txt):
    """ what the operator printed: time step (ps), PEC edges, conformal (partly metal) edges """
    out = {}
    m = re.search(r'Timestep \(s\)\s*:\s*([\d.eE+-]+)', txt) or re.search(r'FDTD timestep is:\s*([\d.eE+-]+)\s*s', txt)
    if m:
        out['dt_ps'] = round(float(m.group(1)) * 1e12, 5)
    m = re.search(r'Number of PEC edges\s*:\s*(\d+)', txt)
    if m:
        out['pec_edges'] = int(m.group(1))
    m = re.search(r'Conformal PEC edges\s*:\s*(\d+)', txt)
    if m:
        out['conformal_edges'] = int(m.group(1))
        if out.get('pec_edges'):
            out['conformal_frac'] = round(out['conformal_edges'] / (out['pec_edges'] + out['conformal_edges']), 4)
    return out


GUARD_FAIL = re.compile(r'MEMORY GUARD FAILED|non-finite field energy|GPU memory error')


def gpu_fault(r, sim):
    """ True when a run died of a GPU fault: the engine's memory guard, or the watchdog's explosion verdict on
        the CUDA engine (a bit error on a worn GPU blows the field up within a few hundred steps). """
    if r.get('variant') != 'cuda':
        return False
    if (r.get('watchdog') or '').startswith('unstable'):
        return True
    try:
        with open(os.path.join(sim, 'stdout.txt'), errors='replace') as fh:
            txt = fh.read()
    except OSError:
        return False
    if r.get('rc', 0) != 0 and GUARD_FAIL.search(txt):
        return True
    # an engine without the memory guard reports a blown-up run as "end-criteria ... reached (-nandB)"
    return bool(re.search(r'RunFDTD:.*\(-?nan', txt))


def run_engine(xml, sim, variant=None, threads=4, timeout=6 * 3600, progress=None, watchdog=None, retries=None):
    """ run_engine(...) -> dict: the engine run, repeated after a GPU fault.

    watchdog: None, or True / a dict of Watchdog options: early stop (see Watchdog); the verdict is in the
    returned dict ('watchdog'), an unstable run comes back with rc=-9 and 'failed' set.
    retries: how often a CUDA run is repeated after a GPU fault (the engine's memory guard, or an explosion
    caught by the watchdog), default Z2_GPU_RETRIES or 2. Old or marginal GPUs corrupt a value now and then
    (a GTX 1080 Ti: 1 run in 8); a repeat costs the same run again, a missed fault costs a wrong result.
    'attempts' in the result counts the runs, 'faults' lists the verdicts of the failed ones. """
    if retries is None:
        retries = int(os.environ.get('Z2_GPU_RETRIES', '2'))
    faults = []
    for attempt in range(retries + 1):
        r = _run_engine_once(xml, sim, variant, threads, timeout, progress, watchdog)
        if not gpu_fault(r, sim) or attempt == retries:
            break
        why = r.get('watchdog') or ('engine memory guard (rc=%s)' % r.get('rc') if r.get('rc', 0) != 0 else 'blown up (NaN energy)')
        faults.append(why)
        print('twin: GPU fault (%s), running again (%d/%d)' % (why, attempt + 1, retries), flush=True)
        try:
            shutil.copy(os.path.join(sim, 'stdout.txt'), os.path.join(sim, '..', 'stdout_fault%d.txt' % (attempt + 1)))
        except OSError:
            pass
    r['attempts'] = len(faults) + 1
    if faults:
        r['faults'] = faults
    return r


def _run_engine_once(xml, sim, variant=None, threads=4, timeout=6 * 3600, progress=None, watchdog=None):
    used = variant or DEFAULT_VARIANT
    if used == 'auto':
        used = 'cuda' if cuda_available() else DEFAULT_VARIANT
    v = VARIANTS[used]
    if os.path.exists(sim):
        shutil.rmtree(sim)
    os.makedirs(sim)
    shutil.copy(xml, os.path.join(sim, 'model.xml'))
    env = dict(os.environ, PATH=os.pathsep.join(v['path']) + os.pathsep + os.environ['PATH'])
    if v.get('libs'):
        env['LD_LIBRARY_PATH'] = os.pathsep.join(v['libs'] + [os.environ.get('LD_LIBRARY_PATH', '')])
    t0 = time.time()
    base = v.get('engine', ['--engine=multithreaded', '--numThreads=%d' % threads])
    extra = os.environ.get('Z2_OPENEMS_ARGS', '').split()        # A/B runs; an option given here replaces the default one (engine too)
    given = {a.split('=')[0] for a in extra}
    cmd = [v['exe'], 'model.xml'] + [a for a in base + engine_default_args(v['exe'], env) if a.split('=')[0] not in given] + extra
    eff_args = ' '.join(cmd[2:])
    wd = None
    if watchdog:
        wd = Watchdog(sim, **(watchdog if isinstance(watchdog, dict) else {}))
        user = progress
        progress = (lambda line: (wd(line), user(line) if user else None))
    with open(os.path.join(sim, 'stdout.txt'), 'w') as fh:
        if progress is None:
            p = subprocess.run(cmd, cwd=sim, env=env, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout)
        else:                       # stream the engine output: progress(line) for every line; an exception
            p = subprocess.Popen(cmd, cwd=sim, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace')
            try:                    # raised by progress() (a cancel) stops the engine
                for line in p.stdout:
                    fh.write(line)
                    try:
                        progress(line.rstrip())
                    except RuntimeError:
                        if wd is None or not (wd.verdict or '').startswith('unstable'):
                            raise
                        fh.write('twin watchdog: %s\n' % wd.verdict)
                        p.kill()
                        break
                    if time.time() - t0 > timeout:
                        p.kill()
                p.wait()
            finally:
                if p.poll() is None:
                    p.kill()
                    p.wait()
    txt = open(os.path.join(sim, 'stdout.txt'), errors='replace').read()
    sp = re.findall(r'Speed:\s*([\d.]+)\s*MCells/s', txt)
    ts = re.findall(r'Time for\s+(\d+)\s+iterations', txt)
    out = dict(rc=p.returncode, wall_s=round(time.time() - t0, 2), mcells_s=float(sp[-1]) if sp else None,
               timesteps=int(ts[-1]) if ts else None, variant=used, args=eff_args)
    out.update(engine_figures(txt))
    if wd is not None:
        out['watchdog'] = wd.verdict
        out['failed'] = bool(wd.verdict and wd.verdict.startswith('unstable'))
        if wd.verdict:
            with open(os.path.join(sim, 'stdout.txt'), 'a') as fh:
                fh.write('twin watchdog: %s\n' % wd.verdict)
    return out


def write_touchstone(path, f, S, z0=50.0, comment=''):
    n = S.shape[1]
    with open(path, 'w') as fh:
        for line in comment.splitlines():
            fh.write('! %s\n' % line)
        fh.write('# Hz S RI R %g\n' % z0)
        for k, fk in enumerate(f):
            vals = []
            # Touchstone order: 2-port S11 S21 S12 S22; N-port row by row
            idx = [(0, 0), (1, 0), (0, 1), (1, 1)] if n == 2 else [(i, j) for i in range(n) for j in range(n)]
            for i, j in idx:
                vals += ['%.9e' % S[k, i, j].real, '%.9e' % S[k, i, j].imag]
            fh.write('%.6e %s\n' % (fk, ' '.join(vals)))


def line_params(f, S, length_mm, z_ref=50.0):
    """ Z0 and eeff of a uniform 2-port line from its S matrix (ABCD conversion). """
    S11, S21, S12, S22 = S[:, 0, 0], S[:, 1, 0], S[:, 0, 1], S[:, 1, 1]
    d = 2 * S21
    A = ((1 + S11) * (1 - S22) + S12 * S21) / d
    B = z_ref * ((1 + S11) * (1 + S22) - S12 * S21) / d
    C = ((1 - S11) * (1 - S22) - S12 * S21) / (d * z_ref)
    Z0 = np.sqrt(B / C)
    Z0 = np.where(Z0.real < 0, -Z0, Z0)
    beta = -np.unwrap(np.angle(S21 * np.exp(0j))) / (length_mm * 1e-3)
    eeff = (beta * C0 / (2 * np.pi * f)) ** 2
    return Z0, eeff


def simulate(board_file, port_refs, fmin, fmax, out, fab=None, variant=None, threads=4, copper='sheet',
             with_mask=False, coarse=1.0, nf=401, log=None, materials=None, ports=None, progress=None, z0=None, fields=None,
             board_obj=None, converge=False, refine=1.0):
    """ Simulate a KiCad board between the given port pads (footprint references), or explicit port
        definitions `ports` (as returned by make_port / snap_port). Returns the result dict (also result.json). """
    log = log or (lambda m: print(m, flush=True))
    _setup_python()
    from CSXCAD.KiCadImport import KiCadBoard
    from CSXCAD import FabStackups
    board = board_obj if board_obj is not None else KiCadBoard(board_file)     # board_obj: pre-processed (e.g. cropped)
    if fab:
        fn, nm = fab.split(':', 1)
        FabStackups.Apply(board, FabStackups.Find(fab=fn, name=nm)[0])
    if ports is None:
        ports = [make_port(board, find_pad(board, r), i + 1) for i, r in enumerate(port_refs)]
    for i, p in enumerate(ports):
        p['number'] = i + 1
        p.setdefault('z0', z0[i] if z0 else 50.0)
    os.makedirs(out, exist_ok=True)
    f = np.linspace(fmin, fmax, nf)
    n = len(ports)
    S = np.zeros((nf, n, n), complex)
    runs = []
    meta = None
    field_res = None
    for e in range(1, n + 1):
        FDTD, CSX, pobj, meta = build(board, ports, fmin, fmax, e, copper=copper, with_mask=with_mask, coarse=coarse, refine=refine, materials=materials,
                                      fields=fields)
        cdir = os.path.join(out, 'exc%d' % e)
        os.makedirs(cdir, exist_ok=True)
        xml = os.path.join(cdir, 'model.xml')
        FDTD.Write2XML(xml)
        q = meta.get('quality', {})
        log('excitation %d/%d: mesh %s (%.2f MCells), engine %s' % (e, n, meta['mesh'], meta['mcells'], variant or DEFAULT_VARIANT))
        if q.get('grade'):
            log('  mesh quality %s: min %s mm, max %s mm, grading %s, aspect %s, %s cells/wavelength' % (
                q['grade'], q['min_mm'], q['max_mm'], q['grading'], q['aspect'], q['cells_per_wavelength']))
        r = run_engine(xml, os.path.join(cdir, 'sim'), variant, threads,
                       progress=(lambda line, e=e: progress(e, n, line)) if progress else None)
        log('  rc=%s %.1f s, %s timesteps, %s MCells/s' % (r['rc'], r['wall_s'], r['timesteps'], r['mcells_s']))
        if r['rc'] != 0:
            raise RuntimeError('openEMS failed, see %s' % os.path.join(cdir, 'sim', 'stdout.txt'))
        for p, pd in zip(pobj, ports):
            p.CalcPort(os.path.join(cdir, 'sim'), f, ref_impedance=pd.get('z0', 50.0))
        a = pobj[e - 1].uf_inc
        for i, p in enumerate(pobj):
            S[:, i, e - 1] = p.uf_ref / a
        runs.append(r)
        if meta.get('dumps'):
            import fields as fmod
            try:
                field_res = fmod.export(os.path.join(cdir, 'sim'), os.path.join(out, 'fields'), board, meta['dumps'])
            except Exception as ex:                    # the S-parameters stay valid without the maps
                field_res = dict(error='%s: %s' % (type(ex).__name__, ex))
            log('  current maps: %s' % ('%d layers' % len(field_res.get('maps', {})) if field_res and 'maps' in field_res else field_res))
    src = os.path.splitext(os.path.basename(board_file))[0]
    write_touchstone(os.path.join(out, '%s.s%dp' % (src, n)), f, S,
                     comment='openEMS-z2labs kicad_sim, %s\nports: %s\nstackup: %s\ncopper: %s, mask: %s' % (
                         os.path.basename(board_file), ', '.join('%d=%s (%s over %s)' % (p['number'], p['ref'], p['layer'], p['ref_layer']) for p in ports),
                         board.stackup_source, copper, with_mask))
    res = dict(board=os.path.basename(board_file), stackup=board.stackup_source, copper=copper, with_mask=with_mask,
               fields=field_res,
               materials={k: v for k, v in dict(MATERIALS, **(materials or {})).items()},
               components=[dict(ref=c['ref'], kind=c['kind'], value=c['value']) for c in getattr(board, 'components', [])],
               ports=ports, mesh=meta['mesh'], mcells=meta['mcells'], res_mm=meta['res_mm'], runs=runs,
               mesh_quality=dict(meta.get('quality', {}), **{k: runs[0][k] for k in ('dt_ps', 'pec_edges', 'conformal_edges', 'conformal_frac') if runs and k in runs[0]}),
               S_re={'S%d%d' % (i + 1, j + 1): S[:, i, j].real.round(6).tolist() for i in range(n) for j in range(n)},
               S_im={'S%d%d' % (i + 1, j + 1): S[:, i, j].imag.round(6).tolist() for i in range(n) for j in range(n)},
               f_GHz=(f / 1e9).round(5).tolist(),
               S_dB={'S%d%d' % (i + 1, j + 1): (20 * np.log10(np.abs(S[:, i, j]) + 1e-15)).round(3).tolist() for i in range(n) for j in range(n)})
    if n == 2:
        L = math.dist(ports[0]['start'][:2], ports[1]['start'][:2]) if ports[0]['normal_axis'] == ports[1]['normal_axis'] else None
        if L:
            Z0, ee = line_params(f, S, L)
            band = (f >= 0.15 * fmax) & (f <= 0.7 * fmax)
            res.update(line_length_mm=L, Z0_ohm=float(np.median(Z0.real[band])), eeff=float(np.median(ee[band])),
                       Z0_curve=Z0.real.round(3).tolist(), eeff_curve=ee.round(4).tolist())
    if converge:
        # M3: the same board with 0.6x cells; the differences say whether this mesh can be trusted
        import convergence
        fine = simulate(board_file, port_refs, fmin, fmax, os.path.join(out, 'fine'), fab=fab, variant=variant, threads=threads,
                        copper=copper, with_mask=with_mask, coarse=coarse * 0.6, nf=nf, log=log, materials=materials, ports=ports,
                        progress=progress, z0=z0, fields=None, board_obj=board_obj, converge=False, refine=refine * 0.6)
        res['convergence'] = convergence.compare(res, fine)
        log('convergence: ' + res['convergence']['verdict'])
    json.dump(res, open(os.path.join(out, 'result.json'), 'w'), indent=1)
    return res


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    ap.add_argument('board')
    ap.add_argument('--port', action='append', help='footprint reference of a port pad (repeat); without it the RF '
                    'connectors found by rf_ports are used')
    ap.add_argument('--list-ports', action='store_true', help='only list the RF port candidates and exit')
    ap.add_argument('--fmin', type=float, default=0.1e9)
    ap.add_argument('--fmax', type=float, default=6e9)
    ap.add_argument('--fab', help='library stackup, e.g. JLCPCB:JLC04161H-7628')
    ap.add_argument('--out')
    ap.add_argument('--variant', default=DEFAULT_VARIANT, choices=sorted(VARIANTS) + ['auto'],
                    help='engine build; auto = the GPU build when it exists, else %s' % DEFAULT_VARIANT)
    ap.add_argument('--dispersion', choices=('debye',), help='Djordjevic-Sarkar laminates (flat tanD, causal eps_r(f))')
    ap.add_argument('--loss', action='store_true', help='copper as conducting sheets (conductor loss)')
    ap.add_argument('--roughness', type=float, help='copper RMS roughness Rq in um (default: from the fab stackup)')
    ap.add_argument('--components', action='store_true', help='2-pad R/L/C footprints as lumped elements')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--copper', default='sheet', choices=('sheet', 'solid'))
    ap.add_argument('--mask', action='store_true', help='include the solder mask layers')
    ap.add_argument('--coarse', type=float, default=1.0, help='mesh coarsening factor (>1: faster, less accurate)')
    ap.add_argument('--converge', action='store_true', help='run again with 0.6x cells and report whether the mesh is converged')
    a = ap.parse_args(argv)
    if not a.port or a.list_ports:
        import rf_ports
        _setup_python()
        from CSXCAD.KiCadImport import KiCadBoard
        board = KiCadBoard(a.board)
        cands = rf_ports.candidates(board, rf_ports.load_project(a.board))
        print(rf_ports.report(board, rf_ports.load_project(a.board)), file=sys.stderr)
        if a.list_ports:
            return
        a.port = [c['ref'] for c in cands if c['confident']]
        if not a.port:
            sys.exit('no confident RF port found: give --port REF')
    out = a.out or os.path.splitext(a.board)[0] + '_sim'
    mat = dict(dispersion=a.dispersion, conductor='loss' if a.loss else 'pec', roughness_um=a.roughness, components=a.components)
    r = simulate(a.board, a.port, a.fmin, a.fmax, out, a.fab, a.variant, a.threads, a.copper, a.mask, a.coarse, materials=mat, converge=a.converge)
    keys = ('board', 'stackup', 'copper', 'mesh', 'mcells', 'mesh_quality', 'Z0_ohm', 'eeff', 'line_length_mm', 'convergence')
    print(json.dumps({k: r[k] for k in keys if k in r}, indent=1))


if __name__ == '__main__':
    main(sys.argv[1:])
