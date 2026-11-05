#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2026 Zoltan Doczi -- GPL-3.0-or-later
"""
sma_model: a simple VRML (IndexedFaceSet) model of a cheap edge-mount SMA jack for the KiCad 3D
viewer: square body straddling the board edge, threaded barrel, PTFE face, round centre pin lying
on the pad, four ground legs. Dimensions in mm; KiCad VRML units are 0.1 inch (2.54 mm).

    write_wrl(path, pin=1.27, lp=3.0, body=6.4, barrel_len=5.5, board=1.6)

The model origin is the footprint origin (centre of the signal pad); the board edge is at x = -lp/2,
the connector extends to -x. Rotate the footprint 180 deg for the opposite edge.
"""
import math

SCALE = 1 / 2.54


def _box(x0, y0, z0, x1, y1, z1):
    v = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0), (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    f = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    return v, f


def _cyl_x(x0, x1, yc, zc, r, n=24):
    v = []
    for x in (x0, x1):
        for i in range(n):
            a = 2 * math.pi * i / n
            v.append((x, yc + r * math.cos(a), zc + r * math.sin(a)))
    f = []
    for i in range(n):
        j = (i + 1) % n
        f.append((i, j, n + j, n + i))
    f.append(tuple(range(n - 1, -1, -1)))
    f.append(tuple(range(n, 2 * n)))
    return v, f


def _shape(verts, faces, rgb):
    pts = ', '.join('%.4f %.4f %.4f' % (x * SCALE, y * SCALE, z * SCALE) for x, y, z in verts)
    idx = ', '.join(' '.join(str(i) for i in f) + ' -1' for f in faces)
    return ('Shape { appearance Appearance { material Material { diffuseColor %.2f %.2f %.2f specularColor 0.6 0.6 0.6 shininess 0.5 } }\n'
            ' geometry IndexedFaceSet { solid FALSE coord Coordinate { point [ %s ] } coordIndex [ %s ] } }\n' % (rgb + (pts, idx)))


def write_wrl(path, pin=1.27, lp=3.0, body=6.4, barrel_len=5.5, board=1.6, leg_w=1.6, leg_len=3.0):
    xe = -lp / 2                                   # board edge
    zc = -board / 2                                # body centred on the board thickness
    out = ['#VRML V2.0 utf8\n']
    v, f = _box(xe - body, -body / 2, zc - body / 2, xe, body / 2, zc + body / 2)
    out.append(_shape(v, f, (0.75, 0.72, 0.55)))                                   # brass body
    v, f = _cyl_x(xe - body - barrel_len, xe - body, 0, zc, 3.15)
    out.append(_shape(v, f, (0.78, 0.75, 0.6)))                                    # threaded barrel
    v, f = _cyl_x(xe - body - barrel_len - 0.01, xe - body - barrel_len + 1.5, 0, zc, 2.1)
    out.append(_shape(v, f, (0.95, 0.95, 0.95)))                                   # PTFE in the barrel
    v, f = _cyl_x(xe - 0.5, lp / 2 - 0.3, 0, pin / 2, pin / 2)
    out.append(_shape(v, f, (0.85, 0.7, 0.3)))                                     # centre pin on the pad
    gy = leg_w / 2 + 0.5 + 1.6 / 2 if lp else 0
    for sy in (-1, 1):
        for z0, z1 in ((0.0, 0.3), (-board - 0.3, -board)):
            v, f = _box(xe, sy * (body / 2) - (leg_w if sy > 0 else 0), z0, xe + leg_len, sy * (body / 2) + (leg_w if sy < 0 else 0), z1)
            out.append(_shape(v, f, (0.75, 0.72, 0.55)))                           # ground legs
    open(path, 'w').write(''.join(out))
    return path


if __name__ == '__main__':
    import sys
    write_wrl(sys.argv[1] if len(sys.argv) > 1 else 'sma_edge.wrl')
