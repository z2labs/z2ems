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
fields: per-copper-layer current density maps for the board viewer.

A frequency-domain rot(H) dump (openEMS dump type 13, the total current density)
is placed on the sheet plane of every copper layer. On a copper sheet rot(H) is
dominated by the surface current, so the magnitude map shows where the signal
and its return current flow -- e.g. a return path forced around a plane split.

The maps are resampled to a uniform pixel grid over the board and written as
PNG heat maps (dB relative to the strongest current of the excitation on any
layer, so the layers compare directly). numpy only; h5py to read the dumps.
"""

import math
import os
import struct
import zlib

import numpy as np


def add_dumps(CSX, board, freqs, layers=None, margin=0.3):
    """ One FD rot(H) dump per copper sheet plane over the board area. Returns [(layer, dump name)]. """
    (x0, x1), (y0, y1) = board.BoundingBox()[:2]
    out = []
    for c in layers or board.copper_layers:
        z = board.SheetZ(c)
        name = 'J_' + c.replace('.', '_')
        d = CSX.AddDump(name, dump_type=13, dump_mode=0, file_type=1, frequency=[float(f) for f in freqs])
        d.AddBox([x0 - margin, y0 - margin, z], [x1 + margin, y1 + margin, z])
        out.append((c, name))
    return out


def _read(path):
    """ (x, y, F, J) of a 2D FD dump: mesh lines (in drawing units) and J[f, comp, ny, nx] complex. """
    import h5py
    try:
        from openEMS.utilities import HDF5Dump            # the fork's reader handles every layout
        with HDF5Dump(path) as d:
            mesh = d.GetMesh()
            x, y = (np.asarray(l, float) / (mesh.get('scaling') or 1.0) for l in mesh['lines'][:2])   # drawing units (mm)
            F = np.asarray(d.frequencies, float)
            J = []
            for k in range(len(F)):
                v = np.asarray(d.GetFieldAtIndex(f_idx=k))            # (3, nx, ny, nz)
                v = v[..., 0] if v.ndim == 4 else v
                J.append(np.transpose(v, (0, 2, 1)))                   # -> (3, ny, nx)
            return x, y, F, np.asarray(J)
    except (ImportError, AttributeError):
        # no fork bindings, or the official openEMS package whose HDF5Dump has no `frequencies`
        # (the official openEMS release wheels lack it: every field export failed with AttributeError)
        pass
    with h5py.File(path, 'r') as h:                                  # legacy layout: f<k>_real (3, nz, ny, nx)
        x, y = np.asarray(h['Mesh/x']), np.asarray(h['Mesh/y'])
        fd = h['FieldData/FD']
        F = np.asarray(fd.attrs.get('frequency', []), float)
        J = []
        k = 0
        while 'f%d_real' % k in fd:
            a = np.asarray(fd['f%d_real' % k]) + 1j * np.asarray(fd['f%d_imag' % k])
            J.append(a[:, 0] if a.ndim == 4 else a)
            k += 1
        # openEMS HDF5 v0.3 (official openEMS releases): complex datasets f<k> with d_order NXYZ, i.e. (3, nx, ny, nz)
        while 'f%d' % k in fd:
            a = np.asarray(fd['f%d' % k])
            if a.ndim == 4:
                a = a[..., 0]
            order = fd['f%d' % k].attrs.get('d_order', 'NXYZ')
            order = order.decode() if isinstance(order, bytes) else str(order)
            if order == 'NXYZ':
                a = np.transpose(a, (0, 2, 1))                                    # -> (3, ny, nx)
            J.append(a)
            k += 1
        if not J:
            raise ValueError('no field data in %s (keys: %s)' % (path, list(fd.keys())[:6]))
        return x, y, F, np.asarray(J)


# a perceptually ordered dark-to-bright map (inferno-like anchors)
_ANCHORS = np.array([[0, 0, 4], [40, 11, 84], [101, 21, 110], [159, 42, 99], [212, 72, 66], [245, 125, 21],
                     [250, 193, 39], [252, 255, 164]], float)


def colormap(v):
    """ v in [0, 1] -> uint8 RGB. """
    v = np.clip(np.asarray(v, float), 0, 1) * (len(_ANCHORS) - 1)
    i = np.minimum(v.astype(int), len(_ANCHORS) - 2)
    t = (v - i)[..., None]
    return (_ANCHORS[i] * (1 - t) + _ANCHORS[i + 1] * t).astype(np.uint8)


def write_png(path, rgb):
    """ Minimal RGB PNG writer (rgb: uint8 [h, w, 3], row 0 = top). """
    h, w = rgb.shape[:2]
    raw = b''.join(b'\x00' + rgb[r].tobytes() for r in range(h))

    def chunk(t, d):
        c = struct.pack('>I', len(d)) + t + d
        return c + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
    with open(path, 'wb') as fh:
        fh.write(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 2, 0, 0, 0))
                 + chunk(b'IDAT', zlib.compress(raw, 6)) + chunk(b'IEND', b''))


def export(sim_dir, out_dir, board, dumps, db_range=40.0, px_per_mm=6.0, max_px=1400):
    """ Turn the dumps of one run into PNG heat maps; returns the description for result.json:
        dict(f_GHz=[...], db_range, bbox=[x0, y0, x1, y1], maps={layer: [png per frequency]}, peak_db={layer: [...]}) """
    os.makedirs(out_dir, exist_ok=True)
    (bx0, bx1), (by0, by1) = board.BoundingBox()[:2]
    W = int(min(max_px, max(64, (bx1 - bx0) * px_per_mm)))
    H = int(min(max_px, max(64, (by1 - by0) * px_per_mm)))
    xs = np.linspace(bx0, bx1, W)
    ys = np.linspace(by1, by0, H)                 # image row 0 = top = max y
    mags, F = {}, None
    for layer, name in dumps:
        p = os.path.join(sim_dir, name + '.h5')
        if not os.path.exists(p):
            continue
        x, y, F, J = _read(p)
        if np.ptp(x) < (bx1 - bx0) / 100:                         # mesh stored in metres
            x, y = x * 1e3, y * 1e3
        m = np.sqrt(np.sum(np.abs(J) ** 2, axis=1))                  # |J| per frequency: (nf, ny, nx)
        ix = np.clip(np.searchsorted(x, xs) - 1, 0, len(x) - 1)
        iy = np.clip(np.searchsorted(y, ys) - 1, 0, len(y) - 1)
        mags[layer] = m[:, iy][:, :, ix]
    if not mags:
        return None
    rel = os.path.basename(out_dir.rstrip('/\\'))
    res = dict(f_GHz=[round(float(f) / 1e9, 4) for f in F], db_range=db_range, bbox=[bx0, by0, bx1, by1], maps={}, peak_db={})
    for k in range(len(F)):
        ref = max(float(m[k].max()) for m in mags.values()) or 1.0  # common reference: layers compare directly
        for layer, m in mags.items():
            db = 20 * np.log10(np.maximum(m[k], 1e-30) / ref)
            img = colormap((db + db_range) / db_range)
            fn = '%s_%dMHz.png' % (layer.replace('.', '_'), round(F[k] / 1e6))
            write_png(os.path.join(out_dir, fn), img)
            res['maps'].setdefault(layer, []).append(rel + '/' + fn)
            res['peak_db'].setdefault(layer, []).append(round(float(db.max()), 1))
    return res
