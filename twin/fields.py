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
fields: surface currents, voltages and E-fields of a board run, for the 3D viewer and KiCad overlays.

Dumps (frequency domain, every excitation):
  J_<layer>  rot(H) on the sheet plane of each copper layer (openEMS dump type 13). On a copper sheet
             rot(H) times the dual cell height is the surface current density Js (A/m): the signal
             current and its return current in the planes.
  E_vol      E (dump type 10) in a box around the whole stack at a reduced resolution, so any cut plane
             can be looked at after the run (field lines from a trace to its reference plane).

Everything is normalised to 1 W incident at the driven port (peak phasors, phase 0 = incident wave at the
port), so the numbers mean something: a matched 50 ohm line carries 0.2 A and 10 V. The voltage of a copper
layer against the next copper layer below or above it is V = -integral(Ez dz) between the two sheets.

Export (fields/ next to result.json):
  J_e<p>_f<k>.f16    float16 [n_layers, ny, nx, 4]   Jx re, Jx im, Jy re, Jy im (A/m / scale)
  E_e<p>_f<k>.f16    float16 [nz, ny, nx, 6]         Ex, Ey, Ez re/im (V/m / scale) on the dump grid
  V_e<p>_f<k>.f16    float16 [n_pairs, ny, nx, 2]    V re/im (V / scale) on the E grid (x, y)
  mask.u8            uint8   [n_layers, ny, nx]      copper of each layer on the J grid
  <layer>_<MHz>.png  |J| heat maps of the first excitation (the KiCad plugin 0.2 overlays)
and the index in result.json "fields" (version 2).

streamlines() and contour_lines() turn a field into polylines (KiCad overlays, numpy only).
"""

import math
import os
import struct
import zlib

import numpy as np

VERSION = 2
POWER_W = 1.0                  # incident power the fields are scaled to


# ---------------------------------------------------------------- dumps
def add_dumps(CSX, board, freqs, layers=None, margin=0.3, volume=True, z_lines=None):
    """ One FD rot(H) dump per copper sheet plane and one E volume dump. Returns a list of dicts
        {kind: 'J'|'E', name, layer?, dz?}. z_lines: the mesh z lines (for the dual cell height). """
    (x0, x1), (y0, y1) = board.BoundingBox()[:2]
    fr = [float(f) for f in freqs]
    out = []
    zl = np.unique(np.asarray(z_lines, float)) if z_lines is not None else None
    for c in layers or board.copper_layers:
        z = board.SheetZ(c)
        name = 'J_' + c.replace('.', '_')
        d = CSX.AddDump(name, dump_type=13, dump_mode=0, file_type=1, frequency=fr)
        d.AddBox([x0 - margin, y0 - margin, z], [x1 + margin, y1 + margin, z])
        dz = None
        if zl is not None and len(zl) > 2:
            i = int(np.argmin(np.abs(zl - z)))
            lo, hi = zl[max(i - 1, 0)], zl[min(i + 1, len(zl) - 1)]
            dz = float(hi - lo) / 2.0 if hi > lo else None
        out.append(dict(kind='J', layer=c, name=name, z=float(z), dz=dz))
    if volume:
        zs = [board.SheetZ(c) for c in board.copper_layers]
        t = max(zs) - min(zs)
        air = min(1.0, max(0.3, 0.5 * t))
        diag = math.hypot(x1 - x0, y1 - y0)
        r = max(0.05, diag / 260.0)                      # ~250 points along the board diagonal
        rz = max(0.01, t / 40.0) if t > 0 else 0.05
        d = CSX.AddDump('E_vol', dump_type=10, dump_mode=1, file_type=1, frequency=fr, opt_resolution=[r, r, rz])
        d.AddBox([x0 - margin, y0 - margin, min(zs) - air], [x1 + margin, y1 + margin, max(zs) + air])
        out.append(dict(kind='E', name='E_vol', res=[r, r, rz]))
    return out


def _read(path):
    """ (x, y, z, F, A) of an FD dump: mesh lines in drawing units and A[f, comp, nz, ny, nx] complex. """
    import h5py
    try:
        from openEMS.utilities import HDF5Dump            # the fork's reader handles every layout
        with HDF5Dump(path) as d:
            mesh = d.GetMesh()
            x, y, z = (np.asarray(l, float) / (mesh.get('scaling') or 1.0) for l in mesh['lines'][:3])
            F = np.asarray(d.frequencies, float)
            A = []
            for k in range(len(F)):
                v = np.asarray(d.GetFieldAtIndex(f_idx=k))            # (3, nx, ny, nz)
                if v.ndim == 3:
                    v = v[..., None]
                A.append(np.transpose(v, (0, 3, 2, 1)))                # -> (3, nz, ny, nx)
            return x, y, z, F, np.asarray(A)
    except (ImportError, AttributeError):
        # no fork bindings, or the official openEMS package whose HDF5Dump has no `frequencies`
        pass
    with h5py.File(path, 'r') as h:
        x, y = np.asarray(h['Mesh/x'], float), np.asarray(h['Mesh/y'], float)
        z = np.asarray(h['Mesh/z'], float) if 'Mesh/z' in h else np.zeros(1)
        fd = h['FieldData/FD']
        F = np.asarray(fd.attrs.get('frequency', []), float)
        A = []
        k = 0
        while 'f%d_real' % k in fd:                                  # legacy layout: f<k>_real (3, nz, ny, nx)
            a = np.asarray(fd['f%d_real' % k]) + 1j * np.asarray(fd['f%d_imag' % k])
            A.append(a if a.ndim == 4 else a[:, None])
            k += 1
        while 'f%d' % k in fd:                                       # HDF5 v0.3: complex (3, nx, ny, nz), d_order NXYZ
            a = np.asarray(fd['f%d' % k])
            order = fd['f%d' % k].attrs.get('d_order', 'NXYZ')
            order = order.decode() if isinstance(order, bytes) else str(order)
            if a.ndim == 3:
                a = a[..., None] if order == 'NXYZ' else a[:, None]
            if order == 'NXYZ':
                a = np.transpose(a, (0, 3, 2, 1))
            A.append(a)
            k += 1
        if not A:
            raise ValueError('no field data in %s (keys: %s)' % (path, list(fd.keys())[:6]))
        if not len(F):
            F = np.array([float(fd['f%d' % i].attrs.get('frequency', 0)) for i in range(len(A))]) if 'f0' in fd else np.zeros(len(A))
        return x, y, z, F, np.asarray(A)


# ---------------------------------------------------------------- grids
def _lin_w(src, dst):
    """ indices and weights of a linear interpolation from sorted src points to dst points (clamped) """
    src = np.asarray(src, float)
    if len(src) == 1:
        return np.zeros(len(dst), int), np.zeros(len(dst), int), np.zeros(len(dst))
    i = np.clip(np.searchsorted(src, dst) - 1, 0, len(src) - 2)
    w = np.clip((dst - src[i]) / np.maximum(src[i + 1] - src[i], 1e-12), 0.0, 1.0)
    return i, i + 1, w


def _resample2(a, x, y, xs, ys):
    """ bilinear resampling of a[..., ny, nx] from axes (x, y) to (xs, ys) """
    i0, i1, wx = _lin_w(x, xs)
    j0, j1, wy = _lin_w(y, ys)
    a0 = a[..., j0, :] * (1 - wy)[:, None] + a[..., j1, :] * wy[:, None]
    return a0[..., i0] * (1 - wx) + a0[..., i1] * wx


def j_grid(board, max_px=900, min_px_mm=0.02, margin=0.3):
    """ uniform export grid over the board: dict(x0, y0, dx, nx, ny) with square pixels """
    (x0, x1), (y0, y1) = board.BoundingBox()[:2]
    x0, y0, x1, y1 = x0 - margin, y0 - margin, x1 + margin, y1 + margin
    dx = max(min_px_mm, max(x1 - x0, y1 - y0) / max_px)
    nx, ny = int(math.ceil((x1 - x0) / dx)) + 1, int(math.ceil((y1 - y0) / dx)) + 1
    return dict(x0=float(x0), y0=float(y0), dx=float(dx), nx=nx, ny=ny)


def copper_mask(board, layer, g, grow=1):
    """ uint8 [ny, nx]: 1 where the layer has copper (grown by `grow` pixels) """
    from matplotlib.path import Path
    nx, ny, dx = g['nx'], g['ny'], g['dx']
    xs = g['x0'] + dx * np.arange(nx)
    ys = g['y0'] + dx * np.arange(ny)
    m = np.zeros((ny, nx), bool)
    for s in board.ShapesOn(layer):
        p = np.asarray(s.pts, float)
        if len(p) < 3:
            continue
        i0, i1 = np.searchsorted(xs, [p[:, 0].min(), p[:, 0].max()])
        j0, j1 = np.searchsorted(ys, [p[:, 1].min(), p[:, 1].max()])
        i0, j0 = max(i0 - 1, 0), max(j0 - 1, 0)
        i1, j1 = min(i1 + 1, nx), min(j1 + 1, ny)
        if i1 <= i0 or j1 <= j0:
            continue
        X, Y = np.meshgrid(xs[i0:i1], ys[j0:j1])
        pts = np.column_stack([X.ravel(), Y.ravel()])
        inside = Path(p).contains_points(pts, radius=1e-9)
        for h in (s.holes or []):
            if len(h) >= 3:
                inside &= ~Path(np.asarray(h, float)).contains_points(pts)
        m[j0:j1, i0:i1] |= inside.reshape(X.shape)
    return _grow(m, grow).astype(np.uint8)


def _grow(m, n):
    """ binary dilation by n pixels (4-neighbour steps), along the last two axes """
    m = np.asarray(m, bool)
    for _ in range(n):
        m = m | np.roll(m, 1, -2) | np.roll(m, -1, -2) | np.roll(m, 1, -1) | np.roll(m, -1, -1)
    return m


def _f16(a, path):
    """ write a float array as float16 with a scale so the largest value is ~1000; returns the scale """
    a = np.asarray(a, float)
    mx = float(np.nanmax(np.abs(a))) if a.size else 0.0
    scale = mx / 1000.0 if mx > 0 else 1.0
    (np.nan_to_num(a / scale)).astype('<f2').tofile(path)
    return scale


def _ri(c):
    """ complex [..., k] -> float [..., 2k] interleaved re, im per component """
    c = np.asarray(c)
    out = np.empty(c.shape[:-1] + (2 * c.shape[-1],), float)
    out[..., 0::2], out[..., 1::2] = c.real, c.imag
    return out


# ---------------------------------------------------------------- export
def export(sim_dir, out_dir, board, dumps, excite=1, norm=None, index=None, db_range=40.0, png=True):
    """ Export the dumps of one excitation. norm(f_Hz) -> complex factor that scales the raw fields to
        1 W incident at the driven port (None: raw). index: the dict returned for an earlier excitation of
        the same run (extended in place). Returns the index for result.json["fields"]. """
    os.makedirs(out_dir, exist_ok=True)
    rel = os.path.basename(out_dir.rstrip('/\\'))
    dumps = [d if isinstance(d, dict) else dict(kind='J', layer=d[0], name=d[1]) for d in dumps]
    idx = index if index is not None else dict(version=VERSION, power_W=POWER_W, normalised=norm is not None, excitations=[],
                                               db_range=db_range, maps={}, peak_db={})
    e = int(excite)
    tag = 'e%d' % e
    # ---- surface currents
    jd = [d for d in dumps if d['kind'] == 'J' and os.path.exists(os.path.join(sim_dir, d['name'] + '.h5'))]
    F = None
    if jd:
        g = idx.get('grid') or j_grid(board)
        xs = g['x0'] + g['dx'] * np.arange(g['nx'])
        ys = g['y0'] + g['dx'] * np.arange(g['ny'])
        layers = [d['layer'] for d in jd]
        if 'grid' not in idx:
            idx['grid'] = g
            idx['bbox'] = [g['x0'], g['y0'], g['x0'] + g['dx'] * (g['nx'] - 1), g['y0'] + g['dx'] * (g['ny'] - 1)]
            idx['layers'] = [dict(name=d['layer'], z=d.get('z')) for d in jd]
            mask = np.stack([copper_mask(board, l, g) for l in layers])
            mask.tofile(os.path.join(out_dir, 'mask.u8'))
            idx['mask'] = rel + '/mask.u8'
        mask = np.fromfile(os.path.join(out_dir, 'mask.u8'), np.uint8).reshape(len(layers), g['ny'], g['nx'])
        # The dump holds the current of each mesh cell at its node; a node just outside a copper edge carries the
        # edge current of the cell it shares with the metal. Zeroing at the geometric edge lost 30 % of the strip
        # current (coupon: 0.183 A instead of 0.20 A for 1 W): only copper-free areas a cell away are zeroed.
        keep = _grow(mask > 0, max(1, int(math.ceil(0.25 / g['dx']))))
        Js = []
        for d in jd:
            x, y, z, F, A = _read(os.path.join(sim_dir, d['name'] + '.h5'))
            if np.ptp(x) < (xs[-1] - xs[0]) / 100:              # mesh stored in metres
                x, y = x * 1e3, y * 1e3
            dz = (d.get('dz') or 0.0) * 1e-3                     # rot(H) [A/m^2] * dual cell height [m] = Js [A/m]
            a = A[:, :2, 0] * (dz if dz > 0 else 1.0)            # (nf, 2, ny, nx)
            Js.append(_resample2(a, x, y, xs, ys))               # (nf, 2, ny_g, nx_g)
        Js = np.stack(Js, 1)                                     # (nf, nl, 2, ny, nx)
        nf = Js.shape[0]
        if 'f_GHz' not in idx:
            idx['f_GHz'] = [round(float(f) / 1e9, 4) for f in F]
        files, scales, peaks = {}, {}, {}
        for k in range(nf):
            s = norm(F[k]) if norm else 1.0
            jk = Js[k] * s * keep[:, None]                        # zero off the copper: compresses, no noise arrows
            fn = 'J_%s_f%d.f16' % (tag, k)
            scales[str(k)] = _f16(_ri(np.moveaxis(jk, 1, -1)), os.path.join(out_dir, fn))
            files[str(k)] = rel + '/' + fn
            peaks[str(k)] = float(np.sqrt((np.abs(jk) ** 2).sum(1)).max())
        idx.setdefault('J', dict(unit='A/m', files={}, scale={}, peak={}))
        idx['J']['files'][tag], idx['J']['scale'][tag], idx['J']['peak'][tag] = files, scales, peaks
        if png and not idx['maps']:                              # legacy heat maps (plugin 0.2): first excitation
            mag = np.sqrt((np.abs(Js) ** 2).sum(2))              # (nf, nl, ny, nx)
            for k in range(nf):
                ref = float(mag[k].max()) or 1.0
                for li, layer in enumerate(layers):
                    db = 20 * np.log10(np.maximum(mag[k, li][::-1], 1e-30) / ref)    # row 0 = top = max y
                    fn = '%s_%dMHz.png' % (layer.replace('.', '_'), round(F[k] / 1e6))
                    write_png(os.path.join(out_dir, fn), colormap((db + db_range) / db_range))
                    idx['maps'].setdefault(layer, []).append(rel + '/' + fn)
                    idx['peak_db'].setdefault(layer, []).append(round(float(db.max()), 1))
    # ---- E volume and layer voltages
    ed = [d for d in dumps if d['kind'] == 'E' and os.path.exists(os.path.join(sim_dir, d['name'] + '.h5'))]
    if ed:
        x, y, z, F2, A = _read(os.path.join(sim_dir, ed[0]['name'] + '.h5'))
        sc = 1e3 if np.ptp(x) < 1.0 else 1.0                     # metres -> mm
        x, y, z = x * sc, y * sc, z * sc
        if F is None:
            F = F2
            idx.setdefault('f_GHz', [round(float(f) / 1e9, 4) for f in F])
        if 'E' not in idx:
            idx['E'] = dict(unit='V/m', x=[round(float(v), 5) for v in x], y=[round(float(v), 5) for v in y],
                            z=[round(float(v), 5) for v in z], files={}, scale={}, peak={})
            cu = list(board.copper_layers)
            zs = {c: float(board.SheetZ(c)) for c in cu}
            pairs = []
            for i in range(len(cu) - 1):
                pairs.append([cu[i], cu[i + 1]])                 # V of the upper layer against the lower one
            idx['V'] = dict(unit='V', pairs=pairs, z={c: zs[c] for c in cu}, files={}, scale={}, peak={})
        files, scales, peaks = {}, {}, {}
        vfiles, vscales, vpeaks = {}, {}, {}
        zs = idx['V']['z']
        for k in range(A.shape[0]):
            s = norm(F2[k]) if norm else 1.0
            ek = A[k] * s                                         # (3, nz, ny, nx) V/m
            fn = 'E_%s_f%d.f16' % (tag, k)
            scales[str(k)] = _f16(_ri(np.moveaxis(ek, 0, -1)), os.path.join(out_dir, fn))
            files[str(k)] = rel + '/' + fn
            peaks[str(k)] = float(np.sqrt((np.abs(ek) ** 2).sum(0)).max())
            vs = []
            for up, lo in idx['V']['pairs']:
                vs.append(layer_voltage(ek[2], z, zs[lo], zs[up]))
            if vs:
                vk = np.stack(vs)[..., None]                      # (npairs, ny, nx, 1)
                fn = 'V_%s_f%d.f16' % (tag, k)
                vscales[str(k)] = _f16(_ri(vk), os.path.join(out_dir, fn))
                vfiles[str(k)] = rel + '/' + fn
                vpeaks[str(k)] = float(np.abs(vk).max())
        idx['E']['files'][tag], idx['E']['scale'][tag], idx['E']['peak'][tag] = files, scales, peaks
        idx['V']['files'][tag], idx['V']['scale'][tag], idx['V']['peak'][tag] = vfiles, vscales, vpeaks
    if e not in idx['excitations']:
        idx['excitations'].append(e)
    return idx


def layer_voltage(Ez, z, z_lo, z_hi):
    """ V(x, y) = phi(z_hi) - phi(z_lo) = -integral_{z_lo}^{z_hi} Ez dz (trapezoid on the dump z lines,
        the interval clipped to the dump). Ez: complex [nz, ny, nx]. """
    z = np.asarray(z, float)
    eps = 1e-4 * max(z_hi - z_lo, 1e-6)
    inner = np.where((z > z_lo + eps) & (z < z_hi - eps))[0]
    if len(inner) == 0:                                   # no dump line inside: Ez at the middle
        i0, i1, w = _lin_w(z, np.array([(z_lo + z_hi) / 2]))
        return -(Ez[i0[0]] * (1 - w[0]) + Ez[i1[0]] * w[0]) * (z_hi - z_lo) * 1e-3
    # The node-interpolated Ez on a copper sheet averages the field below and above the metal (about half of
    # it): the end values come from the nearest inner line instead (coupon: -15 % voltage otherwise).
    zz = np.concatenate([[z_lo], z[inner], [z_hi]])
    e = np.concatenate([Ez[inner[:1]], Ez[inner], Ez[inner[-1:]]])
    trap = getattr(np, 'trapezoid', None) or np.trapz
    return -trap(e, zz, axis=0) * 1e-3


def load(run_dir, idx, quantity, excite, k):
    """ The exported arrays of one run: quantity 'J' -> complex [nl, ny, nx, 2] (A/m), 'E' -> complex
        [nz, ny, nx, 3] (V/m), 'V' -> complex [npairs, ny, nx] (V). """
    q = idx[quantity]
    tag, ks = 'e%d' % int(excite), str(int(k))
    path = os.path.join(run_dir, q['files'][tag][ks])
    a = np.fromfile(path, '<f2').astype(float) * float(q['scale'][tag][ks])
    if quantity == 'J':
        g = idx['grid']
        a = a.reshape(len(idx['layers']), g['ny'], g['nx'], 4)
        return a[..., 0::2] + 1j * a[..., 1::2]
    if quantity == 'E':
        a = a.reshape(len(q['z']), len(q['y']), len(q['x']), 6)
        return a[..., 0::2] + 1j * a[..., 1::2]
    E = idx['E']
    a = a.reshape(len(q['pairs']), len(E['y']), len(E['x']), 2)
    return a[..., 0] + 1j * a[..., 1]


# ---------------------------------------------------------------- polylines for overlays
def instant(c, phase_deg):
    """ instantaneous value Re(c e^{j phase}) of a phasor array """
    return (np.asarray(c) * np.exp(1j * math.radians(phase_deg))).real


def peak_phase(c):
    """ phase (deg) at which sum |Re(c e^{j phase})|^2 is largest: the instant of the strongest field """
    c = np.asarray(c).ravel()
    s2 = np.sum(c * c)                        # sum |Re(c e^jp)|^2 = (sum|c|^2 + Re(s2 e^{2jp})) / 2
    return float(math.degrees(-np.angle(s2) / 2.0)) % 360.0


def streamlines(vx, vy, mask, x0, y0, dx, dsep_px=8.0, dtest=0.5, min_rel=0.03, max_len=4000, step=0.5):
    """ Evenly spaced streamlines (Jobard-Lefer) of a 2D vector field on a uniform grid.
        vx, vy: [ny, nx] real; mask: where lines may run. Seeds are tried from the strongest field down,
        a line stops where the field is below min_rel of the peak, leaves the mask or comes closer than
        dtest*dsep to another line. Returns [{pts: [[x, y]...], mag: [...]}] in board units. """
    ny, nx = vx.shape
    mag = np.hypot(vx, vy)
    mx = float(mag.max())
    if mx <= 0:
        return []
    ok = (mask > 0) & (mag > min_rel * mx)
    cs = max(dsep_px * dtest, 1.0)                 # occupancy cells
    gw, gh = int(nx / cs) + 1, int(ny / cs) + 1
    occ = -np.ones((gh, gw), int)

    def field(px, py):
        if px < 0 or py < 0 or px > nx - 1.001 or py > ny - 1.001:
            return None
        i, j = int(px), int(py)
        fx, fy = px - i, py - j
        if not (ok[j, i] or ok[min(j + 1, ny - 1), min(i + 1, nx - 1)]):
            return None
        u = (vx[j, i] * (1 - fx) + vx[j, i + 1] * fx) * (1 - fy) + (vx[j + 1, i] * (1 - fx) + vx[j + 1, i + 1] * fx) * fy
        v = (vy[j, i] * (1 - fx) + vy[j, i + 1] * fx) * (1 - fy) + (vy[j + 1, i] * (1 - fx) + vy[j + 1, i + 1] * fx) * fy
        m = math.hypot(u, v)
        if m < min_rel * mx:
            return None
        return u / m, v / m, m

    def free(px, py, line_id, d):
        ci, cj = int(px / cs), int(py / cs)
        r = int(math.ceil(d / cs))
        for jj in range(max(cj - r, 0), min(cj + r + 1, gh)):
            for ii in range(max(ci - r, 0), min(ci + r + 1, gw)):
                o = occ[jj, ii]
                if o >= 0 and o != line_id:
                    return False
        return True

    def trace(px, py, sgn, line_id):
        pts, mags = [], []
        for _ in range(max_len):
            a = field(px, py)
            if a is None:
                break
            hx, hy = px + sgn * step * a[0], py + sgn * step * a[1]
            b = field(hx, hy)
            if b is None:
                break
            px, py = px + sgn * step * (a[0] + b[0]), py + sgn * step * (a[1] + b[1])
            if not free(px, py, line_id, dsep_px * dtest):
                break
            if pts and math.hypot(px - pts[0][0], py - pts[0][1]) < step and len(pts) > 8:
                break                                          # closed loop
            pts.append((px, py))
            mags.append(a[2])
        return pts, mags

    order = np.argsort(-np.where(ok, mag, 0).ravel())
    out = []
    stride = max(int(dsep_px / 2), 1)
    for flat in order[::stride]:
        j, i = divmod(int(flat), nx)
        if not ok[j, i] or not free(i, j, -2, dsep_px):
            continue
        lid = len(out)
        fw, fm = trace(float(i), float(j), 1, lid)
        bw, bm = trace(float(i), float(j), -1, lid)
        pts = bw[::-1] + [(float(i), float(j))] + fw
        ms = bm[::-1] + [float(mag[j, i])] + fm
        if len(pts) < 4:
            continue
        for (px, py) in pts:
            occ[min(int(py / cs), gh - 1), min(int(px / cs), gw - 1)] = lid
        out.append(dict(pts=[[x0 + px * dx, y0 + py * dx] for px, py in pts], mag=[float(m) for m in ms]))
        if len(out) > 3000:
            break
    return out


def contour_lines(Z, x, y, levels):
    """ iso-lines of Z[ny, nx] on axes x, y: [{level, pts}] (contourpy, which comes with matplotlib) """
    import contourpy
    gen = contourpy.contour_generator(np.asarray(x, float), np.asarray(y, float), np.asarray(Z, float))
    out = []
    for lv in levels:
        for seg in gen.lines(float(lv)):
            if len(seg) >= 2:
                out.append(dict(level=float(lv), pts=np.asarray(seg, float).round(4).tolist()))
    return out


# ---------------------------------------------------------------- images
# a perceptually ordered dark-to-bright map (inferno-like anchors)
_ANCHORS = np.array([[0, 0, 4], [40, 11, 84], [101, 21, 110], [159, 42, 99], [212, 72, 66], [245, 125, 21],
                     [250, 193, 39], [252, 255, 164]], float)
# diverging map for signed quantities (blue - white - red)
_DIV = np.array([[33, 102, 172], [103, 169, 207], [209, 229, 240], [247, 247, 247], [253, 219, 199], [239, 138, 98],
                 [178, 24, 43]], float)


def colormap(v, anchors=None):
    """ v in [0, 1] -> uint8 RGB. """
    an = _ANCHORS if anchors is None else anchors
    v = np.clip(np.asarray(v, float), 0, 1) * (len(an) - 1)
    i = np.minimum(v.astype(int), len(an) - 2)
    t = (v - i)[..., None]
    return (an[i] * (1 - t) + an[i + 1] * t).astype(np.uint8)


def write_png(path, rgb):
    """ Minimal RGB / RGBA PNG writer (uint8 [h, w, 3|4], row 0 = top). """
    h, w, ch = rgb.shape
    raw = b''.join(b'\x00' + rgb[r].tobytes() for r in range(h))

    def chunk(t, d):
        c = struct.pack('>I', len(d)) + t + d
        return c + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
    with open(path, 'wb') as fh:
        fh.write(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 6 if ch == 4 else 2, 0, 0, 0))
                 + chunk(b'IDAT', zlib.compress(raw, 6)) + chunk(b'IEND', b''))
