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
launch_study: edge-SMA launch variants on the JLC04161H-7628 stackup.

A 20 mm 50 Ohm microstrip with the same launch at both ends. The connector is modelled
as its centre pin (a metal cylinder lying on the signal pad) and a short 50 Ohm PTFE coax
section outside the board edge whose flange touches the top and bottom GND; the lumped
ports sit at the far end of the coax sections. The board itself (pads, GND pads, vias,
planes, the In1 cut-out under the pad) comes from a generated .kicad_pcb through
CSXCAD.KiCadImport, so the same footprint can go straight onto the test panel.

    python launch_study.py out_dir [--only NAME]
"""

import json
import math
import os
import sys
import uuid

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kicad_sim as ks                                       # noqa: E402
from make_coupon import STACKUP                              # noqa: E402

C0 = 299792458.0
W50 = 0.3726            # 50 Ohm microstrip on 7628 (Hammerstad-Jensen with t = 35 um)
H = 14.0                # board width (y)
L = 20.0                # board length (x)


def U():
    return str(uuid.uuid4())


# ---------------------------------------------------------------- launch variants
# pin: centre-pin diameter; wp/lp: signal pad width/length; g: gap to the top GND pads;
# cut: In1 cut-out under the pad (length from the edge, width), taper: pad-to-line taper length
VARIANTS = {
    'thick_plain':    dict(pin=1.27, wp=1.30, lp=3.0, g=0.50, cut=None, taper=0.0),
    'thick_cut':      dict(pin=1.27, wp=1.30, lp=3.0, g=0.50, cut=(3.0, 2.3), taper=0.0),
    'thick_cut_tap':  dict(pin=1.27, wp=1.30, lp=2.5, g=0.50, cut=(3.5, 2.3), taper=1.0),
    'thick_narrow':   dict(pin=1.27, wp=0.80, lp=3.0, g=0.40, cut=(3.0, 1.6), taper=0.6),
    'thin_plain':     dict(pin=0.76, wp=0.80, lp=2.0, g=0.40, cut=None, taper=0.5),
}


def board_text(v):
    """ KiCad board: 20 mm line, launch footprint at both ends, GND on In1/In2/B.Cu, stitching vias. """
    y = H / 2
    wp, lp, g, tap = v['wp'], v['lp'], v['g'], v['taper']
    fps, vias, segs = '', '', ''
    for ref, x0, sgn in (('J1', 0.0, 1), ('J2', L, -1)):
        cx = x0 + sgn * lp / 2
        gw = 1.6                                       # top GND pad width (connector legs)
        gy = wp / 2 + g + gw / 2
        fps += '''  (footprint "z2labs:SMA_Edge_Launch" (layer "F.Cu") (uuid "%s") (at %g %g)
    (property "Reference" "%s" (at 0 -4 0) (layer "F.SilkS") (uuid "%s") (effects (font (size 0.8 0.8) (thickness 0.12))))
    (property "Value" "SMA edge pin %.2f" (at 0 4 0) (layer "F.Fab") (uuid "%s") (effects (font (size 0.8 0.8) (thickness 0.12))))
    (attr smd)
    (pad "1" smd rect (at 0 0) (size %g %g) (layers "F.Cu" "F.Mask") (net 2 "RF") (uuid "%s"))
    (pad "2" smd rect (at 0 %g) (size %g %g) (layers "F.Cu" "F.Mask") (net 1 "GND") (uuid "%s"))
    (pad "2" smd rect (at 0 %g) (size %g %g) (layers "F.Cu" "F.Mask") (net 1 "GND") (uuid "%s"))
    (pad "2" smd rect (at 0 0) (size %g %g) (layers "B.Cu" "B.Mask") (net 1 "GND") (uuid "%s"))
  )
''' % (U(), cx, y, ref, U(), v['pin'], U(), lp, wp, U(), -gy, lp, gw, U(), gy, lp, gw, U(), lp, wp + 2 * (g + gw), U())
        for k in range(3):                             # GND vias in the legs, next to the pad
            vx = x0 + sgn * (0.6 + k * 1.0)
            for vy in (y - gy, y + gy):
                vias += '  (via (at %g %g) (size 0.6) (drill 0.3) (layers "F.Cu" "B.Cu") (net 1) (uuid "%s"))\n' % (vx, vy, U())
        # taper from the pad width to the line width (a short trapezoid as a graphic copper polygon)
        xa = x0 + sgn * lp
        xb = xa + sgn * tap
        if tap > 0:
            segs += '''  (gr_poly (pts (xy %g %g) (xy %g %g) (xy %g %g) (xy %g %g)) (stroke (width 0) (type solid)) (fill solid) (layer "F.Cu") (net 2) (uuid "%s"))
''' % (xa, y - wp / 2, xb, y - W50 / 2, xb, y + W50 / 2, xa, y + wp / 2, U())
    xs, xe = lp + v['taper'], L - lp - v['taper']
    segs += '  (segment (start %g %g) (end %g %g) (width %g) (layer "F.Cu") (net 2) (uuid "%s"))\n' % (xs, y, xe, y, W50, U())
    x = 4.0
    while x < L - 3.9:
        for yy in (y - 2.5, y + 2.5):
            vias += '  (via (at %g %g) (size 0.6) (drill 0.3) (layers "F.Cu" "B.Cu") (net 1) (uuid "%s"))\n' % (x, yy, U())
        x += 2.0
    # In1: full plane, with a notch from each board edge under the signal pad when the variant cuts it
    if v['cut']:
        cl, cw = v['cut']
        a, b = y - cw / 2, y + cw / 2
        pts = [(0, 0), (L, 0), (L, a), (L - cl, a), (L - cl, b), (L, b), (L, H), (0, H), (0, b), (cl, b), (cl, a), (0, a)]
    else:
        pts = [(0, 0), (L, 0), (L, H), (0, H)]
    zone_pts = ' '.join('(xy %g %g)' % p for p in pts)
    rect = '(xy 0 0) (xy %g 0) (xy %g %g) (xy 0 %g)' % (L, L, H, H)
    zones = ''
    for layer, pp in (('In1.Cu', zone_pts), ('In2.Cu', rect), ('B.Cu', rect)):
        zones += '''  (zone (net 1) (net_name "GND") (layer "%s") (uuid "%s") (hatch edge 0.5)
    (connect_pads (clearance 0.2)) (min_thickness 0.2) (filled_areas_thickness no) (fill yes)
    (polygon (pts %s))
    (filled_polygon (layer "%s") (pts %s)))
''' % (layer, U(), pp, layer, pp)
    return '''(kicad_pcb (version 20241229) (generator "pcbnew") (generator_version "9.0")
  (general (thickness 1.6) (legacy_teardrops no))
  (paper "A4")
  (layers (0 "F.Cu" signal) (4 "In1.Cu" signal) (6 "In2.Cu" signal) (2 "B.Cu" signal)
    (5 "F.SilkS" user) (7 "B.SilkS" user) (1 "F.Mask" user) (3 "B.Mask" user) (25 "Edge.Cuts" user)
    (31 "F.CrtYd" user) (35 "F.Fab" user))
  (setup
%s
    (pad_to_mask_clearance 0))
  (net 0 "") (net 1 "GND") (net 2 "RF")
%s%s  (gr_rect (start 0 0) (end %g %g) (stroke (width 0.05) (type default)) (fill no) (layer "Edge.Cuts") (uuid "%s"))
%s%s)
''' % (STACKUP, fps, segs, L, H, U(), vias, zones)


# ---------------------------------------------------------------- connector model
COAX_D = 2.6            # bore diameter of the strip-in-coax feed: 50 Ohm (coax_test: D 2.4 -> 43-49, 2.8 -> 51-67 Ohm)
COAX_LEN = 8.0          # coax from the board edge to the PML boundary
COAX_FEED = 3.0         # excitation plane from the boundary (outside the 2 mm PML)
COAX_MEAS = 5.5         # measurement plane from the boundary
BODY = 'shell_plate'    # outer-conductor shell + flange plate at the edge. A solid flange block running into the
                        # PML made the FDTD unstable whenever an air launch had GND vias beside the pad (not understood)


def connector(CSX, board, v, x_edge, sgn, yc=None, D=None, er_ptfe=2.05):
    """ Connector model: a 50 Ohm strip-in-coax feed outside the board edge (flat strip of the pin width as
        the inner conductor, PTFE, metal body that touches the board edge), the round pin lying on the pad from
        1 mm outside the edge to 0.3 mm before the pad end. The port is an MSL-type port on the strip (voltage
        strip-to-bore wall, current around the strip) with the excitation 3 mm and the measurement plane
        5.5 mm inside the boundary; the coax runs into an 8-cell PML on that side. A lumped port inside a coax
        does not work (it only loads the bottom of the annulus) and a port on a PEC end wall is unstable.
        Returns (port definition for make_ports, X, Y, Z mesh lines). """
    d = v['pin']
    D = D or v.get('D') or COAX_D * d / 1.27                     # the strip-in-coax feed scales with the pin (same geometry, same Z0)
    zt = board.SheetZ(board.copper_layers[0])
    zb = board.SheetZ(board.copper_layers[-1])
    fy, fz = v.get('flange') or (D / 2 + 1.5, D / 2 + 1.5)       # flange plate half-sizes (y, z) around the pin axis
    zc = zt + d / 2                                              # pin axis, bore axis
    yc = -H / 2 if yc is None else yc
    x_in = x_edge + sgn * (v.get('pin_len') or v['lp'] - 0.3)    # pin tip: the part's pin length, else 0.3 mm short of the pad end
    x_end = x_edge - sgn * COAX_LEN                              # PML boundary
    tag = 'L' if sgn > 0 else 'R'
    pin = CSX.AddMetal('pin_' + tag)
    pin.AddCylinder([x_edge - sgn * 1.0, yc, zc], [x_in, yc, zc], d / 2, priority=20)
    ptfe = CSX.AddMaterial('ptfe_' + tag, epsilon=er_ptfe)
    ptfe.AddCylinder([x_end, yc, zc], [x_edge, yc, zc], D / 2, priority=15)
    body = CSX.AddMetal('body_' + tag)                           # flange + outer conductor (priority below the PTFE)
    if BODY == 'shell':                                          # experiment: thin outer conductor only, no flange block
        body.AddCylindricalShell([x_end, yc, zc], [x_edge, yc, zc], D / 2 + 0.25, 0.5, priority=10)
    elif BODY == 'shell_plate':                                  # outer conductor shell + 0.5 mm flange plate at the edge
        body.AddCylindricalShell([x_end, yc, zc], [x_edge, yc, zc], D / 2 + 0.25, 0.5, priority=10)
        body.AddBox([x_edge - sgn * 0.5, yc - fy, zc - fz], [x_edge, yc + fy, zc + fz], priority=10)
        if v.get('legs'):                                        # ground legs clamping the board (y centre, width, length, thickness)
            ly, lw, ll, lt = v['legs']
            for sy in (-1, 1):
                body.AddBox([x_edge, yc + sy * ly - lw / 2, zt], [x_edge + sgn * ll, yc + sy * ly + lw / 2, zt + lt], priority=10)
                body.AddBox([x_edge, yc + sy * ly - lw / 2, zb - lt], [x_edge + sgn * ll, yc + sy * ly + lw / 2, zb], priority=10)
    elif BODY == 'box_gap':                                      # experiment: flange block 0.05 mm short of the board edge
        body.AddBox([x_end, yc - D / 2 - 1.5, zc - D / 2 - 1.5], [x_edge - sgn * 0.05, yc + D / 2 + 1.5, zc + D / 2 + 1.5], priority=10)
    elif BODY == 'box_top':                                      # experiment: flange block not below the board bottom
        body.AddBox([x_end, yc - D / 2 - 1.5, 0.0], [x_edge, yc + D / 2 + 1.5, zc + D / 2 + 1.5], priority=10)
    else:
        body.AddBox([x_end, yc - D / 2 - 1.5, zc - D / 2 - 1.5], [x_edge, yc + D / 2 + 1.5, zc + D / 2 + 1.5], priority=10)
    strip = CSX.AddMetal('strip_' + tag)
    port = dict(metal=strip, start=[x_end, yc - d / 2, zt], stop=[x_edge, yc + d / 2, zc - D / 2],
                feed=COAX_FEED, meas=COAX_MEAS, sgn=sgn)
    X = [x_end, x_end + sgn * COAX_FEED, x_end + sgn * COAX_MEAS, x_edge - sgn * 1.0, x_edge, x_in] + \
        list(x_end + sgn * np.arange(0.25, 2.6, 0.25))            # fine cells in the PML (8 cells = 2 mm)
    Y = [yc - fy, yc - D / 2, yc - d / 2, yc, yc + d / 2, yc + D / 2, yc + fy]
    Z = [zc - fz, zc - D / 2, zt, zc, zc + d / 2, zc + D / 2, zc + fz]
    if v.get('legs'):
        ly, lw, ll, lt = v['legs']
        Y += [yc + sy * ly + e for sy in (-1, 1) for e in (-lw / 2, lw / 2)]
        Z += [zt + lt, zb - lt]
        X.append(x_edge + sgn * ll)
    return port, X, Y, Z


def make_ports(FDTD, pdefs, excite=1):
    """ MSL-type ports on the coax strips (call after the mesh lines are set). """
    out = []
    for i, p in enumerate(pdefs):
        out.append(FDTD.AddMSLPort(i + 1, p['metal'], p['start'], p['stop'], 'x', 'z', excite=1 if i + 1 == excite else 0,
                                   FeedShift=p['feed'], MeasPlaneShift=p['meas'], priority=20))
    return out


BOUNDARY = ['PML_8'] * 6     # was PML x + MUR y/z: MUR next to the connector flange was unstable (2026-10-21)


def run_variant(name, v, out, fmax=10e9, variant=None, threads=4, log=print):
    ks._setup_python()
    from CSXCAD import ContinuousStructure
    from CSXCAD.KiCadImport import KiCadBoard
    from openEMS import openEMS
    os.makedirs(out, exist_ok=True)
    pcb = os.path.join(out, '%s.kicad_pcb' % name)
    open(pcb, 'w').write(board_text(v))
    board = KiCadBoard(pcb)
    fmin = 0.2e9
    f = np.linspace(fmin, fmax, 401)
    S = np.zeros((len(f), 2, 2), complex)
    meta = {}
    for e in (1, 2):
        er = 4.43
        t_max = max(3e-9, 15 * math.hypot(L, H) * 1e-3 * math.sqrt(er) / C0)
        FDTD = openEMS(EndCriteria=1e-4, NrTS=300000, MaxTime=t_max)
        FDTD.SetGaussExcite(0.5 * fmax + 0.15e9, 0.5 * fmax - 0.15e9)      # no DC: the floating strips would keep a static charge
        FDTD.SetBoundaryCond(BOUNDARY)
        CSX = ContinuousStructure()
        FDTD.SetCSX(CSX)
        board.ToCSX(CSX, frequency=0.5 * fmax, copper='sheet')
        # pseudo-port dicts so kicad_sim meshes the RF copper of the board
        pads = [ks.find_pad(board, r) for r in ('J1', 'J2')]
        fake = [dict(net='RF', ref_net='GND', start=[0, 0, 0], stop=[0, 0, 0], normal_axis=0)]
        X, Y, Z, res = ks.mesh_lines(board, fake, fmax, margin=4.0)
        X, Y, Z = list(X), list(Y), list(Z)
        pdefs = []
        for x_edge, sgn in ((0.0, 1), (L, -1)):
            pd, x2, y2, z2 = connector(CSX, board, v, x_edge, sgn)
            pdefs.append(pd)
            X += x2; Y += y2; Z += z2
        g = CSX.GetGrid()
        g.SetDeltaUnit(1e-3)
        for ax, lines in zip('xyz', (X, Y, Z)):
            g.SetLines(ax, np.unique(np.round(ks.smooth(lines, res, merge=0.05), 5)))
        ports = make_ports(FDTD, pdefs, excite=e)
        nl = [len(g.GetLines(a)) for a in 'xyz']
        cdir = os.path.join(out, name, 'exc%d' % e)
        os.makedirs(cdir, exist_ok=True)
        xml = os.path.join(cdir, 'model.xml')
        FDTD.Write2XML(xml)
        log('%s exc%d: mesh %s (%.2f MCells)' % (name, e, nl, np.prod(nl) / 1e6))
        r = ks.run_engine(xml, os.path.join(cdir, 'sim'), variant, threads)
        log('  rc=%s %.0f s, %s TS' % (r['rc'], r['wall_s'], r['timesteps']))
        for p in ports:
            p.CalcPort(os.path.join(cdir, 'sim'), f, ref_impedance=50)
        for i, p in enumerate(ports):
            S[:, i, e - 1] = p.uf_ref / ports[e - 1].uf_inc
        meta = dict(mesh=nl, mcells=float(np.prod(nl) / 1e6))
    ks.write_touchstone(os.path.join(out, '%s.s2p' % name), f, S, comment='launch_study %s %s' % (name, json.dumps(v)))
    s11 = 20 * np.log10(np.abs(S[:, 0, 0]) + 1e-15)
    s21 = 20 * np.log10(np.abs(S[:, 1, 0]) + 1e-15)
    res = dict(name=name, variant=v, f_GHz=(f / 1e9).round(4).tolist(), S11_dB=s11.round(2).tolist(), S21_dB=s21.round(3).tolist(),
               s11_max_6G=float(s11[f <= 6e9].max()), s11_max_10G=float(s11.max()), s21_min_6G=float(s21[f <= 6e9].min()), **meta)
    json.dump(res, open(os.path.join(out, '%s.json' % name), 'w'))
    log('  %s: S11 max %.1f dB (<6 GHz), %.1f dB (<10 GHz)' % (name, res['s11_max_6G'], res['s11_max_10G']))
    return res


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'launch_study_out')
    only = sys.argv[sys.argv.index('--only') + 1].split(',') if '--only' in sys.argv else None
    for name, v in VARIANTS.items():
        if only and name not in only:
            continue
        run_variant(name, v, out, log=lambda m: print(m, flush=True))
