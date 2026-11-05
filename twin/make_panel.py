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
make_panel: edge-SMA launch coupons on JLC04161H-7628 (4 layers, 1.6 mm), GCPW lines.

Coupons of 50 x 12 mm, every coupon a two-port with the same edge-SMA launch at both ends. Coupon GND on F.Cu
(coplanar), In1, In2 and B.Cu, via fences along the line. The launch footprint is parametrised by the mechanisms the
literature identifies (Southwest Microwave 2007/2011, Signal Integrity Journal "Demystifying edge launch connectors",
MDPI Electronics 11(13):1990 and 13(14):2686, Altium "SMA edge connector transitions", Wuerth ANE012, Copper Mountain
fixture note):

    wp, lp   signal pad width / length (the pad reaches the board edge)
    g        coplanar gap around the pad (very sensitive, ~1-2 Ohm per 10 um)
    cut      reference-plane cut-out under the pad: layers, extra length beyond the pad, extra width
    taper    pad-to-line taper length (version-2 taper of Southwest: never narrower than the line)
    vias     GND via pairs next to the pad (tie the coplanar ground to the planes, kill the pad resonance)
    pull     copper pull-back from the SMA edge (0 = flush; 0.3 = JLC default copper-to-edge rule)
    air      "air" launch: the board is routed away under the hanging pin, air = pad length at the pin tip

COUPONS holds the BWSMA-KE-P001 air launch A0 and two comparison launches for the same connector. Every coupon is
also written standalone and sim_coupon() runs it with the connector model of launch_study.py (centre pin on the pad,
50 Ohm PTFE coax stub, metal flange, lumped port at the coax end).

    python make_panel.py out_dir                  # panel + standalone coupons + table
    python make_panel.py out_dir --sim 1,29       # also simulate those coupons
"""

import json
import math
import os
import sys
import time
import uuid

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_coupon import STACKUP                              # noqa: E402

H_SUB, ER = 0.2104, 4.4             # outer prepreg (7628) and its permittivity
PITCH = 12.0                        # coupon height (V-cut pitch)
LC = 50.0                           # coupon length
VCUT_CLEAR = 0.4                    # copper to V-cut line (JLC)
W50 = 0.3726                        # 50 Ohm microstrip
GCPW = {0.15: 0.3555, 0.20: 0.3758, 0.25: 0.3889, 0.35: 0.4054, 0.50: 0.4199}   # gap -> width, 50 Ohm (msl_ref.gcpw)
S0 = 0.25                           # baseline GCPW gap
LPF_HI, LPF_LO = 0.15, 2.4
FDTD_KW = {}                        # extra openEMS() keywords (TimeStepMethod, TimeStepFactor) for experiments
if os.environ.get('Z2_TSF'):        # A/B switches of the stability study (2026-10-21)
    FDTD_KW['TimeStepFactor'] = float(os.environ['Z2_TSF'])
if os.environ.get('Z2_TSM'):
    FDTD_KW['TimeStepMethod'] = int(os.environ['Z2_TSM'])
MARGIN = float(os.environ.get('Z2_MARGIN', '4.0'))      # air around the board (mm); the PML takes the outer 8 cells
LPF_SECT = [('C', 1.305), ('L', 9.859), ('C', 4.224), ('L', 9.859), ('C', 1.305)]


def U():
    return str(uuid.uuid4())


def model_name(v):
    return 'sma_edge_p%.2f_l%.1f.wrl' % (v['pin'], v['lp'])


def write_models(b, out):
    from sma_model import write_wrl
    for pin, lp in b.models:
        write_wrl(os.path.join(out, 'sma_edge_p%.2f_l%.1f.wrl' % (pin, lp)), pin=pin, lp=lp)
    import shutil
    if b.step and os.path.exists(WURTH_STEP):
        shutil.copy(WURTH_STEP, os.path.join(out, WURTH['step']))
    if b.bwsma:
        for k in ('step', 'wrl'):
            src = os.path.join(HERE, 'models', BWSMA[k])
            if os.path.exists(src):
                shutil.copy(src, os.path.join(out, BWSMA[k]))


# ---------------------------------------------------------------- launch variants
def launch(pin=1.27, wp=1.3, lp=3.0, g=0.5, cut=('In1.Cu',), cut_dl=0.5, cut_dw=0.6, taper=1.0, vias=3, pull=0.0, air=0.0, via_dy=None, air_layers=None, air_w=None, fp='generic'):
    """ air > 0: the pin hangs over a routed notch (length lp - air from the edge, width air_w, default wp + 2g);
        the pad is only the last `air` mm (the overlap), so the launch sees air instead of the laminate and the pad
        capacitance. The coplanar ground runs along both notch walls to the board edge on every layer (copper
        0.25 mm from the milled edge), so the pin over the notch still has its return beside it.
        fp='wurth': the Wuerth WR-SMA 60312202114509 end launch (pin 0.76, 1.6 mm board) with its own footprint:
        signal pad 0.9 x 5.6 to the edge, the two ground legs on 1.65 x 5.6 pads at y = +-4.125 top and bottom, its
        STEP model; wp/lp are fixed by the footprint, g (gap to the coplanar ground), cut, taper, vias, air stay free. """
    if fp == 'wurth':
        pin, wp, lp = WURTH['pin'], WURTH['wp'], WURTH['lp']
    elif fp == 'bwsma':                                # the pin is fixed; the pad may run past the pin tip (lp >= pin_len + 0.1)
        pin, lp = BWSMA['pin'], max(lp, BWSMA['pin_len'] + 0.1)
        if air and air < lp - BWSMA['pin_len'] + 0.1:  # the pad starts at lp - air: the pin tip must overlap it
            print('WARNING bwsma air launch: ov %.2f leaves the pin tip %.2f mm short of the pad (open circuit)' % (air, lp - air - BWSMA['pin_len']), file=sys.stderr)
    return dict(pin=pin, wp=wp, lp=lp, g=g, cut=tuple(cut) if cut else (), cut_dl=cut_dl, cut_dw=cut_dw,
                taper=taper, vias=vias, pull=pull, air=air, via_dy=via_dy, air_layers=air_layers, air_w=air_w, fp=fp)


# Wuerth WR-SMA 60312202114509 (datasheet + KiCad library): pin 0.76 mm, pads to the edge 5.6 long, flange 9.52 x 7.92,
# ground legs 1.6 wide at +-4.125 clamping the 1.6 mm board; VSWR 1.2 to 12.4 GHz. The STEP sits in twin/models.
WURTH = dict(pin=0.76, wp=0.9, lp=5.6, leg_y=4.125, leg_w=1.65, flange=(4.76, 3.96),
             step='J_Wurth_WR-SMA_60312202114509.step')
WURTH_STEP = os.path.join(HERE, 'models', WURTH['step'])

# BAT Wireless BWSMA-KE-P001 (LCSC C496550, datasheet July-2-2025-A01): straight edge SMA with an offset pin, the
# cheap LCSC part the panel is built for. Pin d 0.95 lying on the top copper, 4 mm beyond the flange; four 0.9 x 0.9
# legs 4 mm long at y = +-2.55, two in the pin plane on the top copper, two on the bottom copper, 1.7 mm slot
# (1.6 mm board); flange 6.5 x 6.5 x 1.5 centred on the pin axis. Model from bwsma_model.py (twin/models).
BWSMA = dict(pin=0.95, pin_len=4.0, leg_y=2.55, leg=0.9, leg_len=4.0, leg_pad=(4.6, 1.3), flange=(3.25, 3.25),
             step='BWSMA-KE-P001.step', wrl='BWSMA-KE-P001.wrl')

T0 = launch()
N0 = launch(pin=0.76, wp=0.8, lp=2.0, g=0.35, taper=0.5)
W0 = launch(fp='wurth', g=0.5, cut=(), taper=0.5, vias=3)
# BWSMA baseline B0: the T0 compensation scaled to the 0.95 pin (pad 1.0 x 4.3, In1 cut); BN = the naive footprint
# a datasheet-less layout gets (wide pad, small gap, no cut, no taper), the bad reference
B0 = launch(fp='bwsma', wp=1.0, g=0.5, cut=('In1.Cu',), cut_dl=0.5, cut_dw=0.6, taper=0.5, vias=3)


def B(**kw):
    """ A one-factor variant of B0. """
    p = dict(B0)
    p.update(kw)
    p.pop('pin')
    return launch(**{k: p[k] for k in ('wp', 'lp', 'g', 'cut', 'cut_dl', 'cut_dw', 'taper', 'vias', 'pull', 'air', 'via_dy', 'air_layers', 'air_w', 'fp')})


BN = launch(fp='bwsma', wp=1.5, g=0.3, cut=(), taper=0.0, vias=1)
# A0: the air launch found by a genetic search, with taper 0 (the best of the staircase sweep):
# the pin hangs over a 2.8 x 1.9 mm routed notch, pad 0.9 x 1.55 at the tip, one via pair, no plane cut
A0 = launch(fp='bwsma', wp=0.9, lp=4.35, g=0.5, cut=(), taper=0.0, vias=1, air=1.55)


def A(**kw):
    """ A one-factor variant of A0. """
    p = dict(A0)
    p.update(kw)
    p.pop('pin')
    return launch(**{k: p[k] for k in ('wp', 'lp', 'g', 'cut', 'cut_dl', 'cut_dw', 'taper', 'vias', 'pull', 'air', 'via_dy', 'air_layers', 'air_w', 'fp')})


# (id, label, launch, line kind, line params): the BWSMA-KE-P001 air launch A0 and two comparison launches for the
# same connector (the older In1-cut footprint and the naive footprint a datasheet-less layout gets)
COUPONS = [
    (1, 'A0 air ov1.55', A0, 'gcpw', {}),
    (28, 'B0 old In1 cut', B0, 'gcpw', {}),
    (29, 'BN naive', BN, 'gcpw', {}),
]


class Board:
    """ Collects items in absolute coordinates (KiCad, y down); GND zones are coupon rectangles minus holes. """

    def __init__(self):
        self.nets = ['GND']
        self.items, self.coupons, self.edge, self.user = [], [], [], []
        self.holes = {l: [] for l in ('F.Cu', 'In1.Cu', 'In2.Cu', 'B.Cu')}
        self.notches = []           # (x_edge, sgn, y_lo, y_hi, depth): routed notches in the board edge
        self.models = set()         # (pin, lp) of the 3D models used
        self.step = False           # a Wuerth footprint is used: copy its STEP next to the board
        self.bwsma = False          # a BWSMA-KE-P001 footprint is used: copy its STEP/WRL next to the board
        self.zone_lim = {}          # coupon index -> {layer: [(x0, x1), ...]} where the zone exists (default: whole)

    def net(self, name):
        if name not in self.nets:
            self.nets.append(name)
        return self.nets.index(name) + 1

    def seg(self, x0, y0, x1, y1, w, net, layer='F.Cu'):
        self.items.append('  (segment (start %g %g) (end %g %g) (width %g) (layer "%s") (net %d) (uuid "%s"))'
                          % (x0, y0, x1, y1, w, layer, self.net(net), U()))

    def via(self, x, y, net='GND', size=0.6, drill=0.3):
        # GND vias of the launch and of the fences can land on each other (short pads, lp 2): keep the first
        # one, skip a later GND via whose hole would sit closer than 0.5 mm to an existing one (JLC hole-to-hole)
        if net == 'GND':
            gv = self.__dict__.setdefault('_gnd_vias', [])
            if any((x - a) ** 2 + (y - b) ** 2 < (drill + 0.5) ** 2 for a, b in gv):
                return
            gv.append((x, y))
        self.items.append('  (via (at %g %g) (size %g) (drill %g) (layers "F.Cu" "B.Cu") (net %d) (uuid "%s"))'
                          % (x, y, size, drill, self.net(net), U()))

    def poly(self, pts, net, layer='F.Cu'):
        self.items.append('  (gr_poly (pts %s) (stroke (width 0) (type solid)) (fill solid) (layer "%s") (net %d) (uuid "%s"))'
                          % (' '.join('(xy %g %g)' % p for p in pts), layer, self.net(net), U()))

    def rect_cu(self, x0, y0, x1, y1, net, layer='F.Cu'):
        self.poly([(x0, y0), (x1, y0), (x1, y1), (x0, y1)], net, layer)

    def text(self, s, x, y, size=0.8, layer='F.SilkS'):
        self.items.append('  (gr_text "%s" (at %g %g 0) (layer "%s") (uuid "%s") (effects (font (size %g %g) (thickness %g)) (justify left)))'
                          % (s, x, y, layer, U(), size, size, size * 0.15))

    def hole(self, layers, x0, y0, x1, y1):
        for l in layers:
            self.holes[l].append((min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)))

    # ------------------------------------------------------------ launch footprint
    def launch(self, ref, x_edge, sgn, y, v, net, w_line, s_line):
        """ Edge-SMA footprint at a board edge. Returns (x where the line starts, half-width of the F.Cu clearance). """
        wp, lp, g, pull = v['wp'], v['lp'], v['g'], v['pull']
        x0 = x_edge + sgn * pull                       # copper starts here
        air = v.get('air', 0.0)
        pad_len = air if air > 0 else lp               # with the air notch the pad is only the overlap at the pin tip
        cx = x0 + sgn * (lp - pad_len / 2)
        n = self.net(net)
        if v.get('fp') == 'wurth':
            self._wurth(ref, x0, sgn, y, v, net, pad_len)
        elif v.get('fp') == 'bwsma':
            self._bwsma(ref, x_edge, sgn, y, v, net, (cx - x_edge) * sgn, pad_len)
        else:
            self.items.append('''  (footprint "z2labs:SMA_Edge_Launch" (layer "F.Cu") (uuid "%s") (at %g %g %d)
    (property "Reference" "%s" (at 0 -3 0) (layer "F.Fab") (uuid "%s") (effects (font (size 0.6 0.6) (thickness 0.1))))
    (property "Value" "SMA edge pin %.2f" (at 0 3 0) (layer "F.Fab") (uuid "%s") (effects (font (size 0.6 0.6) (thickness 0.1))))
    (attr smd)
    (pad "1" smd rect (at 0 0) (size %g %g) (layers "F.Cu" "F.Mask") (net %d "%s") (uuid "%s"))
    (model "${KIPRJMOD}/%s" (offset (xyz 0 0 0)) (scale (xyz 1 1 1)) (rotate (xyz 0 0 0)))
  )''' % (U(), cx, y, 0 if sgn > 0 else 180, ref, U(), v['pin'], U(), pad_len, wp, n, net, U(), model_name(v)))
            self.models.add((v['pin'], lp))
        # coplanar clearance around the pad and the taper (the F.Cu GND zone is everything else)
        xa = x0 + sgn * lp
        xb = xa + sgn * v['taper']
        hw = wp / 2 + g
        hn = (v.get('air_w') or (wp + 2 * g)) / 2 if air > 0 else hw   # notch half-width (the GND wall is 0.25 mm further out)
        self.hole(['F.Cu'], x_edge, y - hw, xa, y + hw)
        if v['taper'] > 0:
            # the taper overlaps the pad by 0.2 mm: polygons that only touch along an edge are not
            # reliably connected on the mesh (no mesh line sits exactly on a pad edge, thirds rule)
            xo = xa - sgn * min(0.2, pad_len / 2)
            self.poly([(xo, y - wp / 2), (xa, y - wp / 2), (xb, y - w_line / 2), (xb, y + w_line / 2), (xa, y + wp / 2), (xo, y + wp / 2)], net)
            self.hole(['F.Cu'], xa, y - hw, xb, y + hw)
        else:
            xb = xa - sgn * min(0.2, pad_len / 2)          # the line starts inside the pad
            # the pad ends where the narrower line gap starts: keep the coplanar gap g around the pad end too,
            # otherwise the pad corners touch the F.Cu GND fill of the line (short, KiCad DRC clearance)
            self.hole(['F.Cu'], xa, y - hw, xa + sgn * g, y + hw)
        # first line fence via must clear the pad clearance window (pad corner + via radius + 0.2 mm)
        self.fence_from = xa + sgn * (g if v['taper'] <= 0 else v['taper']) + sgn * 0.5
        # GND via pairs beside the pad (in the coplanar ground)
        vy = max(hw, hn) + (v.get('via_dy') or (0.7 if air > 0 else 0.45))   # clear of the routed notch (copper 0.25 mm from the edge)
        for k in range(v['vias']):
            vx = x0 + sgn * (0.6 + k * 1.0)
            self.via(vx, y - vy)
            self.via(vx, y + vy)
        # reference-plane cut-out under the pad
        if v['cut']:
            cl = lp + v['cut_dl']
            cw = wp + 2 * g + v['cut_dw']
            self.hole(list(v['cut']), x_edge, y - cw / 2, x_edge + sgn * cl, y + cw / 2)
        if pull > 0:                                   # all copper back from the edge
            self.hole(['F.Cu', 'In1.Cu', 'In2.Cu', 'B.Cu'], x_edge, y - PITCH, x0, y + PITCH)
        if air > 0:                                    # routed notch under the pin, copper 0.25 mm from the milled edge
            depth = lp - air
            self.notches.append((x_edge, sgn, y - hn, y + hn, depth))
            self.hole(list(v.get('air_layers') or ('F.Cu', 'In1.Cu', 'In2.Cu', 'B.Cu')), x_edge, y - hn - 0.25, x_edge + sgn * (depth + 0.25), y + hn + 0.25)
        return xb, hw

    def _wurth(self, ref, x0, sgn, y, v, net, pad_len):
        """ The Wuerth 60312202114509 footprint at the edge x0, barrel pointing outward (local +y of the library
            footprint is the barrel side, the board edge is at local y = 2.8). Pad angles in the file are absolute. """
        W = WURTH
        rot = 270 if sgn > 0 else 90                   # local +y -> global -sgn * x
        pa = (rot + 90) % 360
        cx = x0 + sgn * W['lp'] / 2
        n = self.net(net)
        ys = -W['lp'] / 2 + pad_len / 2                # the signal pad, shortened at the pin tip when air > 0
        legs = ''.join('    (pad "2" smd rect (at %g 0 %d) (size %g %g) (layers "%s" "%s") (net 1 "GND") (uuid "%s"))\n'
                       % (lx, pa, W['lp'], W['leg_w'], l, m, U()) for lx in (-W['leg_y'], W['leg_y']) for l, m in (('F.Cu', 'F.Mask'), ('B.Cu', 'B.Mask')))
        fab = ''.join('    (fp_line (start %g %g) (end %g %g) (stroke (width 0.1) (type default)) (layer "F.Fab") (uuid "%s"))\n' % (a, b, c, d, U())
                      for a, b, c, d in ((-4.76, -2.58, 4.76, -2.58), (-4.76, 2.8, -4.76, -2.58), (-4.76, 2.8, 4.76, 2.8), (4.76, -2.58, 4.76, 2.8),
                                         (-4.76, 4.07, -4.76, 2.8), (4.76, 2.8, 4.76, 4.07), (-4.76, 4.07, -3.175, 4.07), (3.175, 4.07, 4.76, 4.07),
                                         (-3.175, 4.07, -3.175, 11.69), (3.175, 4.07, 3.175, 11.69), (3.175, 11.69, -3.175, 11.69)))
        self.items.append('''  (footprint "Connector_SMD_Wurth:J_Wurth_WR-SMA_60312202114509" (layer "F.Cu") (uuid "%s") (at %g %g %d)
    (descr "PCB End Launch, https://www.we-online.com/en/components/products/datasheet/60312202114509")
    (property "Reference" "%s" (at 0 -4.445 0) (layer "F.Fab") (uuid "%s") (effects (font (size 0.6 0.6) (thickness 0.1))))
    (property "Value" "WR-SMA 60312202114509" (at 0 12.7 0) (layer "F.Fab") (uuid "%s") (effects (font (size 0.6 0.6) (thickness 0.1))))
    (attr smd)
%s    (pad "1" smd rect (at 0 %g %d) (size %g %g) (layers "F.Cu" "F.Mask") (net %d "%s") (uuid "%s"))
%s    (model "${KIPRJMOD}/%s" (offset (xyz 0 0 0)) (scale (xyz 1 1 1)) (rotate (xyz 0 0 0)))
  )''' % (U(), cx, y, rot, ref, U(), U(), fab, ys, pa, pad_len, W['wp'], n, net, U(), legs, W['step']))
        self.step = True

    def _bwsma(self, ref, x_edge, sgn, y, v, net, pad_x, pad_len):
        """ The BAT Wireless BWSMA-KE-P001 footprint, origin at the board edge on the pin axis, local +x into the
            board (rotation 0 on the left edge, 180 on the right). Signal pad wp x pad_len centred at local pad_x;
            ground leg pads 1.3 x 4.6 at y = +-2.55 on F.Cu (legs in the pin plane) and B.Cu (legs under the board). """
        W = BWSMA
        rot = 0 if sgn > 0 else 180
        n = self.net(net)
        ll, lw = W['leg_pad']
        legs = ''.join('    (pad "2" smd rect (at %g %g %d) (size %g %g) (layers "%s" "%s" "%s") (net 1 "GND") (uuid "%s"))\n'
                       % (ll / 2, sy * W['leg_y'], rot, ll, lw, l, p, m, U())
                       for sy in (-1, 1) for l, p, m in (('F.Cu', 'F.Paste', 'F.Mask'), ('B.Cu', 'B.Paste', 'B.Mask')))
        fab = ''.join('    (fp_line (start %g %g) (end %g %g) (stroke (width 0.1) (type default)) (layer "F.Fab") (uuid "%s"))\n' % (a, b, c, d, U())
                      for a, b, c, d in ((-1.5, -3.25, 0, -3.25), (0, -3.25, 0, 3.25), (0, 3.25, -1.5, 3.25), (-1.5, 3.25, -1.5, -3.25),
                                         (-1.5, -3.125, -9.5, -3.125), (-9.5, -3.125, -9.5, 3.125), (-9.5, 3.125, -1.5, 3.125),
                                         (0, -0.475, W['pin_len'], -0.475), (W['pin_len'], -0.475, W['pin_len'], 0.475), (W['pin_len'], 0.475, 0, 0.475)))
        crt = '    (fp_rect (start -9.75 -3.5) (end %g 3.5) (stroke (width 0.05) (type default)) (fill no) (layer "F.CrtYd") (uuid "%s"))\n' % (ll + 0.25, U())
        self.items.append('''  (footprint "z2labs:SMA_BAT_Wireless_BWSMA-KE-P001_EdgeMount" (layer "F.Cu") (uuid "%s") (at %g %g %d)
    (descr "SMA jack edge mount, offset pin, 1.6 mm board, BAT Wireless BWSMA-KE-P001 (LCSC C496550)")
    (property "Reference" "%s" (at 2 -4.2 %d) (layer "F.Fab") (uuid "%s") (effects (font (size 0.6 0.6) (thickness 0.1))))
    (property "Value" "BWSMA-KE-P001" (at 2 4.2 %d) (layer "F.Fab") (uuid "%s") (effects (font (size 0.6 0.6) (thickness 0.1))))
    (property "LCSC" "C496550" (at 0 0 %d) (layer "F.Fab") (hide yes) (uuid "%s") (effects (font (size 0.6 0.6) (thickness 0.1))))
    (attr smd)
%s%s    (pad "1" smd rect (at %g 0 %d) (size %g %g) (layers "F.Cu" "F.Paste" "F.Mask") (net %d "%s") (uuid "%s"))
%s    (model "${KIPRJMOD}/%s" (offset (xyz 0 0 0)) (scale (xyz 1 1 1)) (rotate (xyz 0 0 0)))
  )''' % (U(), x_edge, y, rot, ref, rot, U(), rot, U(), rot, U(), fab, crt, pad_x, rot, pad_len, v['wp'], n, net, U(), legs, W['step']))
        self.bwsma = True

    # ------------------------------------------------------------ zones
    def _fills(self, layer, cx0, y0, y1, L, lim=None):
        """ Rectangle [cx0, cx0+L] x [y0, y1] minus the holes of the layer, as merged x slabs. """
        x_lo, x_hi = cx0, cx0 + L
        hs = [h for h in self.holes[layer] if h[2] > x_lo and h[0] < x_hi and h[3] > y0 and h[1] < y1]
        if lim:                                        # the zone exists only inside these x intervals
            inv = []
            edges = sorted(set([x_lo, x_hi] + [c for a, b in lim for c in (a, b)]))
            for a, b in zip(edges[:-1], edges[1:]):
                if not any(la <= a and b <= lb for la, lb in lim):
                    inv.append((a, y0, b, y1))
            hs += inv
        xs = sorted({x_lo, x_hi} | {min(max(c, x_lo), x_hi) for h in hs for c in (h[0], h[2])})
        out, prev = [], None
        for xa, xb in zip(xs[:-1], xs[1:]):
            if xb - xa < 1e-9:
                continue
            cut = sorted((max(h[1], y0), min(h[3], y1)) for h in hs if h[0] < xb - 1e-9 and h[2] > xa + 1e-9)
            iv, y = [], y0
            for a, b in cut:
                if a > y + 1e-9:
                    iv.append((y, a))
                y = max(y, b)
            if y < y1 - 1e-9:
                iv.append((y, y1))
            if prev and prev[2] == iv and abs(prev[1] - xa) < 1e-9:
                prev[1] = xb
            else:
                prev = [xa, xb, iv]
                out.append(prev)
        return [(xa, ya, xb, yb) for xa, xb, iv in out for ya, yb in iv]

    def text_pcb(self, outline, vcuts=()):
        edges = ''.join('  (gr_line (start %g %g) (end %g %g) (stroke (width 0.05) (type default)) (layer "Edge.Cuts") (uuid "%s"))\n'
                        % (a[0], a[1], b[0], b[1], U()) for poly in outline for a, b in zip(poly, poly[1:] + poly[:1]))
        for (xa, ya), (xb, yb) in vcuts:
            edges += '  (gr_line (start %g %g) (end %g %g) (stroke (width 0.1) (type default)) (layer "User.1") (uuid "%s"))\n' % (xa, ya, xb, yb, U())
        zones = ''
        for ci, (cx0, y0, L) in enumerate(self.coupons):
            ya, yb = y0 + VCUT_CLEAR, y0 + PITCH - VCUT_CLEAR
            for layer in ('F.Cu', 'In1.Cu', 'In2.Cu', 'B.Cu'):
                rects = self._fills(layer, cx0, ya, yb, L, self.zone_lim.get(ci, {}).get(layer))
                if not rects:
                    continue
                fills = ''.join('\n    (filled_polygon (layer "%s") (pts (xy %g %g) (xy %g %g) (xy %g %g) (xy %g %g)))'
                                % (layer, a, b, c, b, c, d, a, d) for a, b, c, d in rects)
                zones += '''  (zone (net 1) (net_name "GND") (layer "%s") (uuid "%s") (hatch edge 0.5)
    (connect_pads yes (clearance 0.2)) (min_thickness 0.1) (filled_areas_thickness no) (fill yes (thermal_gap 0.3) (thermal_bridge_width 0.3))
    (polygon (pts (xy %g %g) (xy %g %g) (xy %g %g) (xy %g %g)))%s)
''' % (layer, U(), cx0, ya, cx0 + L, ya, cx0 + L, yb, cx0, yb, fills)
        # the holes as rule areas too, so a refill in KiCad keeps them
        for layer, hs in self.holes.items():
            for a, b, c, d in hs:
                zones += '''  (zone (net 0) (net_name "") (layer "%s") (uuid "%s") (name "nocu") (hatch edge 0.5)
    (connect_pads (clearance 0)) (min_thickness 0.1) (filled_areas_thickness no)
    (keepout (tracks allowed) (vias allowed) (pads allowed) (copperpour not_allowed) (footprints allowed))
    (polygon (pts (xy %g %g) (xy %g %g) (xy %g %g) (xy %g %g))))
''' % (layer, U(), a, b, c, b, c, d, a, d)
        nets = ' '.join('(net %d "%s")' % (i + 1, n) for i, n in enumerate(self.nets))
        return '''(kicad_pcb (version 20241229) (generator "pcbnew") (generator_version "9.0")
  (general (thickness 1.6) (legacy_teardrops no))
  (paper "A4")
  (layers (0 "F.Cu" signal) (4 "In1.Cu" signal) (6 "In2.Cu" signal) (2 "B.Cu" signal)
    (5 "F.SilkS" user) (7 "B.SilkS" user) (1 "F.Mask" user) (3 "B.Mask" user) (25 "Edge.Cuts" user)
    (31 "F.CrtYd" user) (35 "F.Fab" user) (39 "User.1" user "V-Cut"))
  (setup
%s
    (pad_to_mask_clearance 0))
  (net 0 "") %s
%s
%s%s)
''' % (STACKUP, nets, '\n'.join(self.items), edges, zones)


# ---------------------------------------------------------------- coupon builders
def _fence(b, y, xs, xe, pitch=1.0, x0=None):
    x = xs if x0 is None else x0
    while x <= xe + 1e-6:
        b.via(x, y)
        x += pitch


def build_coupon(b, ci, cid, label, v, kind, p, x0, y0):
    """ One coupon at (x0, y0); the SMA edges are x0 and x0 + L. """
    L = p.get('L', LC)
    s = p.get('s', S0)
    w = GCPW[s] if kind != 'msl' else W50
    y = y0 + PITCH / 2
    b.coupons.append((x0, y0, L))
    net = 'RF%d' % cid
    xs, hw = b.launch('J%dA' % cid, x0, 1, y, v, net, w, s)
    fa = max(xs + 0.5, b.fence_from)
    xe, _ = b.launch('J%dB' % cid, x0 + L, -1, y, v, net, w, s)
    fb = min(xe - 0.5, b.fence_from)
    xf = max(4.0, v['lp'] + 0.8)                  # outer via fences start clear of the launch footprint
    b.text('%d %s' % (cid, label), x0 + 6.0, y0 + 1.5, size=0.8)
    if kind in ('gcpw', 'msl'):
        b.seg(xs, y, xe, y, w, net)
    if kind == 'gcpw':
        b.hole(['F.Cu'], xs, y - w / 2 - s, xe, y + w / 2 + s)
        for sg in (-1, 1):
            _fence(b, y + sg * (w / 2 + s + 0.55), fa, fb, 1.0)                 # < lambda/10 at 10 GHz
            _fence(b, y + sg * 3.5, x0 + xf, x0 + L - xf, 2.0)
    elif kind == 'msl':
        # coplanar ground only around the launches (the connector legs), open microstrip elsewhere
        lim = [(x0, xs + 1.5), (xe - 1.5, x0 + L)]
        b.zone_lim[ci] = {'F.Cu': lim}
        b.hole(['F.Cu'], xs, y - hw, xe, y + hw)
        for sg in (-1, 1):
            _fence(b, y + sg * 3.5, x0 + xf, x0 + L - xf, 2.0)
    elif kind == 'via':
        x1, x2 = x0 + L / 3, x0 + 2 * L / 3
        b.seg(xs, y, x1, y, w, net)
        b.seg(x1, y, x2, y, w, net, layer='B.Cu')
        b.seg(x2, y, xe, y, w, net)
        b.hole(['F.Cu'], xs, y - w / 2 - s, xe, y + w / 2 + s)
        b.hole(['B.Cu'], x1 - 0.8, y - w / 2 - s, x2 + 0.8, y + w / 2 + s)   # B.Cu GCPW gap
        for xv in (x1, x2):
            b.via(xv, y, net, size=0.6, drill=0.3)
            b.hole(['In1.Cu', 'In2.Cu'], xv - 0.6, y - 0.6, xv + 0.6, y + 0.6)        # antipad
            b.hole(['F.Cu', 'B.Cu'], xv - 0.6, y - 0.6, xv + 0.6, y + 0.6)
            for dx in (-0.95, 0.95):
                for dy in (-0.95, 0.95):
                    b.via(xv + dx, y + dy)                                           # GND return vias
        for sg in (-1, 1):
            _fence(b, y + sg * (w / 2 + s + 0.55), xs + 0.5, xe - 0.5, 1.0)
            _fence(b, y + sg * 3.5, x0 + xf, x0 + L - xf, 2.0)
    elif kind == 'lpf':
        total = sum(l for _, l in LPF_SECT)
        xf = x0 + (L - total) / 2
        b.seg(xs, y, xf, y, w, net)
        x = xf
        for k, l in LPF_SECT:
            if k == 'C':
                b.rect_cu(x, y - LPF_LO / 2, x + l, y + LPF_LO / 2, net)
            else:
                b.seg(x, y, x + l, y, LPF_HI, net)
            x += l
        b.seg(x, y, xe, y, w, net)
        b.hole(['F.Cu'], xs, y - w / 2 - s, xe, y + w / 2 + s)
        b.hole(['F.Cu'], xf - 1.0, y - LPF_LO / 2 - 2.0, x + 1.0, y + LPF_LO / 2 + 2.0)   # filter sits in open microstrip
        for sg in (-1, 1):
            _fence(b, y + sg * (w / 2 + s + 0.55), xs + 0.5, xf - 1.5, 1.0)
            _fence(b, y + sg * (w / 2 + s + 0.55), x + 1.5, xe - 0.5, 1.0)
            _fence(b, y + sg * 4.0, x0 + xf, x0 + L - xf, 2.0)
    return L


def build_panel(coupons=COUPONS):
    """ One column of coupons, every connector on an outer edge; the right-hand outline follows the
        coupon lengths (short coupons at the end). Returns (board, outline polygons, V-cut lines). """
    b = Board()
    for i, (cid, label, v, kind, p) in enumerate(coupons):
        build_coupon(b, i, cid, label, v, kind, p, 0.0, i * PITCH)
    Hh = len(coupons) * PITCH
    outline = [column_outline(b, [c[4].get('L', LC) for c in coupons])]
    vcuts = [((0, r * PITCH), (LC, r * PITCH)) for r in range(1, len(coupons))]
    b.text('z2labs SMA launch DoE  JLC04161H-7628  V-cut rows', 2.0, Hh - 0.6, size=0.7)
    return b, outline, vcuts


def column_outline(b, lengths):
    """ Outline of one column of coupons (x from 0, rows of PITCH), right side following the lengths,
        with the launch notches cut into both edges; clockwise in KiCad coordinates. """
    Hh = len(lengths) * PITCH
    right = sorted([n for n in b.notches if n[1] < 0], key=lambda n: n[2])
    left = sorted([n for n in b.notches if n[1] > 0], key=lambda n: -n[2])
    pts = [(0.0, 0.0)]
    for i, L in enumerate(lengths):                        # down the right edge
        y0, y1 = i * PITCH, (i + 1) * PITCH
        pts.append((L, y0))
        for xe, sgn, ylo, yhi, d in right:
            if y0 - 1e-9 <= ylo and yhi <= y1 + 1e-9 and abs(xe - L) < 1e-9:
                pts += [(L, ylo), (L - d, ylo), (L - d, yhi), (L, yhi)]
        pts.append((L, y1))
    pts.append((0.0, Hh))
    for xe, sgn, ylo, yhi, d in left:                     # up the left edge
        pts += [(0.0, yhi), (d, yhi), (d, ylo), (0.0, ylo)]
    return [p for k, p in enumerate(pts) if k == 0 or (abs(p[0] - pts[k - 1][0]) > 1e-9 or abs(p[1] - pts[k - 1][1]) > 1e-9)]


def build_single(c):
    b = Board()
    cid, label, v, kind, p = c
    L = build_coupon(b, 0, cid, label, v, kind, p, 0.0, 0.0)
    return b, [column_outline(b, [L])]


def write_all(out, name='sma_doe'):
    os.makedirs(out, exist_ok=True)
    files = {}
    half = (len(COUPONS) + 1) // 2
    for tag, cs in (('A', COUPONS[:half]), ('B', COUPONS[half:])):
        b, ol, vc = build_panel(cs)
        files['panel_' + tag] = os.path.join(out, '%s_%s.kicad_pcb' % (name, tag))
        open(files['panel_' + tag], 'w').write(b.text_pcb(ol, vc))
        write_models(b, out)
    for c in COUPONS:
        b, ol = build_single(c)
        p = os.path.join(out, '%s_%02d.kicad_pcb' % (name, c[0]))
        open(p, 'w').write(b.text_pcb(ol))
        files[c[0]] = p
    rows = ['| # | coupon | pin | wp | lp | g | cut | taper | vias | pull | air ov | notch w | line |', '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for cid, label, v, kind, p in COUPONS:
        cut = '+'.join(l[:-3] for l in v['cut']) + (' %+g/%+g' % (v['cut_dl'], v['cut_dw']) if v['cut'] else '') or 'none'
        line = kind + (' s%g' % p['s'] if 's' in p else '') + (' %gmm' % p['L'] if 'L' in p else '')
        rows.append('| %d | %s | %g | %g | %g | %g | %s | %g | %d | %g | %s | %s | %s |' % (cid, label, v['pin'], v['wp'], v['lp'], v['g'], cut, v['taper'], v['vias'], v['pull'], v.get('air', 0) or '-',
                    ('%g' % (v.get('air_w') or v['wp'] + 2 * v['g'])) if v.get('air') else '-', line))
    open(os.path.join(out, name + '_table.md'), 'w').write('\n'.join(rows) + '\n')
    return files


# ---------------------------------------------------------------- simulation of one coupon
def sim_coupon(c, out, fmax=10e9, variant=None, threads=4, both=False, log=print):
    """ openEMS run of a standalone coupon with the connector model of launch_study at both ends.
        Symmetric coupons need only the port-1 excitation (S11, S21); both=True runs port 2 too. """
    import numpy as np
    import kicad_sim as ks
    import launch_study as ls
    ks._setup_python()
    from CSXCAD import ContinuousStructure
    from CSXCAD.KiCadImport import KiCadBoard
    from openEMS import openEMS
    cid, label, v, kind, p = c
    L = p.get('L', LC)
    os.makedirs(out, exist_ok=True)
    pcb = os.path.join(out, 'c%02d.kicad_pcb' % cid)
    b, ol = build_single(c)
    open(pcb, 'w').write(b.text_pcb(ol))
    board = KiCadBoard(pcb)
    f = np.linspace(0.2e9, fmax, int(os.environ.get('Z2_NFREQ', '491')))     # Z2_NFREQ: finer grid for resonators
    S = np.zeros((len(f), 2, 2), complex)
    meta = {}
    for e in ((1, 2) if both else (1,)):
        t_max = max(3e-9, 15 * math.hypot(L, PITCH) * 1e-3 * math.sqrt(ER) / ls.C0)
        FDTD = openEMS(EndCriteria=1e-4, NrTS=400000, MaxTime=t_max, **FDTD_KW)
        FDTD.SetGaussExcite(0.5 * fmax + 0.15e9, 0.5 * fmax - 0.15e9)      # no DC: the floating strips would keep a static charge
        bc = {'mur': ['MUR'] * 6, 'pml': ['PML_8'] * 6, 'pec': ['PEC'] * 6}.get(os.environ.get('Z2_BC', '').lower(), ls.BOUNDARY)
        FDTD.SetBoundaryCond(bc)
        CSX = ContinuousStructure()
        FDTD.SetCSX(CSX)
        board.ToCSX(CSX, frequency=float(os.environ.get('Z2_LOSS_F', 0.5 * fmax)), copper=os.environ.get('Z2_COPPER', 'sheet'))     # Z2_COPPER=pec for A/B
        fake = [dict(net='RF%d' % cid, ref_net='GND', start=[0, 0, 0], stop=[0, 0, 0], normal_axis=0)]
        X, Y, Z, res = ks.mesh_lines(board, fake, fmax, margin=MARGIN, coarse=float(os.environ.get('Z2_COARSE', '1')))
        X, Y, Z = list(X), list(Y), list(Z)
        pdefs = []
        lv = dict(pin=v['pin'], lp=v['lp'] + v['pull'])
        if v.get('fp') == 'wurth':                           # real flange and ground legs of the Wuerth part
            lv.update(flange=WURTH['flange'], legs=(WURTH['leg_y'], 1.6, WURTH['lp'], 0.5))
        elif v.get('fp') == 'bwsma':                         # real flange, 0.9 mm square legs, 4 mm pin
            lv.update(flange=BWSMA['flange'], legs=(BWSMA['leg_y'], BWSMA['leg'], BWSMA['leg_len'], BWSMA['leg']),
                      pin_len=BWSMA['pin_len'])
        ext_y, ext_z = [], []
        for x_edge, sgn in ((0.0, 1), (L, -1)):
            pd, x2, y2, z2 = ls.connector(CSX, board, lv, x_edge, sgn, yc=-PITCH / 2)
            pdefs.append(pd)
            X += x2; Y += y2; Z += z2
            ext_y += [min(y2), max(y2)]; ext_z += [min(z2), max(z2)]
        # Air margin measured from the connector as well: the flange (5.3 mm above the board) otherwise ends one cell
        # below the z boundary (MUR) or inside the PML, and the run grows exponentially after ~40k steps (2026-10-21,
        # A/B runs: all-PML stable, MUR/PEC walls unstable; PML_8 on all sides is the default now).
        for xe, sgn, ylo, yhi, d in b.notches:              # mesh lines on the notch walls
            X.append(xe + sgn * d)
            Y += [-ylo, -yhi]
        g = CSX.GetGrid()
        g.SetDeltaUnit(1e-3)
        for ax, lines in zip('xyz', (X, Y, Z)):
            sm = np.unique(np.round(ks.smooth(lines, res, merge=0.05), 5))
            if ax == 'y':
                sm = ks.pml_pad(sm, min(ext_y), max(ext_y), n=14)
            elif ax == 'z':
                sm = ks.pml_pad(sm, min(ext_z), max(ext_z), n=14)
            g.SetLines(ax, sm)
        ports = ls.make_ports(FDTD, pdefs, excite=e)
        nl = [len(g.GetLines(a)) for a in 'xyz']
        cdir = os.path.join(out, 'c%02d' % cid, 'exc%d' % e)
        os.makedirs(cdir, exist_ok=True)
        dumps = []
        if os.environ.get('Z2_DUMPS', '0') != '0':             # current-density maps on every copper plane (fields.py)
            import fields as fmod
            dump_f = [float(x) * 1e9 for x in os.environ.get('Z2_DUMP_F', '1,3,6').split(',') if float(x) * 1e9 <= fmax]
            dumps = fmod.add_dumps(CSX, board, dump_f)
        xml = os.path.join(cdir, 'model.xml')
        FDTD.Write2XML(xml)
        log('c%02d %s exc%d: mesh %s (%.2f MCells)' % (cid, label, e, nl, np.prod(nl) / 1e6))
        run_info = dict(id=cid, label=label, kind=kind, launch={k: v2 for k, v2 in v.items() if not isinstance(v2, (list, tuple, dict))},
                        params=p, exc=e, mesh=nl, mcells=float(np.prod(nl) / 1e6), fmax_GHz=fmax / 1e9, pcb=pcb,
                        started=time.strftime('%Y-%m-%d %H:%M:%S'), status='running', variant=variant,
                        env={k: os.environ[k] for k in os.environ if k.startswith('Z2_')})
        json.dump(run_info, open(os.path.join(cdir, 'run.json'), 'w'), indent=1)
        r = ks.run_engine(xml, os.path.join(cdir, 'sim'), variant, threads)
        log('  rc=%s %.0f s, %s TS' % (r['rc'], r['wall_s'], r['timesteps']))
        run_info.update(status='done' if r['rc'] == 0 else 'failed', rc=r['rc'], wall_s=r['wall_s'], timesteps=r['timesteps'],
                        finished=time.strftime('%Y-%m-%d %H:%M:%S'))
        if dumps:
            try:
                import fields as fmod
                run_info['fields'] = fmod.export(os.path.join(cdir, 'sim'), os.path.join(cdir, 'fields'), board, dumps)
            except Exception as ex:
                log('  fields export failed: %s' % ex)
        json.dump(run_info, open(os.path.join(cdir, 'run.json'), 'w'), indent=1)
        for pp in ports:
            pp.CalcPort(os.path.join(cdir, 'sim'), f, ref_impedance=50)
        for i, pp in enumerate(ports):
            S[:, i, e - 1] = pp.uf_ref / ports[e - 1].uf_inc
        meta = dict(mesh=nl, mcells=float(np.prod(nl) / 1e6), wall_s=r['wall_s'])
    if not both:                                       # symmetric two-port
        S[:, 1, 1], S[:, 0, 1] = S[:, 0, 0], S[:, 1, 0]
    ks.write_touchstone(os.path.join(out, 'c%02d.s2p' % cid), f, S, comment='make_panel coupon %d %s %s' % (cid, label, json.dumps(v)))
    s11 = 20 * np.log10(np.abs(S[:, 0, 0]) + 1e-15)
    s21 = 20 * np.log10(np.abs(S[:, 1, 0]) + 1e-15)
    res = dict(id=cid, label=label, launch=v, kind=kind, params=p, f_GHz=(f / 1e9).round(4).tolist(), S11_dB=s11.round(2).tolist(),
               S21_dB=s21.round(3).tolist(), s11_max_3G=float(s11[f <= 3e9].max()), s11_max_6G=float(s11[f <= 6e9].max()),
               s11_max_10G=float(s11.max()), s21_min_6G=float(s21[f <= 6e9].min()), **meta)
    json.dump(res, open(os.path.join(out, 'c%02d.json' % cid), 'w'))
    try:                                                   # run manifest for the monitor (z2mon.py)
        rj = os.path.join(out, 'c%02d' % cid, 'exc1', 'run.json')
        ri = json.load(open(rj))
        ri.update(s11_max_3G=res['s11_max_3G'], s11_max_6G=res['s11_max_6G'], s11_max_10G=res['s11_max_10G'], s21_min_6G=res['s21_min_6G'],
                  s2p=os.path.join(out, 'c%02d.s2p' % cid), result=os.path.join(out, 'c%02d.json' % cid))
        json.dump(ri, open(rj, 'w'), indent=1)
    except Exception:
        pass
    log('  c%02d %s: S11 max %.1f / %.1f / %.1f dB (<3/6/10 GHz), S21 min %.2f dB (<6 GHz)'
        % (cid, label, res['s11_max_3G'], res['s11_max_6G'], res['s11_max_10G'], res['s21_min_6G']))
    return res


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'panel_out')
    files = write_all(out)
    print('\n'.join('%s %s' % kv for kv in files.items()))
    if '--sim' in sys.argv:
        want = [int(x) for x in sys.argv[sys.argv.index('--sim') + 1].split(',')]
        for c in COUPONS:
            if c[0] in want:
                sim_coupon(c, os.path.join(out, 'sim'), variant=os.environ.get('Z2_VARIANT') or None,
                           threads=int(os.environ.get('Z2_THREADS', '4')), log=lambda m: print(m, flush=True))
