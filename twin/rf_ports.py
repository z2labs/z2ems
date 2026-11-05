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
rf_ports: find the RF ports of a KiCad board and propose lumped ports for them.

Shared by the command line (kicad_sim) and the viewer (through Pyodide): the same evidence gives
the same suggestions everywhere. Every candidate is scored from independent clues, each with a
one-line reason, and nothing is forced: a candidate with score >= CONFIDENT is proposed checked,
the rest unchecked.

Clues (points):
  footprint library / name says coaxial connector (SMA, U.FL/IPEX/MHF, MMCX, SMP, BNC, N, coax)   3
  footprint value says the same                                                                  1
  one signal pad plus ground pads (connector-like pad set)                                       1
  the signal pad is within 2.5 mm of the board edge                                              1
  net name looks RF (RF, ANT, COAX, LNA, PA, TX/RX with RF, GPS, LTE, WIFI, BLE, ...)            2
  the net is in a net class whose name or rules say RF / 50 Ohm / impedance (from .kicad_pro)    2

    python rf_ports.py board.kicad_pcb [board.kicad_pro]
"""

import fnmatch
import json
import os
import re
import sys

import numpy as np

CONFIDENT = 4

FP_RE = re.compile(r'(^|[^a-z])(sma|u\.?fl|ufl|ipex|i-pex|mhf\d?|mmcx|mcx|smp|smpm|bnc|n[-_ ]?type|tnc|coax|rf[-_ ]?conn\w*|rf[-_ ]?(jack|plug|antenna))(?![a-z])', re.I)
LIB_RE = re.compile(r'coax|rf', re.I)
NET_RE = re.compile(r'(^|[^a-z])(rf|ant|antenna|coax|lna|pa_?out|pa_?in|tx_?rf|rx_?rf|rfin|rfout|gps|gnss|lte|wifi|wlan|ble|lora|ism|nfc)', re.I)
GROUND = re.compile(r'^(gnd|ground|agnd|dgnd|vss|0v)', re.I)
CLASS_RE = re.compile(r'rf|50\s*(r|ohm|Ω)|coax|imped', re.I)
REF_RE = re.compile(r'^(J|X|P|CON|CN|ANT|RF|SMA|U\.?FL|XS|K)\d', re.I)


def _netclass_of(project, net):
    """ Net class name of `net` from a .kicad_pro dict (explicit assignment, then patterns), else None. """
    ns = (project or {}).get('net_settings') or {}
    a = ns.get('netclass_assignments') or {}       # null in projects without assignments
    if net in a:
        return a[net]
    for p in ns.get('netclass_patterns') or []:
        if fnmatch.fnmatch(net, p.get('pattern', '')):
            nc = p.get('netclass')
            return nc[0] if isinstance(nc, list) and nc else nc       # KiCad 10: a list of classes
    for c in ns.get('classes') or []:
        for pat in c.get('nets', []) or []:
            if fnmatch.fnmatch(net, pat):
                return c.get('name')
    return None


def _rf_class(project, name):
    if isinstance(name, list):                    # KiCad 10 netclass assignments are lists
        return any(_rf_class(project, n) for n in name)
    if not name or not project or not isinstance(name, str):
        return False
    if CLASS_RE.search(name):
        return True
    for c in project.get('net_settings', {}).get('classes', []):
        if c.get('name') == name:
            return any(CLASS_RE.search(str(c.get(k, ''))) for k in ('description', 'pcb_color', 'schematic_color')) or \
                any(CLASS_RE.search(str(v)) for v in c.values() if isinstance(v, str))
    return False


def load_project(board_path):
    """ The .kicad_pro next to the board, as a dict (None when absent or unreadable). """
    p = os.path.splitext(board_path)[0] + '.kicad_pro'
    try:
        return json.load(open(p, encoding='utf-8'))
    except Exception:
        return None


def _edge_distance(board, pts):
    (bx0, bx1), (by0, by1) = board.BoundingBox()[:2]
    cx, cy = pts.mean(0)
    return min(cx - bx0, bx1 - cx, cy - by0, by1 - cy)


def candidates(board, project=None):
    """ RF port candidates: list of dicts sorted by score (ref, pad, net, score, reasons, confident). """
    outer = (board.copper_layers[0], board.copper_layers[-1])
    by_ref = {}
    for s in board.shapes:
        if s.kind == 'pad' and s.layer in outer:
            by_ref.setdefault(s.ref, []).append(s)
    out = []
    for ref, pads in by_ref.items():
        F = board.footprints.get(ref)
        sig = [p for p in pads if p.net and not GROUND.match(p.net)]
        gnd = [p for p in pads if p.net and GROUND.match(p.net)]
        nets = sorted({p.net for p in sig})
        if not sig or len(nets) > 2:
            continue                                    # no signal, or clearly not a connector
        score, why = 0, []
        conn_ref = bool(REF_RE.match(ref))                  # J, X, P, CON, CN, ANT, RF ...: a connector by reference
        passive = bool(re.match(r'(R|C|L|D|FB|Q|U|Y|LED|TP)\d', ref, re.I))
        lib = (F.lib_id or '') if F else ''
        libname, fpname = (lib.split(':') + [''])[:2]
        if F and not passive and (FP_RE.search(fpname) or LIB_RE.search(libname)):
            score += 3; why.append('coaxial connector footprint %s' % F.lib_id)
        if F and not passive and FP_RE.search(F.value or ''):
            score += 2 if re.search(r'conn|jack|plug|antenna', F.value or '', re.I) else 1
            why.append('value "%s"' % F.value)
        if len(nets) == 1 and (len(gnd) >= 2 or (gnd and conn_ref)):
            score += 1; why.append('one signal pad and %d ground pad(s)' % len(gnd))
        # the signal pad nearest the board edge is the port pad
        pad = min(sig, key=lambda p: _edge_distance(board, p.pts))
        de = _edge_distance(board, pad.pts)
        if de <= 2.5 and (conn_ref or len(gnd) >= 2):
            score += 1; why.append('signal pad %.1f mm from the board edge' % de)
        if NET_RE.search(pad.net):
            score += 2; why.append('net name "%s"' % pad.net)
        nc = _netclass_of(project, pad.net)
        if _rf_class(project, nc):
            score += 2; why.append('net class "%s"' % nc)
        if score < 2:
            continue
        out.append(dict(ref=ref, pad=pad, net=pad.net, layer=pad.layer, score=score, reasons=why, confident=score >= CONFIDENT,
                        lib_id=F.lib_id if F else '', value=F.value if F else '', edge_mm=round(float(de), 2)))
    out.sort(key=lambda c: (-c['score'], c['ref']))
    return out


def snap(board, x, y, layer=None, exclude_ground=True):
    """ Magnetic snap for a hand-placed port: the nearest outer-layer pad to (x, y) (board mm, y up) and
        its port definition from kicad_sim.make_port. Returns (pad, port dict) or (None, None). """
    outer = (board.copper_layers[0], board.copper_layers[-1])
    best, bd = None, 1e9
    for s in board.shapes:
        if s.kind != 'pad' or s.layer not in outer or (layer and s.layer != layer):
            continue
        if exclude_ground and GROUND.match(s.net or ''):
            continue
        d = float(np.hypot(*(s.pts.mean(0) - (x, y))))
        if d < bd:
            best, bd = s, d
    if best is None:
        return None, None
    import kicad_sim as ks
    return best, ks.make_port(board, best, 1)


def report(board, project=None):
    lines = []
    for c in candidates(board, project):
        lines.append('%s %-6s %-14s score %d  %s  [%s]' % ('[x]' if c['confident'] else '[ ]', c['ref'], c['net'], c['score'],
                                                            c['lib_id'] or c['value'], '; '.join(c['reasons'])))
    return '\n'.join(lines) or '(no RF port candidates)'


if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import kicad_sim as ks
    ks._setup_python()
    from CSXCAD.KiCadImport import KiCadBoard
    b = KiCadBoard(sys.argv[1])
    pro = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else load_project(sys.argv[1])
    print(report(b, pro))
