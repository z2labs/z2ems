"""
coherence: phase-match check of channel paths on a KiCad board (fast, geometric).

A spec groups pad-to-pad copper hops into channels:

    {"f0": 2.44e9, "tol_deg": 1.0,
     "groups": [{"name": "noise_dist", "kind": "absolute",            # absolute: phase must match at f0
                 "channels": [{"name": "CH1", "hops": [{"net": "Noise1", "a": ["R107", "2"], "b": ["IC302", "2"]}]}]},
                {"name": "rf_chain", "kind": "delay", "band_hz": 83.5e6, ...}]}   # delay: calibrated at one
                                                                                     # frequency, slope over the band

Tolerance: absolute groups pass when the channel-to-channel phase spread at f0 is <= tol_deg; delay groups
when the delay spread keeps the phase change over band_hz <= tol_deg (spread_ps <= tol_deg / (360 * band)).
The spec is stored next to the board as <board>.coherence.json.
"""
import json
import os

import pathlen as pl
from kipcb_lite import Board

_cache = {}


def board(path):
    m = os.path.getmtime(path)
    if path not in _cache or _cache[path][0] != m:
        b = Board(path)
        _cache[path] = (m, b, pl.Stack(b.stackup()))
    return _cache[path][1], _cache[path][2]


def spec_path(pcb):
    return pcb[:-len('.kicad_pcb')] + '.coherence.json'


def load_spec(pcb):
    p = spec_path(pcb)
    if os.path.exists(p):
        return json.load(open(p, encoding='utf-8'))
    return dict(f0=2.44e9, tol_deg=1.0, groups=[])


def save_spec(pcb, spec):
    with open(spec_path(pcb), 'w', encoding='utf-8') as f:
        json.dump(spec, f, indent=1)


def compute(pcb, spec):
    b, st = board(pcb)
    f0 = float(spec.get('f0', 2.44e9))
    tol = float(spec.get('tol_deg', 1.0))
    out = []
    MAXCH = 20
    for g in spec.get('groups', []):
        if len(g.get('channels', [])) > MAXCH:
            raise ValueError('group %s: at most %d channels' % (g.get('name'), MAXCH))
        rows = []
        for ch in g.get('channels', []):
            r = dict(name=ch.get('name', ''), length_mm=0.0, delay_ps=0.0, vias=0, layers={}, hops=[], problems=[], draw=[])
            for h in ch.get('hops', []):
                a, bb = tuple(h['a']), tuple(h['b'])
                p = pl.path(b, h['net'], a, bb, st)
                if p and p['length_mm'] is None:
                    p = pl.path(b, h['net'], a, bb, st, gap_tol=0.5)
                if not p or p['length_mm'] is None:
                    r['problems'].append('%s: %s.%s -> %s.%s %s' % (h['net'], a[0], a[1], bb[0], bb[1],
                                                                  'pad not on this net' if not p else 'no copper path'))
                    continue
                for gp in p.get('gaps', []):
                    r['problems'].append('%s: copper gap %.3f mm on %s at %s' % (h['net'], gp['gap_mm'], gp['layer'], gp['a']))
                r['length_mm'] += p['length_mm']; r['delay_ps'] += p['delay_ps']; r['vias'] += p['vias']
                for l, v in p['layers'].items():
                    r['layers'][l] = r['layers'].get(l, 0.0) + v
                r['hops'].append(dict(net=h['net'], a=list(a), b=list(bb), length_mm=p['length_mm'], delay_ps=p['delay_ps'], vias=p['vias']))
                # copper path for the picture: [x1, y1, x2, y2, layer] in KiCad coordinates (y down), vias as points
                for ly, q0, q1 in p.get('pieces', []):
                    if q0 and q1:
                        r['draw'].append([round(q0[0], 3), round(q0[1], 3), round(q1[0], 3), round(q1[1], 3), ly])
            r['phase_deg'] = r['delay_ps'] * 1e-12 * f0 * 360
            rows.append(r)
        ds = [r['delay_ps'] for r in rows if not r['problems']]
        spread = (max(ds) - min(ds)) if len(ds) > 1 else 0.0
        kind = g.get('kind', 'absolute')
        if kind == 'delay':
            band = float(g.get('band_hz', 83.5e6))
            limit_ps = tol / (360 * band) * 1e12
            metric = spread * 1e-12 * band * 360
        else:
            limit_ps = tol / (360 * f0) * 1e12
            metric = spread * 1e-12 * f0 * 360
        ref = max(ds) if ds else 0.0
        for r in rows:
            r['to_longest_ps'] = ref - r['delay_ps']
            # extra length on the dominant layer of the path that would bring it to the longest one
            lay = max(r['layers'], key=r['layers'].get) if r['layers'] else st.cu[0]
            w = 0.2
            r['extend_mm'] = r['to_longest_ps'] / (pl.math.sqrt(st.eps_eff(lay, w)) / pl.C0 * 1e9)
            r['extend_layer'] = lay
        out.append(dict(name=g.get('name', ''), kind=kind, band_hz=g.get('band_hz'), channels=rows, spread_ps=spread,
                        spread_deg=metric, limit_ps=limit_ps, ok=(spread <= limit_ps) and all(not r['problems'] for r in rows)))
    xs = [q[0] for es in b.edges.values() for e in es for q in (e[1], e[2])]
    ys = [q[1] for es in b.edges.values() for e in es for q in (e[1], e[2])]
    bbox = [min(xs), min(ys), max(xs), max(ys)] if xs else [0, 0, 1, 1]
    return dict(f0=f0, tol_deg=tol, stackup_note=_stack_note(b), groups=out, bbox=bbox)


def pad_near(pcb, x, y, r_max=2.0):
    """ the pad under / nearest to (x, y) in KiCad coordinates (y down): (ref, pad, net) """
    b, _ = board(pcb)
    best = None
    for p in b.padgeo:
        d = pl.math.hypot(p['x'] - x, p['y'] - y) - p['r']
        if best is None or d < best[0]:
            best = (d, p)
    if best is None or best[0] > r_max:
        return None
    p = best[1]
    return dict(ref=p['ref'], pad=p['pad'], net=p['net'], x=p['x'], y=p['y'], layers=p['layers'])


def _stack_note(b):
    st = b.stackup()
    d = [l for l in st if 'dielectric' in l['name']]
    ers = sorted({l.get('epsilon_r') for l in d})
    ts = sorted({l.get('thickness') for l in d})
    if ers == ['4.5'] and len(ts) == 1:
        return 'stackup looks like KiCad defaults (eps_r 4.5, equal dielectrics): ps and degrees scale with sqrt(eps_r); mm differences do not'
    return ''
