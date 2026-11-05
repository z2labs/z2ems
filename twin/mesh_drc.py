#!/usr/bin/env python3
""" mesh_drc.py — connectivity check of the staircase mesh against the KiCad nets ("DRC on the mesh").

    The staircase FDTD grid turns a copper sheet into the grid edges whose both end nodes lie inside the copper.
    A gap narrower than the local cell can vanish (two nets shorted: via pad to the GND plane across the antipad)
    and a thin non-axis-aligned conductor can break (diagonal track, arc, ring). Both turned real coupons into
    a short or an open on 2026-10-21. This check labels the staircase copper islands per layer and reports

      SHORT  two different nets in one island
      OPEN   one net of one layer split into islands that its KiCad copper does not split (heuristic)

        python mesh_drc.py <model.xml> <board.kicad_pcb> [--png out.png]
    exit code 1 when a SHORT is found.
"""
import re, sys, json
import numpy as np
from matplotlib.path import Path
from scipy import ndimage

sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import runviz

LAYER_PROP = {'F.Cu': 'F_Cu', 'In1.Cu': 'In1_Cu', 'In2.Cu': 'In2_Cu', 'In3.Cu': 'In3_Cu', 'In4.Cu': 'In4_Cu', 'B.Cu': 'B_Cu'}


def kicad_points(pcb):
    """ representative copper points per (layer, net): segment midpoints and ends, via centres, filled-zone interior points """
    s = open(pcb, encoding='utf-8', errors='replace').read()
    nets = {a: b for a, b in re.findall(r'\(net (\d+) "([^"]*)"\)', s)}
    pts = []                          # (layer, net, x, y, kind)
    for m in re.finditer(r'\(segment\s*\(start ([-\d.]+) ([-\d.]+)\)\s*\(end ([-\d.]+) ([-\d.]+)\)\s*\(width ([\d.]+)\)\s*\(layer "?([^")]+)"?\)\s*\(net (\d+)\)', s):
        x0, y0, x1, y1 = map(float, m.groups()[:4]); L, n = m.group(6), m.group(7)
        for t in (0.15, 0.5, 0.85):
            pts.append((L, n, x0 + t * (x1 - x0), y0 + t * (y1 - y0), 'seg'))
    for m in re.finditer(r'\(via\s*\(at ([-\d.]+) ([-\d.]+)\)\s*\(size ([\d.]+)\)\s*\(drill ([\d.]+)\)\s*\(layers ([^)]*)\).*?\(net (\d+)\)', s, re.S):
        x, y = float(m.group(1)), float(m.group(2))
        for L in LAYER_PROP:                                        # through via: pad on every layer
            pts.append((L, m.group(6), x, y, 'via'))
    for m in re.finditer(r'\(zone\s*\(net (\d+)\).*?(?=\n  \(zone|\n  \(gr_|\n\)\s*$)', s, re.S):
        n = m.group(1)
        for fm in re.finditer(r'\(filled_polygon\s*\(layer "?([^")]+)"?\)\s*\(pts(.*?)\)\s*\)', m.group(0), re.S):
            P = np.array(re.findall(r'\(xy ([-\d.]+) ([-\d.]+)\)', fm.group(2)), float)
            if len(P) < 3:
                continue
            c = P.mean(0)
            if Path(P).contains_point(c):
                pts.append((fm.group(1), n, c[0], c[1], 'zone'))
    return nets, pts


def check(model_xml, pcb, png=None):
    prims, mesh = runviz.load_model(model_xml)
    mx, my = mesh['x'], mesh['y']
    nets, kp = kicad_points(pcb)
    polys = [p for p in prims if p['kind'] == 'poly' and p['axis'] == 2 and p['typ'] in ('Metal', 'ConductingSheet')]
    gx, gy = np.meshgrid(mx, my, indexing='ij'); G = np.c_[gx.ravel(), gy.ravel()]
    report = dict(shorts=[], opens=[], layers={})
    labs = {}
    for L, prop in LAYER_PROP.items():
        P = [p for p in polys if p['prop'] == prop]
        if not P:
            continue
        ins = np.zeros(len(G), bool)
        for p in P:
            ins |= Path(p['pts']).contains_points(G, radius=1e-7)
        ins = ins.reshape(gx.shape)
        lab, n = ndimage.label(ins)
        labs[L] = lab
        owners = {}                                                 # island -> nets
        per_net = {}
        for (l, net, x, y, kind) in kp:
            if l != L:
                continue
            y_ = -y                                                 # KiCad y down, model y up
            i, j = np.argmin(abs(mx - x)), np.argmin(abs(my - y_))
            k = lab[i, j]
            if k == 0:                                              # snap: nearest metal node within one cell
                ii = slice(max(i - 1, 0), i + 2); jj = slice(max(j - 1, 0), j + 2)
                w = lab[ii, jj][lab[ii, jj] > 0]
                k = int(w[0]) if len(w) else 0
            if k == 0:
                continue
            owners.setdefault(k, set()).add(net)
            per_net.setdefault(net, set()).add(k)
        for k, ns in owners.items():
            ns = {x for x in ns if x != '0'}
            if len(ns) > 1:
                ii, jj = np.where(lab == k)
                report['shorts'].append(dict(layer=L, nets=sorted(nets.get(x, x) for x in ns), island=int(k),
                                             bbox=[float(mx[ii].min()), float(mx[ii].max()), float(-my[jj].max()), float(-my[jj].min())]))
        for net, ks in per_net.items():
            if net != '0' and nets.get(net, '') != 'GND' and len(ks) > 1:
                report['opens'].append(dict(layer=L, net=nets.get(net, net), islands=len(ks)))
        report['layers'][L] = int(n)
    if png:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        Ls = list(labs)
        fig, axs = plt.subplots(1, len(Ls), figsize=(4.2 * len(Ls), 4.2), squeeze=False)
        bad = {(s_['layer'], s_['island']) for s_ in report['shorts']}
        for ax, L in zip(axs[0], Ls):
            lab = labs[L]
            img = np.where(lab > 0, lab % 9 + 1, 0).astype(float)
            for (l, k) in bad:
                if l == L:
                    img[lab == k] = 11
            ax.pcolormesh(mx, -my, img.T, cmap='tab20', shading='nearest', vmin=0, vmax=12)
            ax.set_aspect('equal'); ax.invert_yaxis(); ax.set_title(L + (' SHORT' if any(l == L for l, _ in bad) else ''), color='red' if any(l == L for l, _ in bad) else 'k')
        fig.suptitle('staircase copper islands per layer (KiCad coordinates); red title = nets shorted by the mesh')
        fig.tight_layout(); fig.savefig(png, dpi=110); plt.close(fig)
    return report


if __name__ == '__main__':
    png = sys.argv[sys.argv.index('--png') + 1] if '--png' in sys.argv else None
    r = check(sys.argv[1], sys.argv[2], png)
    print(json.dumps(r, indent=1))
    sys.exit(1 if r['shorts'] else 0)
