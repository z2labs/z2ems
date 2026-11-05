"""
kisch: KiCad 6..10 schematic (.kicad_sch) reader for the twin viewer.

Flattens the sheet hierarchy into per-instance drawings (symbols transformed to sheet coordinates,
wires, labels, sheet boxes) and resolves the net of every wire group and pin. Net names come from
the PCB (pad nets by reference and pin number) when a board is given, so a click in the schematic maps
1:1 onto the copper in the 3D view; otherwise from labels.

    python kisch.py project.kicad_sch [board.kicad_pcb] > scene.json
"""
import json
import math
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kipcb_lite import parse_sexpr, find, first          # noqa: E402

GRID = 100.0       # point matching resolution: 0.01 mm


def _f(v):
    return float(v)


def _at(node):
    a = first(node, 'at')
    if a is None:
        return 0.0, 0.0, 0.0
    return _f(a[1]), _f(a[2]), _f(a[3]) if len(a) > 3 else 0.0


def _hidden(node):
    if first(node, 'hide') is not None and (len(first(node, 'hide')) < 2 or first(node, 'hide')[1] == 'yes'):
        return True
    eff = first(node, 'effects')
    if eff is not None:
        if 'hide' in eff[1:]:
            return True
        h = first(eff, 'hide')
        if h is not None and (len(h) < 2 or h[1] == 'yes'):
            return True
    return False


def _font(node):
    eff = first(node, 'effects')
    size, just = 1.27, []
    if eff is not None:
        fnt = first(eff, 'font')
        if fnt is not None and first(fnt, 'size') is not None:
            size = _f(first(fnt, 'size')[2])
        j = first(eff, 'justify')
        if j is not None:
            just = j[1:]
    return size, just


class Xform:
    """library (y up) -> sheet (y down): flip y, rotate by the symbol angle, then mirror"""

    def __init__(self, x, y, angle, mirror):
        self.x, self.y, self.mirror = x, y, mirror
        a = math.radians(angle)
        self.c, self.s = round(math.cos(a), 9), round(math.sin(a), 9)

    def __call__(self, px, py):
        u, v = px, -py
        # KiCad orientation: rotation (counter-clockwise on screen), then mirror about the symbol axis
        u, v = u * self.c + v * self.s, -u * self.s + v * self.c
        if self.mirror == 'x':
            v = -v
        elif self.mirror == 'y':
            u = -u
        return self.x + u, self.y + v


def _pts(node):
    p = first(node, 'pts')
    return [(_f(q[1]), _f(q[2])) for q in find(p, 'xy')] if p is not None else []


def _fill(node):
    f = first(node, 'fill')
    t = first(f, 'type') if f is not None else None
    return t[1] if t is not None else 'none'


class LibSymbol:
    def __init__(self, node, libs):
        self.name = node[1]
        self.power = first(node, 'power') is not None
        self.units = defaultdict(list)       # (unit, style) -> graphics/pins
        ext = first(node, 'extends')
        if ext is not None and ext[1] in libs:
            self.units = libs[ext[1]].units
        for sub in find(node, 'symbol'):
            m = re.match(r'.*_(\d+)_(\d+)$', sub[1])
            key = (int(m.group(1)), int(m.group(2))) if m else (0, 1)
            self.units[key] += [c for c in sub[1:] if isinstance(c, list)]

    def items(self, unit, style=1):
        out = []
        for (u, s), its in self.units.items():
            if u in (0, unit) and s in (0, style):
                out += its
        return out


def _graphics(items, X):
    prims, pins = [], []
    for it in items:
        k = it[0]
        if k == 'polyline':
            prims.append(dict(k='pl', p=[X(*q) for q in _pts(it)], f=_fill(it)))
        elif k == 'rectangle':
            a, b = first(it, 'start'), first(it, 'end')
            c = [X(_f(a[1]), _f(a[2])), X(_f(b[1]), _f(a[2])), X(_f(b[1]), _f(b[2])), X(_f(a[1]), _f(b[2]))]
            prims.append(dict(k='pl', p=c + [c[0]], f=_fill(it)))
        elif k == 'circle':
            c, r = first(it, 'center'), _f(first(it, 'radius')[1])
            prims.append(dict(k='c', p=X(_f(c[1]), _f(c[2])), r=r, f=_fill(it)))
        elif k == 'arc':
            a, m, b = first(it, 'start'), first(it, 'mid'), first(it, 'end')
            if m is not None:
                prims.append(dict(k='a', p=[X(_f(a[1]), _f(a[2])), X(_f(m[1]), _f(m[2])), X(_f(b[1]), _f(b[2]))], f=_fill(it)))
        elif k == 'text':
            x, y, ang = _at(it)
            size, _ = _font(it)
            if not _hidden(it):
                prims.append(dict(k='t', t=it[1], p=X(x, y), s=size))
        elif k == 'pin':
            x, y, ang = _at(it)
            L = _f(first(it, 'length')[1]) if first(it, 'length') is not None else 2.54
            a = math.radians(ang)
            ex, ey = x + L * math.cos(a), y + L * math.sin(a)
            num = first(it, 'number')[1] if first(it, 'number') is not None else ''
            nm = first(it, 'name')[1] if first(it, 'name') is not None else ''
            hide = 'hide' in it[1:] or (first(it, 'hide') is not None)
            pins.append(dict(n=num, name=nm, p=X(x, y), e=X(ex, ey), hide=hide, type=it[1]))
    return prims, pins


def _bbox(pts):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return [min(xs), min(ys), max(xs), max(ys)] if xs else None


class UF:
    def __init__(self):
        self.p = {}

    def f(self, a):
        self.p.setdefault(a, a)
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def u(self, a, b):
        a, b = self.f(a), self.f(b)
        if a != b:
            self.p[a] = b


def _k(p):
    return (round(p[0] * GRID), round(p[1] * GRID))


def _on_seg(p, a, b, tol=0.01):
    (px, py), (ax, ay), (bx, by) = p, a, b
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    if L2 == 0:
        return False
    t = ((px - ax) * dx + (py - ay) * dy) / L2
    if t <= 1e-6 or t >= 1 - 1e-6:
        return False
    return math.hypot(ax + t * dx - px, ay + t * dy - py) <= tol


class Schematic:
    def __init__(self, root_file, pad_nets=None):
        self.dir = os.path.dirname(os.path.abspath(root_file))
        self.pad_nets = pad_nets or {}          # (ref, pin) -> net
        self.files = {}
        self.sheets = []
        root = self._load(os.path.basename(root_file))
        ruid = first(root, 'uuid')[1]
        self._walk(os.path.basename(root_file), '/' + ruid, 'Root', None, '1')

    def _load(self, name):
        if name not in self.files:
            self.files[name] = parse_sexpr(open(os.path.join(self.dir, name), encoding='utf-8').read())
        return self.files[name]

    def _walk(self, fname, path, name, parent, page):
        node = self._load(fname)
        idx = len(self.sheets)
        sh = dict(path=path, name=name, file=fname, parent=parent, page=page, children=[])
        self.sheets.append(sh)
        self._draw(node, sh)
        for s in find(node, 'sheet'):
            props = {p[1]: p[2] for p in find(s, 'property')}
            sname = props.get('Sheetname', props.get('Sheet name', '?'))
            sfile = props.get('Sheetfile', props.get('Sheet file', ''))
            cpath = path + '/' + first(s, 'uuid')[1]
            cpage = ''
            for inst in find(s, 'instances'):
                for proj in find(inst, 'project'):
                    for pth in find(proj, 'path'):
                        if pth[1] == path and first(pth, 'page') is not None:
                            cpage = first(pth, 'page')[1]
            if not os.path.exists(os.path.join(self.dir, sfile)):
                continue
            sh['children'].append(len(self.sheets))
            self._walk(sfile, cpath, sname, idx, cpage)

    def _ref(self, sym, path):
        for inst in find(sym, 'instances'):
            for proj in find(inst, 'project'):
                for pth in find(proj, 'path'):
                    if pth[1] == path:
                        r = first(pth, 'reference')
                        u = first(pth, 'unit')
                        return (r[1] if r is not None else None), (int(u[1]) if u is not None else None)
        props = {p[1]: p[2] for p in find(sym, 'property')}
        return props.get('Reference'), None

    def _draw(self, node, sh):
        libs = {}
        ls = first(node, 'lib_symbols')
        if ls is not None:
            for s in find(ls, 'symbol'):
                libs[s[1]] = LibSymbol(s, libs)
        paper = first(node, 'paper')
        sizes = dict(A4=(297, 210), A3=(420, 297), A2=(594, 420), A1=(841, 594), A0=(1189, 841), A=(279.4, 215.9),
                     B=(431.8, 279.4), C=(558.8, 431.8), D=(863.6, 558.8), E=(1117.6, 863.6))
        w, h = sizes.get(paper[1], (297, 210)) if paper is not None else (297, 210)
        if paper is not None and paper[1] == 'User' and len(paper) > 3:
            w, h = _f(paper[2]), _f(paper[3])
        if paper is not None and 'portrait' in paper[1:]:
            w, h = h, w
        sh['paper'] = [w, h]
        wires = [_pts(wn)[:2] for wn in find(node, 'wire')]
        buses = [_pts(wn)[:2] for wn in find(node, 'bus')]
        lines = [_pts(wn) for wn in find(node, 'polyline')]
        junctions = [list(_at(j)[:2]) for j in find(node, 'junction')]
        noconn = [list(_at(j)[:2]) for j in find(node, 'no_connect')]
        labels = []
        for kind in ('label', 'global_label', 'hierarchical_label'):
            for lb in find(node, kind):
                x, y, a = _at(lb)
                size, just = _font(lb)
                shape = first(lb, 'shape')
                labels.append(dict(t=lb[1], p=[x, y], a=a, k=kind[0], s=size, j=just, shape=shape[1] if shape else ''))
        texts = []
        for tx in find(node, 'text'):
            x, y, a = _at(tx)
            size, just = _font(tx)
            texts.append(dict(t=tx[1], p=[x, y], a=a, s=size, j=just))
        symbols = []
        for sym in find(node, 'symbol'):
            lid = first(sym, 'lib_id')
            if lid is None:
                continue
            lib = libs.get(first(sym, 'lib_name')[1] if first(sym, 'lib_name') is not None else lid[1]) or libs.get(lid[1])
            x, y, a = _at(sym)
            mir = first(sym, 'mirror')
            unit = int(first(sym, 'unit')[1]) if first(sym, 'unit') is not None else 1
            style = int(first(sym, 'body_style')[1]) if first(sym, 'body_style') is not None else 1
            ref, iu = self._ref(sym, sh['path'])
            if iu:
                unit = iu
            X = Xform(x, y, a, mir[1] if mir is not None else None)
            prims, pins = _graphics(lib.items(unit, style), X) if lib else ([], [])
            props = []
            pv = {}
            for p in find(sym, 'property'):
                pv[p[1]] = p[2]
                if p[1] in ('Reference', 'Value') and not _hidden(p):
                    px, py, pa = _at(p)
                    size, just = _font(p)
                    txt = ref if p[1] == 'Reference' and ref else p[2]
                    props.append(dict(t=txt, p=[px, py], a=pa, s=size, j=just, k=p[1][0]))
            allp = [q for pr in prims for q in (pr['p'] if pr['k'] in ('pl', 'a') else [pr['p']])] + [q for pn in pins for q in (pn['p'], pn['e'])]
            power = bool(lib and lib.power) or (ref or '').startswith('#')
            symbols.append(dict(ref=ref, uid=first(sym, 'uuid')[1] if first(sym, 'uuid') is not None else '', value=pv.get('Value', ''), lib=lid[1], unit=unit, power=power,
                                fp=pv.get('Footprint', ''), dnp=first(sym, 'dnp') is not None and first(sym, 'dnp')[1] == 'yes',
                                prims=prims, pins=pins, props=props, bbox=_bbox(allp)))
        boxes = []
        for s in find(node, 'sheet'):
            x, y, _ = _at(s)
            sz = first(s, 'size')
            props = {p[1]: p[2] for p in find(s, 'property')}
            spins = []
            for pn in find(s, 'pin'):
                px, py, pa = _at(pn)
                spins.append(dict(t=pn[1], type=pn[2], p=[px, py], a=pa))
            boxes.append(dict(p=[x, y], s=[_f(sz[1]), _f(sz[2])], name=props.get('Sheetname', props.get('Sheet name', '')),
                              file=props.get('Sheetfile', props.get('Sheet file', '')), uuid=first(s, 'uuid')[1],
                              path=sh['path'] + '/' + first(s, 'uuid')[1], pins=spins))
        # ---- connectivity inside this sheet instance
        uf = UF()
        for a, b in wires:
            uf.u(_k(a), _k(b))
        ends = [q for wv in wires for q in wv]
        anchor = ends + [tuple(j) for j in junctions] + [tuple(pn['p']) for s in symbols for pn in s['pins']] + \
            [tuple(l['p']) for l in labels] + [tuple(pn['p']) for bx in boxes for pn in bx['pins']]
        for q in anchor:
            for a, b in wires:
                if _on_seg(q, a, b):
                    uf.u(_k(q), _k(a))
        names = defaultdict(set)
        for s in symbols:
            for pn in s['pins']:
                k = uf.f(_k(pn['p']))
                net = None
                if s['power']:
                    net = s['value']
                elif s['ref'] and (s['ref'], pn['n']) in self.pad_nets:
                    net = self.pad_nets[(s['ref'], pn['n'])]
                if net:
                    names[k].add(net)
        for l in labels:
            k = uf.f(_k(l['p']))
            if not names[k]:
                names[k].add(('/' + sh['name'] + '/' if l['k'] == 'h' else '') + l['t'])
        def net_of(p):
            n = names.get(uf.f(_k(p)))
            return sorted(n)[0] if n else ''
        for s in symbols:
            for pn in s['pins']:
                if not s['power'] and s['ref'] and (s['ref'], pn['n']) in self.pad_nets:
                    pn['net'] = self.pad_nets[(s['ref'], pn['n'])]
                else:
                    pn['net'] = net_of(pn['p'])
        sh['wires'] = [dict(p=[list(a), list(b)], net=net_of(a)) for a, b in wires]
        sh['buses'] = [[list(a), list(b)] for a, b in buses]
        sh['lines'] = [[list(q) for q in l] for l in lines]
        sh['junctions'] = junctions
        sh['noconn'] = noconn
        for l in labels:
            l['net'] = net_of(l['p'])
        sh['labels'] = labels
        sh['texts'] = texts
        for bx in boxes:
            for pn in bx['pins']:
                pn['net'] = net_of(pn['p'])
        sh['boxes'] = boxes
        sh['symbols'] = symbols

    def scene(self):
        return _round(dict(version=1, sheets=self.sheets))


def _round(o):
    if isinstance(o, float):
        return round(o, 3)
    if isinstance(o, (list, tuple)):
        return [_round(v) for v in o]
    if isinstance(o, dict):
        return {k: _round(v) for k, v in o.items()}
    return o


def pad_nets_from_board(pcb):
    from kipcb_lite import Board
    b = Board(pcb)
    return {(r, p): n for r, f in b.footprints.items() for p, n in f['pads'].items()}


def find_root(folder_or_pcb):
    """root schematic of a project: <stem>.kicad_sch next to the board"""
    if folder_or_pcb.endswith('.kicad_pcb'):
        cand = folder_or_pcb[:-len('.kicad_pcb')] + '.kicad_sch'
        if os.path.exists(cand):
            return cand
        folder_or_pcb = os.path.dirname(folder_or_pcb)
    pros = [f for f in os.listdir(folder_or_pcb) if f.endswith('.kicad_pro')]
    for p in pros:
        c = os.path.join(folder_or_pcb, p[:-4] + '_sch')
        if os.path.exists(c):
            return c
    return None


def stats(sc):
    """connectivity sanity: share of symbol pins that touch a wire, label or other pin"""
    tot = hit = 0
    for sh in sc.sheets:
        pts = set()
        for w in sh['wires']:
            for q in w['p']:
                pts.add(_k(q))
        for s in sh['symbols']:
            for pn in s['pins']:
                if s['power']:
                    continue
                tot += 1
                k = _k(pn['p'])
                if k in pts or any(_on_seg(pn['p'], w['p'][0], w['p'][1]) for w in sh['wires']):
                    hit += 1
    return hit, tot


if __name__ == '__main__':
    root = sys.argv[1]
    pn = pad_nets_from_board(sys.argv[2]) if len(sys.argv) > 2 else None
    sc = Schematic(root, pn)
    if '--stats' in sys.argv:
        h, t = stats(sc)
        print('sheets', len(sc.sheets), 'symbols', sum(len(s['symbols']) for s in sc.sheets), 'pins on wires %d/%d' % (h, t))
    else:
        json.dump(sc.scene(), sys.stdout, separators=(',', ':'))
