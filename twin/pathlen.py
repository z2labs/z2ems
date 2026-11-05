"""Pad-to-pad routed path length, delay and phase on a KiCad board (fast, geometric).

Builds a per-net copper graph (track/arc ends, vias, pads), finds the shortest copper path between two pads
and converts it to delay with a per-layer effective permittivity (Hammerstad-Jensen microstrip on outer
layers, stripline eps_r inside). Zones/pours are not followed: a path that only closes through a pour is
reported as unrouted.
"""
import math, heapq
from collections import defaultdict

C0 = 299792458.0


def eps_eff_msl(w, h, er, t=0.035):
    if w <= 0 or h <= 0:
        return er
    u = w / h
    if t > 0:   # thickness correction (Hammerstad)
        du = t / (math.pi * h) * math.log(1 + 4 * math.e / (t / h * (1 / math.tanh(math.sqrt(6.517 * u)))**2))
        u += 0.5 * (1 + 1 / math.cosh(math.sqrt(er - 1))) * du
    a = 1 + 1/49 * math.log((u**4 + (u/52)**2) / (u**4 + 0.432)) + 1/18.7 * math.log(1 + (u/18.1)**3)
    b = 0.564 * ((er - 0.9) / (er + 3))**0.053
    return (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u)**(-a * b)


class Stack:
    def __init__(self, stackup):
        self.cu = [l['name'] for l in stackup if l.get('type') == 'copper']
        self.z = {}
        self.h_below = {}
        self.er = {}
        z = 0.0
        prev = None
        diel = []
        for l in stackup:
            t = float(l.get('thickness', 0) or 0)
            if l.get('type') == 'copper':
                self.z[l['name']] = z
                prev = l['name']
            elif 'dielectric' in l['name'] or l.get('type') in ('core', 'prepreg'):
                diel.append((z, z + t, float(l.get('epsilon_r', 4.5))))
            z += t
        self.diel = diel
        self.thickness = z

    def er_near(self, layer):
        z = self.z[layer]
        best = min(self.diel, key=lambda d: min(abs(d[0] - z), abs(d[1] - z)))
        return best[2]

    def h_to_plane(self, layer):
        """dielectric height to the adjacent inner copper layer (assumed reference plane)"""
        i = self.cu.index(layer)
        j = i + 1 if i == 0 else i - 1
        return abs(self.z[self.cu[j]] - self.z[layer]) - 0.035

    def eps_eff(self, layer, w):
        er = self.er_near(layer)
        if layer in (self.cu[0], self.cu[-1]):
            return eps_eff_msl(w, self.h_to_plane(layer), er)
        return er

    def via_span(self, la, lb):
        return abs(self.z[la] - self.z[lb])


def _key(layer, p):
    return (layer, round(p[0], 3), round(p[1], 3))


def _proj(p, a, b):
    ax, ay = a; bx, by = b
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((p[0] - ax) * dx + (p[1] - ay) * dy) / L2))
    q = (ax + t * dx, ay + t * dy)
    return t, q, math.hypot(p[0] - q[0], p[1] - q[1])


def net_graph(board, net, stack, gap_tol=0.0):
    from kipcb_lite import arc_points
    g = defaultdict(list)
    pts = defaultdict(list)                   # layer -> [(x, y, key)]
    def add(a, b, L, delay, info):
        g[a].append((b, L, delay, info)); g[b].append((a, L, delay, info))
    # pieces: straight copper pieces (arcs as polylines) with their own nodes
    pieces = []
    ends = defaultdict(list)
    for (layer, p0, p1, L, w, mid) in board.edges.get(net, []):
        poly = arc_points(p0, mid, p1) if mid is not None else [p0, p1]
        ee = stack.eps_eff(layer, w)
        keys = [_key(layer, q) for q in poly]
        for (qa, ka), (qb, kb) in zip(zip(poly, keys), zip(poly[1:], keys[1:])):
            l = math.hypot(qb[0] - qa[0], qb[1] - qa[1])
            pieces.append((layer, qa, qb, ka, kb, w, ee))
            add(ka, kb, l, l * math.sqrt(ee) / C0 * 1e-3, ('trk', layer))
        pts[layer] += [(p0[0], p0[1], keys[0]), (p1[0], p1[1], keys[-1])]
        ends[layer] += [(p0[0], p0[1], 0.5 * w, keys[0]), (p1[0], p1[1], 0.5 * w, keys[-1])]
    # copper overlap: track end touching another track end or body (T junction)
    for layer, lst in ends.items():
        for (xi, yi, ri, ki) in lst:
            for (ly, qa, qb, ka, kb, w, ee) in pieces:
                if ly != layer or ki in (ka, kb):
                    continue
                t, q, d = _proj((xi, yi), qa, qb)
                rr = ri + 0.5 * w
                if d <= rr + gap_tol:
                    kq = _key(layer, q)
                    info = ('trk', layer) if d <= rr else ('gap', layer, (round(xi, 3), round(yi, 3)), (round(q[0], 3), round(q[1], 3)), round(d - rr, 3))
                    if kq != ki:
                        add(ki, kq, d, d * math.sqrt(ee) / C0 * 1e-3, info)
                    if kq not in (ka, kb):
                        la, lb = math.hypot(q[0] - qa[0], q[1] - qa[1]), math.hypot(qb[0] - q[0], qb[1] - q[1])
                        add(ka, kq, la, la * math.sqrt(ee) / C0 * 1e-3, ('trk', layer))
                        add(kq, kb, lb, lb * math.sqrt(ee) / C0 * 1e-3, ('trk', layer))
    def near(layer, x, y, r):
        return [k for (px, py, k) in pts.get(layer, []) if math.hypot(px - x, py - y) <= r + 1e-3]
    for (x, y, lys) in board.vias.get(net, []):
        i0, i1 = stack.cu.index(lys[0]), stack.cu.index(lys[1])
        span = stack.cu[min(i0, i1):max(i0, i1) + 1]
        nodes = []
        for l in span:
            k = ('via', round(x, 3), round(y, 3), l)
            nodes.append((l, k))
            for e in near(l, x, y, 0.3):
                add(k, e, math.hypot(e[1] - x, e[2] - y), 0.0, ('pad', l))
        for (la, ka), (lb, kb) in zip(nodes, nodes[1:]):
            dz = stack.via_span(la, lb)
            add(ka, kb, dz, dz * math.sqrt(stack.er_near(la)) / C0 * 1e-3, ('via', la + '>' + lb, (round(x, 3), round(y, 3))))
    pads = {}
    for p in board.padgeo:
        if p['net'] != net:
            continue
        k = ('pad', p['ref'], p['pad'])
        pads[(p['ref'], p['pad'])] = k
        lys = stack.cu if any(l == '*.Cu' for l in p['layers']) else [l for l in p['layers'] if l in stack.z]
        for l in lys:
            for e in near(l, p['x'], p['y'], p['r']):
                d = math.hypot(e[1] - p['x'], e[2] - p['y'])
                add(k, e, d, d * math.sqrt(stack.eps_eff(l, 2 * p['r'])) / C0 * 1e-3, ('pad', l))
    return g, pads


def path(board, net, a, b, stack, gap_tol=0.0):
    """shortest copper path between pads a=(ref, pad) and b on one net: dict(length_mm, delay_ps, vias, layers)"""
    g, pads = net_graph(board, net, stack, gap_tol)
    if a not in pads or b not in pads:
        return None
    s, t = pads[a], pads[b]
    dist = {s: (0.0, 0.0)}
    prev = {}
    q = [(0.0, 0.0, 0, s)]
    cnt = 1
    while q:
        d, L, _, u = heapq.heappop(q)
        if u == t:
            break
        if d > dist[u][0]:
            continue
        for v, l, dl, info in g[u]:
            nd = d + dl
            if v not in dist or nd < dist[v][0] - 1e-18:
                dist[v] = (nd, L + l); prev[v] = (u, info, l)
                heapq.heappush(q, (nd, L + l, cnt, v)); cnt += 1
    if t not in dist:
        return dict(length_mm=None, delay_ps=None, vias=None, layers={}, routed=False)
    layers, vias, gaps = defaultdict(float), set(), []
    padxy = {('pad', q['ref'], q['pad']): (q['x'], q['y']) for q in board.padgeo if q['net'] == net}
    def xy(k):
        if k[0] == 'pad':
            return padxy.get(k)
        if k[0] == 'via':
            return (k[1] / 1.0, k[2] / 1.0)
        return (k[1], k[2])
    pieces = []
    v = t
    while v != s:
        u, info, l = prev[v]
        pieces.append((info[1] if info[0] in ('trk', 'gap') else info[0], xy(u), xy(v)))
        if info[0] == 'trk':
            layers[info[1]] += l
        elif info[0] == 'via':
            vias.add(info[2])
        elif info[0] == 'gap':
            gaps.append(dict(layer=info[1], a=info[2], b=info[3], gap_mm=info[4]))
            layers[info[1]] += l
        v = u
    return dict(length_mm=dist[t][1], delay_ps=dist[t][0] * 1e12, vias=len(vias), layers=dict(layers), routed=not gaps, gaps=gaps, pieces=pieces[::-1])
