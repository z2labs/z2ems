"""Light KiCad PCB reader for routing analysis: footprints (ref, value, position, pad nets),
track segments, arcs and vias per net. Pure Python, no KiCad needed."""
import re, math
from collections import defaultdict


def _tok(s):
    return re.findall(r'\(|\)|"(?:[^"\\]|\\.)*"|[^\s()]+', s)


def parse_sexpr(s):
    stack, cur = [], []
    for t in _tok(s):
        if t == '(':
            stack.append(cur); cur = []
        elif t == ')':
            done = cur; cur = stack.pop(); cur.append(done)
        else:
            cur.append(t[1:-1] if t.startswith('"') else t)
    return cur[0]


def track_layer(t):
    """ (layer "F.Cu"), or KiCad 10 tracks with a mask opening: (layers "F.Cu" "F.Mask") """
    l = first(t, 'layer')
    if l is not None:
        return l[1]
    ls = first(t, 'layers')
    cu = [x for x in (ls[1:] if ls else []) if isinstance(x, str) and x.endswith('.Cu')]
    return cu[0] if cu else 'None'


def find(node, key):
    return [c for c in node[1:] if isinstance(c, list) and c and c[0] == key]


def first(node, key, default=None):
    r = find(node, key)
    return r[0] if r else default


class Board:
    def __init__(self, path):
        self.root = parse_sexpr(open(path, encoding='utf-8').read())
        self.nets = {}
        for n in find(self.root, 'net'):
            self.nets[n[1]] = n[2] if len(n) > 2 else n[1]
        self.footprints = {}
        self.padgeo = []
        for fp in find(self.root, 'footprint'):
            props = {p[1]: p[2] for p in find(fp, 'property')}
            at = first(fp, 'at')
            pads = {}
            fx, fy = float(at[1]), float(at[2])
            frot = float(at[3]) if len(at) > 3 else 0.0
            for p in find(fp, 'pad'):
                net = first(p, 'net')
                if net is not None:
                    pads.setdefault(p[1], self._netname(net))
                    pat = first(p, 'at'); size = first(p, 'size')
                    px, py = float(pat[1]), float(pat[2])
                    # KiCad: footprint rotation is counter-clockwise in a y-down frame
                    c, s_ = math.cos(math.radians(-frot)), math.sin(math.radians(-frot))
                    ax, ay = fx + px * c - py * s_, fy + px * s_ + py * c
                    lay = first(p, 'layers')[1:]
                    self.padgeo.append(dict(ref=props.get('Reference', '?'), pad=p[1], net=self._netname(net), x=ax, y=ay,
                                            r=0.5 * max(float(size[1]), float(size[2])), layers=lay, kind=p[2]))
            ref = props.get('Reference', '?')
            self.footprints[ref] = dict(ref=ref, value=props.get('Value', ''), lib=fp[1], layer=first(fp, 'layer')[1],
                                        at=(float(at[1]), float(at[2])), pads=pads, sheet=first(fp, 'sheetname', [None, ''])[1])
        self.segments = defaultdict(list)     # net -> [(layer, length, width)]
        self.edges = defaultdict(list)        # net -> [(layer, p0, p1, length, width, arc_mid or None)]
        self.vias = defaultdict(list)          # net -> [(x, y, layers)]
        for seg in find(self.root, 'segment'):
            a, b = first(seg, 'start'), first(seg, 'end')
            L = math.hypot(float(b[1]) - float(a[1]), float(b[2]) - float(a[2]))
            n = self._netname(first(seg, 'net'))
            self.segments[n].append((track_layer(seg), L, float(first(seg, 'width')[1])))
            self.edges[n].append((track_layer(seg), (float(a[1]), float(a[2])), (float(b[1]), float(b[2])), L, float(first(seg, 'width')[1]), None))
        for arc in find(self.root, 'arc'):
            if first(arc, 'net') is None:
                continue
            p = [first(arc, k) for k in ('start', 'mid', 'end')]
            p = [(float(q[1]), float(q[2])) for q in p]
            L = _arc_len(*p)
            n = self._netname(first(arc, 'net'))
            self.segments[n].append((track_layer(arc), L, float(first(arc, 'width')[1])))
            self.edges[n].append((track_layer(arc), p[0], p[2], L, float(first(arc, 'width')[1]), p[1]))
        for v in find(self.root, 'via'):
            at = first(v, 'at'); ly = first(v, 'layers')
            self.vias[self._netname(first(v, 'net'))].append((float(at[1]), float(at[2]), ly[1:]))

    def _netname(self, net):
        if net is None:
            return ''
        v = net[1]
        if len(net) > 2:
            return net[2]
        return self.nets.get(v, v)

    def stackup(self):
        st = first(first(self.root, 'setup'), 'stackup')
        out = []
        for l in find(st, 'layer'):
            d = {c[0]: c[1] for c in l[2:] if isinstance(c, list) and len(c) > 1}
            out.append(dict(name=l[1], **d))
        return out


def _arc_len(a, m, b):
    ax, ay = a; mx, my = m; bx, by = b
    d = 2 * (ax * (my - by) + mx * (by - ay) + bx * (ay - my))
    if abs(d) < 1e-12:
        return math.hypot(bx - ax, by - ay)
    ux = ((ax*ax + ay*ay) * (my - by) + (mx*mx + my*my) * (by - ay) + (bx*bx + by*by) * (ay - my)) / d
    uy = ((ax*ax + ay*ay) * (bx - mx) + (mx*mx + my*my) * (ax - bx) + (bx*bx + by*by) * (mx - ax)) / d
    r = math.hypot(ax - ux, ay - uy)
    t = [math.atan2(p[1] - uy, p[0] - ux) for p in (a, m, b)]
    def ccw(x, y): return (y - x) % (2 * math.pi)
    sweep = ccw(t[0], t[2])
    if ccw(t[0], t[1]) > sweep:
        sweep = 2 * math.pi - sweep
    return r * sweep


def arc_points(a, m, b, n=16):
    """polyline approximation of the arc a-m-b (n pieces)"""
    ax, ay = a; mx, my = m; bx, by = b
    d = 2 * (ax * (my - by) + mx * (by - ay) + bx * (ay - my))
    if abs(d) < 1e-12:
        return [a, b]
    ux = ((ax*ax + ay*ay) * (my - by) + (mx*mx + my*my) * (by - ay) + (bx*bx + by*by) * (ay - my)) / d
    uy = ((ax*ax + ay*ay) * (bx - mx) + (mx*mx + my*my) * (ax - bx) + (bx*bx + by*by) * (mx - ax)) / d
    r = math.hypot(ax - ux, ay - uy)
    t0, tm, t1 = [math.atan2(p[1] - uy, p[0] - ux) for p in (a, m, b)]
    def ccw(x, y): return (y - x) % (2 * math.pi)
    sweep = ccw(t0, t1)
    if ccw(t0, tm) > sweep:
        sweep = sweep - 2 * math.pi
    pts = [(ux + r * math.cos(t0 + sweep * k / n), uy + r * math.sin(t0 + sweep * k / n)) for k in range(n + 1)]
    pts[0], pts[-1] = a, b
    return pts
