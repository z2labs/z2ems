"""Silabs WES0110-01 Gerbers -> per-board copper geometry JSON (reference path, no PADS importer).
Parses RS-274X (D01/D02/D03, circular/rect apertures, G02/G03 arcs, G36/37 regions), builds the copper union
per layer with shapely, crops to each board's outline (windows from split_panel), reads the Excellon drills.
Output: kicad/gerber_ref/<board>.json  {outline, F_Cu:[polys], B_Cu:[polys], drills:[x,y,d]}
Coordinates are in the split boards' KiCad frame (board bbox top-left at (10,10), y down), so the SMA
footprint of the split board can be dropped in unchanged.
"""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import re, math, json, os
from shapely.geometry import Polygon, Point, LineString, box
from shapely.ops import unary_union
from shapely import affinity
G = os.path.join(ROOT, 'MSC-AMS915/WES0110-01-gerbers')
OUT = os.path.join(ROOT, 'kicad/gerber_ref'); os.makedirs(OUT, exist_ok=True)
OX, OY = 201.903822, 223.514406          # PADS mm -> KiCad panel frame: x-OX, OY-y
BOARDS = {   # KiCad panel-frame windows recovered from the split boards (edge bbox +0.05 mm line width)
 'WES0111-01-APL915M_medium-ILA': ('rect', -201.9, 111.51, -161.9, 156.51),
 'WES0112-01-ACM915D_ceramic': ('rect', -201.9, 166.51, -161.9, 223.51),
 'WES0113-01-AWH915S_small-helical': ('rect', -151.9, 133.51, -111.9, 173.51),
 'WES0114-01-AWH915M_medium-helical': ('rect', -151.9, 183.51, -111.9, 223.51),
 'WES0115-01-APF915P_panic-IFA': ('circle', -86.9, 168.51, 15.0),
 'WES0116-01-APL915P_panic-ILA': ('circle', -86.95, 208.52, 15.0),
 'WES0117-01-APN915D_meander': ('rect', -151.9, 111.51, -63.9, 122.01),
 'WES0118-01-APL915S_small-ILA': ('rect', -98.9, 133.01, -74.9, 143.51),
 'WES0119-01-APB915D_BIFA': ('rect', -60.9, 133.01, -20.9, 195.01),
}

def parse(path):
    txt = open(path, errors='replace').read()
    fmt = re.search(r'%FSLAX(\d)(\d)Y(\d)(\d)', txt); xd = int(fmt.group(2))
    unit = 25.4 if '%MOIN' in txt else 1.0; sc = unit / 10 ** xd
    aps = {}
    for m in re.finditer(r'%ADD(\d+)([CRO]),([\d.]+)(?:X([\d.]+))?', txt):
        aps[int(m.group(1))] = (m.group(2), float(m.group(3)) * unit, float(m.group(4)) * unit if m.group(4) else None)
    body = re.sub(r'%[^%]*%', '', txt).replace('\r', '').replace('\n', '')
    x = y = 0.0; ap = None; mode = 1; region = None; lastop = 2; geoms = []; multi = False   # G74 single-quadrant is the default
    def flash(X, Y, a):
        if a[0] == 'R':
            return box(X - a[1] / 2, Y - (a[2] or a[1]) / 2, X + a[1] / 2, Y + (a[2] or a[1]) / 2)
        return Point(X, Y).buffer(a[1] / 2, 16)
    for c in body.split('*'):
        if not c or c.startswith('G04'):
            continue
        m = re.match(r'^(?:G54)?D(\d+)$', c)
        if m and int(m.group(1)) >= 10:
            ap = int(m.group(1)); continue
        if c.startswith('G74'):
            multi = False; continue
        if c.startswith('G75'):
            multi = True; continue
        if c.startswith('G36'):
            region = []; continue
        if c.startswith('G37'):
            if region and len(region) > 2:
                geoms.append(Polygon(region).buffer(0))
            region = None; continue
        g = re.search(r'G0([123])', c)
        if g:
            mode = int(g.group(1))
        mx = re.search(r'X(-?\d+)', c); my = re.search(r'Y(-?\d+)', c); mi = re.search(r'I(-?\d+)', c); mj = re.search(r'J(-?\d+)', c)
        if not (mx or my):
            continue
        nx = float(mx.group(1)) * sc if mx else x; ny = float(my.group(1)) * sc if my else y
        d = re.search(r'D0([123])$', c); op = int(d.group(1)) if d else (1 if lastop == 1 else None)
        if d:
            lastop = op
        if region is not None:
            if op == 2 and region:
                geoms.append(Polygon(region).buffer(0)); region = []
            region.append((nx, ny)); x, y = nx, ny; continue
        if op == 1 and ap in aps:
            w = aps[ap][1]
            if mode == 1 or not (mi or mj):
                geoms.append(LineString([(x, y), (nx, ny)]).buffer(w / 2, 16))
            else:
                i = float(mi.group(1)) * sc if mi else 0; j = float(mj.group(1)) * sc if mj else 0
                def sweep(cx, cy):
                    a0 = math.atan2(y - cy, x - cx); a1 = math.atan2(ny - cy, nx - cx)
                    if mode == 2:
                        while a1 >= a0: a1 -= 2 * math.pi
                    else:
                        while a1 <= a0: a1 += 2 * math.pi
                    return a0, a1
                if multi:
                    cx, cy = x + i, y + j
                else:   # single quadrant: I, J unsigned; pick the centre whose radii agree and whose sweep is <= 90 deg
                    best = None
                    for si in (1, -1):
                        for sj in (1, -1):
                            ccx, ccy = x + si * abs(i), y + sj * abs(j)
                            r0 = math.hypot(x - ccx, y - ccy); r1 = math.hypot(nx - ccx, ny - ccy)
                            b0, b1 = sweep(ccx, ccy); sw = abs(b1 - b0)
                            score = abs(r0 - r1) + (0 if sw <= math.pi / 2 + 1e-6 else 1e3)
                            if best is None or score < best[0]: best = (score, ccx, ccy)
                    cx, cy = best[1], best[2]
                r = math.hypot(x - cx, y - cy)
                a0, a1 = sweep(cx, cy)
                n = max(4, int(abs(a1 - a0) * r * 4) + 2)
                pts = [(cx + r * math.cos(a0 + (a1 - a0) * k / n), cy + r * math.sin(a0 + (a1 - a0) * k / n)) for k in range(n + 1)]
                geoms.append(LineString(pts).buffer(w / 2, 16))
        elif op == 3 and ap in aps:
            geoms.append(flash(nx, ny, aps[ap]))
        x, y = nx, ny
    return unary_union(geoms)

def drills(path):
    out = []; tools = {}; cur = None
    txt = open(path, errors='replace').read().replace('\r', '')
    inch = 'INCH' in txt or 'M72' in txt or ('METRIC' not in txt)
    print('  drill: inch' if inch else '  drill: metric')
    for l in txt.split('\n'):
        m = re.match(r'^T(\d+)C([\d.]+)', l)
        if m:
            tools[int(m.group(1))] = float(m.group(2)) * (25.4 if inch else 1); cur = int(m.group(1)); continue
        m = re.match(r'^T(\d+)\s*$', l)
        if m:
            cur = int(m.group(1)); continue
        m = re.match(r'^X(-?[\d.]+)Y(-?[\d.]+)', l)
        if m and cur in tools:
            def cv(s):
                if '.' in s:
                    return float(s) * (25.4 if inch else 1)
                return float(s[:2] + '.' + s[2:]) * 25.4 if inch else float(s[:3] + '.' + s[3:])   # 2.4 inch, leading zeros kept
            out.append((cv(m.group(1)), cv(m.group(2)), tools[cur]))
    return out

print('parsing gerbers...')
top = parse(os.path.join(G, 'WES0110-01_top_wire.pho')); bot = parse(os.path.join(G, 'WES0110-01_bot_wire.pho'))
print(' top area %.0f mm2, bot %.0f mm2' % (top.area, bot.area))
dr = drills(os.path.join(G, 'WES0110-01_nc_drill.drl')); print(' drills', len(dr), dr[:2])
GX, GY = 328.9, 350.5          # Gerber mm -> KiCad panel frame (from the gerber/kicad raster alignment): x-GX, GY-y
topk = affinity.affine_transform(top, [1, 0, 0, -1, -GX, GY]); botk = affinity.affine_transform(bot, [1, 0, 0, -1, -GX, GY])
drk = [(x - GX, GY - y, d) for x, y, d in dr]
for name, spec in BOARDS.items():
    if spec[0] == 'rect':
        x0, y0, x1, y1 = spec[1:]; win = box(x0, y0, x1, y1); ox, oy = x0, y0
        outline = [10, 10, 10 + x1 - x0, 10 + y1 - y0]
    else:
        cx, cy, r = spec[1:]; win = Point(cx, cy).buffer(r, 64); ox, oy = cx - r, cy - r
        outline = [10 + r, 10 + r, r]
    def plist(g):
        g = affinity.translate(g.intersection(win), 10 - ox, 10 - oy)
        gs = [g] if g.geom_type == 'Polygon' else [p for p in getattr(g, 'geoms', []) if p.geom_type == 'Polygon']
        return [dict(ext=[list(c) for c in p.exterior.coords], holes=[[list(c) for c in h.coords] for h in p.interiors]) for p in gs if p.area > 0.001]
    out = dict(name=name, outline=spec[0], outline_param=outline, F_Cu=plist(topk), B_Cu=plist(botk),
               drills=[(x - ox + 10, y - oy + 10, d) for x, y, d in drk if win.contains(Point(x, y))])
    json.dump(out, open(os.path.join(OUT, name + '.json'), 'w'))
    print(' %-36s F %2d polys  B %2d polys  drills %3d  Fcu area %.0f' % (
        name, len(out['F_Cu']), len(out['B_Cu']), len(out['drills']), sum(Polygon(p['ext'], p['holes']).area for p in out['F_Cu'])))
