#!/usr/bin/env python3
""" runviz: the pictures of one simulation run, for a human to look at (z2mon.py and the KiCad plugin's result view
    both use these, so the two views grow together).

      geometry(model_xml, out_prefix)   top view, side cut on the port axis, cross-section, isometric 3D - from the CSX
                                        primitives of model.xml (what the solver meshes), with the mesh lines
      sparams(s2p_or_json, out_prefix)  |S11| |S21| (dB), Smith chart of S11, group delay / phase
      progress(stdout_txt)              dict(ts, t_s, speed, energy_dB) parsed from the engine output (for the live chart)
      progress_png(stdout_txt, out)     energy / speed vs timestep

    Only numpy + matplotlib (no CSXCAD, no browser libraries): it runs on the offline laptop as well.
"""
import glob, json, math, os, re, sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
try:                                                   # prefer the mpl_toolkits shipped with this matplotlib (mixed installs)
    import mpl_toolkits
    _mt = os.path.join(os.path.dirname(matplotlib.__file__), '..', 'mpl_toolkits')
    if os.path.isdir(_mt) and os.path.abspath(_mt) not in mpl_toolkits.__path__:
        mpl_toolkits.__path__.insert(0, os.path.abspath(_mt))
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    if '3d' not in matplotlib.projections.get_projection_names():
        matplotlib.projections.register_projection(Axes3D)
except Exception:
    import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MPoly, Rectangle, Circle

# ---------------------------------------------------------------- model.xml primitives
_F = r'(-?[\d.eE+-]+)'


def load_model(xml):
    """ primitives of model.xml: list of dict(prop, typ, kind, ...), and the mesh lines {x,y,z} """
    s = open(xml, encoding='utf-8', errors='replace').read()
    prims = []
    for m in re.finditer(r'<(Metal|Material|ConductingSheet|LumpedElement|Excitation|ProbeBox|DumpBox|LorentzMaterial|DebyeMaterial)\s+ID="\d+"\s+Name="([^"]+)"(.*?)</\1>', s, re.S):
        typ, name, body = m.group(1), m.group(2), m.group(3)
        if typ in ('ProbeBox', 'DumpBox', 'Excitation'):
            continue
        for b in re.finditer(r'<Box\s+Priority="(-?\d+)">\s*<P1 X="%s" Y="%s" Z="%s"\s*/>\s*<P2 X="%s" Y="%s" Z="%s"' % ((_F,) * 6), body):
            v = [float(x) for x in b.groups()[1:]]
            prims.append(dict(prop=name, typ=typ, kind='box', pri=int(b.group(1)), lo=np.minimum(v[:3], v[3:]), hi=np.maximum(v[:3], v[3:])))
        for c in re.finditer(r'<Cylinder\s+Priority="(-?\d+)"\s+Radius="%s">\s*<P1 X="%s" Y="%s" Z="%s"\s*/>\s*<P2 X="%s" Y="%s" Z="%s"' % ((_F,) * 7), body):
            g = [float(x) for x in c.groups()[1:]]
            prims.append(dict(prop=name, typ=typ, kind='cyl', pri=int(c.group(1)), r=g[0], s=np.array(g[1:4]), e=np.array(g[4:7])))
        for c in re.finditer(r'<CylindricalShell\s+Priority="(-?\d+)"\s+Radius="%s"\s+ShellWidth="%s">\s*<P1 X="%s" Y="%s" Z="%s"\s*/>\s*<P2 X="%s" Y="%s" Z="%s"' % ((_F,) * 8), body):
            g = [float(x) for x in c.groups()[1:]]
            prims.append(dict(prop=name, typ=typ, kind='shell', pri=int(c.group(1)), r=g[0], w=g[1], s=np.array(g[2:5]), e=np.array(g[5:8])))
        for p in re.finditer(r'<(LinPoly|Polygon)\s+Priority="(-?\d+)"\s+Elevation="%s"\s+NormDir="(\d)"\s+QtyVertices="\d+"(?:\s+Length="%s")?\s*>(.*?)</\1>' % (_F, _F), body, re.S):
            pts = np.array([[float(a), float(b_)] for a, b_ in re.findall(r'<Vertex X1="%s" X2="%s"' % (_F, _F), p.group(6))])
            el = float(p.group(3)); ln = float(p.group(5)) if p.group(5) else 0.0
            prims.append(dict(prop=name, typ=typ, kind='poly', pri=int(p.group(2)), pts=pts, axis=int(p.group(4)), z0=el, z1=el + ln))
    mesh = {}
    for ax in 'XYZ':
        m = re.search(r'<%sLines[^>]*>([^<]*)</%sLines>' % (ax, ax), s)
        mesh[ax.lower()] = np.array([float(v) for v in m.group(1).replace(',', ' ').split()]) if m else np.array([])
    return prims, mesh


def style(p):
    n, t = p['prop'].lower(), p['typ'].lower()
    if t in ('metal', 'conductingsheet'):
        if 'pin' in n or 'strip' in n:
            return '#d62728', 0.95, 'pin'
        if 'body' in n or 'flange' in n or 'leg' in n:
            return '#7f7f7f', 0.55, 'body'
        if 'hole' in n or 'via' in n:
            return '#444444', 0.9, 'via'
        return '#e6a817', 0.9, 'cu'
    if 'ptfe' in n or 'coax' in n:
        return '#1f77b4', 0.25, 'ptfe'
    if 'diel' in n or 'sub' in n or 'fr4' in n or t == 'material':
        return '#2ca02c', 0.18, 'diel'
    return '#999999', 0.4, 'other'


def _bbox(prims):
    lo = np.array([np.inf] * 3); hi = -lo
    for p in prims:
        if p['kind'] == 'box':
            lo = np.minimum(lo, p['lo']); hi = np.maximum(hi, p['hi'])
        elif p['kind'] in ('cyl', 'shell'):
            lo = np.minimum(lo, np.minimum(p['s'], p['e']) - p['r']); hi = np.maximum(hi, np.maximum(p['s'], p['e']) + p['r'])
        elif p['kind'] == 'poly' and p['axis'] == 2:
            lo = np.minimum(lo, [p['pts'][:, 0].min(), p['pts'][:, 1].min(), min(p['z0'], p['z1'])])
            hi = np.maximum(hi, [p['pts'][:, 0].max(), p['pts'][:, 1].max(), max(p['z0'], p['z1'])])
    return lo, hi


def _top(ax, prims, mesh, xr, yr, title, mesh_lines=True):
    if mesh_lines and len(mesh['x']):
        for x in mesh['x']:
            if xr[0] <= x <= xr[1]:
                ax.axvline(x, color='#000', lw=0.15, alpha=0.35)
        for y in mesh['y']:
            if yr[0] <= y <= yr[1]:
                ax.axhline(y, color='#000', lw=0.15, alpha=0.35)
    for p in sorted(prims, key=lambda q: q['pri']):
        col, al, cls = style(p)
        if cls == 'diel' or cls == 'other':
            continue
        if p['kind'] == 'box':
            ax.add_patch(Rectangle((p['lo'][0], p['lo'][1]), p['hi'][0] - p['lo'][0], p['hi'][1] - p['lo'][1], fc=col, ec='k', lw=0.2, alpha=al))
        elif p['kind'] == 'poly' and p['axis'] == 2:
            ax.add_patch(MPoly(p['pts'], closed=True, fc=col, ec='k', lw=0.2, alpha=al))
        elif p['kind'] in ('cyl', 'shell'):
            if abs(p['s'][2] - p['e'][2]) > 1e-9 and abs(p['s'][0] - p['e'][0]) < 1e-9:      # vertical: circle
                ax.add_patch(Circle((p['s'][0], p['s'][1]), p['r'], fc=col, ec='k', lw=0.2, alpha=al))
            else:                                                                              # along x: strip
                x0, x1 = sorted([p['s'][0], p['e'][0]])
                ax.add_patch(Rectangle((x0, p['s'][1] - p['r']), x1 - x0, 2 * p['r'], fc=col, ec='k', lw=0.2, alpha=al))
    ax.set_xlim(*xr); ax.set_ylim(*yr); ax.set_aspect('equal'); ax.set_title(title, fontsize=9); ax.tick_params(labelsize=7)


def _section_xz(ax, prims, mesh, y0, xr, zr, title):
    """ side cut at y = y0 """
    for x in mesh['x']:
        if xr[0] <= x <= xr[1]:
            ax.axvline(x, color='#000', lw=0.15, alpha=0.3)
    for z in mesh['z']:
        if zr[0] <= z <= zr[1]:
            ax.axhline(z, color='#000', lw=0.15, alpha=0.3)
    for p in sorted(prims, key=lambda q: q['pri']):
        col, al, cls = style(p)
        if p['kind'] == 'box' and p['lo'][1] - 1e-9 <= y0 <= p['hi'][1] + 1e-9:
            ax.add_patch(Rectangle((p['lo'][0], p['lo'][2]), p['hi'][0] - p['lo'][0], max(p['hi'][2] - p['lo'][2], 0.02), fc=col, ec='k', lw=0.2, alpha=al))
        elif p['kind'] == 'poly' and p['axis'] == 2 and _inside(p['pts'], y0, xr):
            for x0, x1 in _inside(p['pts'], y0, xr):
                ax.add_patch(Rectangle((x0, min(p['z0'], p['z1'])), x1 - x0, max(abs(p['z1'] - p['z0']), 0.02), fc=col, ec='k', lw=0.2, alpha=al))
        elif p['kind'] in ('cyl', 'shell') and abs(p['s'][1] - y0) <= p['r'] + 1e-9:
            if abs(p['s'][2] - p['e'][2]) < 1e-9:                      # along x
                h = math.sqrt(max(p['r'] ** 2 - (p['s'][1] - y0) ** 2, 0))
                x0, x1 = sorted([p['s'][0], p['e'][0]])
                ax.add_patch(Rectangle((x0, p['s'][2] - h), x1 - x0, 2 * h, fc=col, ec='k', lw=0.2, alpha=al))
            else:                                                        # vertical via
                z0, z1 = sorted([p['s'][2], p['e'][2]])
                ax.add_patch(Rectangle((p['s'][0] - p['r'], z0), 2 * p['r'], z1 - z0, fc=col, ec='k', lw=0.2, alpha=al))
    ax.set_xlim(*xr); ax.set_ylim(*zr); ax.set_aspect('equal'); ax.set_title(title, fontsize=9); ax.tick_params(labelsize=7)


def _inside(pts, y0, xr):
    """ x intervals where the horizontal line y = y0 is inside the polygon """
    xs = []
    n = len(pts)
    for i in range(n):
        (x1, y1), (x2, y2) = pts[i], pts[(i + 1) % n]
        if (y1 <= y0 < y2) or (y2 <= y0 < y1):
            xs.append(x1 + (y0 - y1) * (x2 - x1) / (y2 - y1))
    xs.sort()
    return [(max(xs[i], xr[0]), min(xs[i + 1], xr[1])) for i in range(0, len(xs) - 1, 2) if xs[i + 1] > xr[0] and xs[i] < xr[1]]


def _section_yz(ax, prims, mesh, x0, yr, zr, title):
    for y in mesh['y']:
        if yr[0] <= y <= yr[1]:
            ax.axvline(y, color='#000', lw=0.15, alpha=0.3)
    for z in mesh['z']:
        if zr[0] <= z <= zr[1]:
            ax.axhline(z, color='#000', lw=0.15, alpha=0.3)
    for p in sorted(prims, key=lambda q: q['pri']):
        col, al, cls = style(p)
        if p['kind'] == 'box' and p['lo'][0] - 1e-9 <= x0 <= p['hi'][0] + 1e-9:
            ax.add_patch(Rectangle((p['lo'][1], p['lo'][2]), p['hi'][1] - p['lo'][1], max(p['hi'][2] - p['lo'][2], 0.02), fc=col, ec='k', lw=0.2, alpha=al))
        elif p['kind'] == 'poly' and p['axis'] == 2:
            for y0, y1 in _inside(p['pts'][:, ::-1], x0, yr):
                ax.add_patch(Rectangle((y0, min(p['z0'], p['z1'])), y1 - y0, max(abs(p['z1'] - p['z0']), 0.02), fc=col, ec='k', lw=0.2, alpha=al))
        elif p['kind'] in ('cyl', 'shell'):
            if abs(p['s'][2] - p['e'][2]) < 1e-9 and min(p['s'][0], p['e'][0]) <= x0 <= max(p['s'][0], p['e'][0]):
                ax.add_patch(Circle((p['s'][1], p['s'][2]), p['r'], fc=col, ec='k', lw=0.2, alpha=al))
            elif abs(p['s'][0] - x0) <= p['r']:
                z0, z1 = sorted([p['s'][2], p['e'][2]])
                ax.add_patch(Rectangle((p['s'][1] - p['r'], z0), 2 * p['r'], z1 - z0, fc=col, ec='k', lw=0.2, alpha=al))
    ax.set_xlim(*yr); ax.set_ylim(*zr); ax.set_aspect('equal'); ax.set_title(title, fontsize=9); ax.tick_params(labelsize=7)


def _iso(ax, prims, lo, hi):
    """ isometric 3D of the metal: boxes and polygons as filled faces (matplotlib 3D) """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    faces, cols = [], []
    for p in sorted(prims, key=lambda q: q['pri']):
        col, al, cls = style(p)
        if cls in ('diel', 'other', 'ptfe'):
            continue
        if p['kind'] == 'box':
            (x0, y0, z0), (x1, y1, z1) = p['lo'], p['hi']
            v = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0], [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
            for f in ([0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4]):
                faces.append(v[f]); cols.append((*matplotlib.colors.to_rgb(col), al))
        elif p['kind'] == 'poly' and p['axis'] == 2:
            for z in (p['z0'], p['z1']) if abs(p['z1'] - p['z0']) > 1e-9 else (p['z0'],):
                faces.append(np.column_stack([p['pts'], np.full(len(p['pts']), z)])); cols.append((*matplotlib.colors.to_rgb(col), al))
        elif p['kind'] in ('cyl', 'shell'):
            n = 16; a = np.linspace(0, 2 * np.pi, n, endpoint=False)
            if abs(p['s'][2] - p['e'][2]) > 1e-9:                                     # vertical
                ring0 = np.column_stack([p['s'][0] + p['r'] * np.cos(a), p['s'][1] + p['r'] * np.sin(a), np.full(n, p['s'][2])])
                ring1 = ring0.copy(); ring1[:, 2] = p['e'][2]
            else:
                ring0 = np.column_stack([np.full(n, p['s'][0]), p['s'][1] + p['r'] * np.cos(a), p['s'][2] + p['r'] * np.sin(a)])
                ring1 = ring0.copy(); ring1[:, 0] = p['e'][0]
            for i in range(n):
                faces.append(np.array([ring0[i], ring0[(i + 1) % n], ring1[(i + 1) % n], ring1[i]])); cols.append((*matplotlib.colors.to_rgb(col), al))
    pc = Poly3DCollection(faces, facecolors=cols, edgecolors=(0, 0, 0, 0.15), linewidths=0.2)
    ax.add_collection3d(pc)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect((hi[0] - lo[0], hi[1] - lo[1], max(hi[2] - lo[2], 1.0) * 2))
    ax.view_init(elev=28, azim=-60); ax.tick_params(labelsize=6)
    ax.set_title('3D (metal), z stretched 2x', fontsize=9)


def geometry(model_xml, out_prefix, port_y=None, zoom=8.0):
    """ 4 pages: full top view, launch zoom (top + side cut + cross-section), 3D """
    prims, mesh = load_model(model_xml)
    lo, hi = _bbox([p for p in prims if style(p)[2] in ('cu', 'pin', 'body', 'via')])
    # the port axis: the pin cylinders along x
    pins = [p for p in prims if p['kind'] == 'cyl' and style(p)[2] == 'pin' and abs(p['s'][2] - p['e'][2]) < 1e-9]
    yc = port_y if port_y is not None else (pins[0]['s'][1] if pins else 0.5 * (lo[1] + hi[1]))
    zc = pins[0]['s'][2] if pins else 0.5 * (lo[2] + hi[2])
    cu = [p for p in prims if p['kind'] == 'poly' and style(p)[2] == 'cu']
    xl = min(p['pts'][:, 0].min() for p in cu) if cu else (min(p['s'][0] for p in pins) if pins else lo[0])   # board edge
    out = {}
    fig, ax = plt.subplots(figsize=(14, 4 + 10 * (hi[1] - lo[1]) / max(hi[0] - lo[0], 1)))
    _top(ax, prims, mesh, (lo[0] - 1, hi[0] + 1), (lo[1] - 1, hi[1] + 1), 'top view of the simulated model (CSX primitives + mesh lines)', mesh_lines=True)
    fig.tight_layout(); out['top'] = out_prefix + '_geom_top.png'; fig.savefig(out['top'], dpi=110); plt.close(fig)
    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    xr = (xl - 3.5, xl + zoom)
    _top(ax[0, 0], prims, mesh, xr, (yc - 4.5, yc + 4.5), 'launch zoom, top view')
    _section_xz(ax[0, 1], prims, mesh, yc, xr, (zc - 3.0, zc + 3.5), 'side cut on the pin axis (y = %.2f)' % yc)
    _section_yz(ax[1, 0], prims, mesh, xl + 1.0, (yc - 4.5, yc + 4.5), (zc - 3.0, zc + 3.5), 'cross-section 1 mm inside the edge (x = %.2f)' % (xl + 1.0))
    _section_yz(ax[1, 1], prims, mesh, xl + zoom * 0.6, (yc - 4.5, yc + 4.5), (zc - 3.0, zc + 3.5), 'cross-section on the line (x = %.2f)' % (xl + zoom * 0.6))
    fig.tight_layout(); out['launch'] = out_prefix + '_geom_launch.png'; fig.savefig(out['launch'], dpi=110); plt.close(fig)
    try:
        fig = plt.figure(figsize=(12, 7)); ax = fig.add_subplot(111, projection='3d')
        _iso(ax, prims, lo, hi)
        fig.tight_layout(); out['iso'] = out_prefix + '_geom_3d.png'; fig.savefig(out['iso'], dpi=110); plt.close(fig)
    except Exception as e:                                     # the 3D toolkit is not always importable
        out['iso_error'] = str(e)
    out['mesh'] = dict(n=[len(mesh[a]) for a in 'xyz'], cells=float(np.prod([max(len(mesh[a]) - 1, 1) for a in 'xyz'])) / 1e6,
                       min_mm=[float(np.diff(mesh[a]).min()) if len(mesh[a]) > 1 else None for a in 'xyz'])
    return out


# ---------------------------------------------------------------- S-parameters
def read_touchstone(path):
    f, S = [], []
    fmt, unit = 'RI', 1e9
    for line in open(path):
        t = line.strip()
        if not t or t.startswith('!'):
            continue
        if t.startswith('#'):
            u = t.upper().split()
            unit = {'HZ': 1, 'KHZ': 1e3, 'MHZ': 1e6, 'GHZ': 1e9}[u[1]]
            fmt = u[3]
            continue
        v = [float(x) for x in t.split()]
        n = (len(v) - 1) // 2
        c = []
        for i in range(n):
            a, b = v[1 + 2 * i], v[2 + 2 * i]
            if fmt == 'RI':
                c.append(a + 1j * b)
            elif fmt == 'MA':
                c.append(a * np.exp(1j * np.radians(b)))
            else:
                c.append(10 ** (a / 20) * np.exp(1j * np.radians(b)))
        f.append(v[0] * unit); S.append(c)
    S = np.array(S); n = int(round(math.sqrt(S.shape[1])))
    return np.array(f), S.reshape(len(f), n, n)


def sparams(src, out_prefix, title=''):
    """ src: .s2p (or .s1p) -> _s.png (dB), _smith.png, _gd.png; returns the dict of files + key numbers """
    f, S = read_touchstone(src)
    fg = f / 1e9
    n = S.shape[1]
    db = lambda z: 20 * np.log10(np.abs(z) + 1e-15)
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.plot(fg, db(S[:, 0, 0]), color='#c8531d', lw=1.4, label='|S11|')
    if n > 1:
        ax.plot(fg, db(S[:, 1, 0]), color='#1f5fa8', lw=1.4, label='|S21|')
        ax.plot(fg, db(S[:, 1, 1]), color='#e8a985', lw=0.8, ls='--', label='|S22|')
    ax.axhline(-10, color='#999', lw=0.7, ls=':'); ax.axhline(-20, color='#999', lw=0.7, ls=':')
    ax.set_xlim(fg[0], fg[-1]); ax.set_ylim(min(-50, float(db(S[:, 0, 0]).min()) - 3), max(3, float(db(S[:, 0, 0]).max()) + 2)); ax.grid(color='#ddd'); ax.legend(fontsize=8)
    if db(S[:, 0, 0]).max() > 0:
        ax.text(0.02, 0.92, '|S11| > 0 dB: unstable / non-physical run', transform=ax.transAxes, color='#c8531d', fontsize=11, weight='bold')
    ax.set_xlabel('f [GHz]'); ax.set_ylabel('dB'); ax.set_title(title or os.path.basename(src), fontsize=10)
    fig.tight_layout(); fig.savefig(out_prefix + '_s.png', dpi=110); plt.close(fig)
    # Smith
    fig, ax = plt.subplots(figsize=(5.4, 5.4))
    t = np.linspace(0, 2 * np.pi, 361)
    ax.plot(np.cos(t), np.sin(t), color='#444', lw=0.8)
    for r in (0.2, 0.5, 1, 2, 5):
        c, rr = r / (1 + r), 1 / (1 + r)
        ax.plot(c + rr * np.cos(t), rr * np.sin(t), color='#bbb', lw=0.5)
    for x in (0.5, 1, 2, 5):
        for sg in (1, -1):
            cx, cy, rr = 1, sg / x, 1 / x
            a = np.linspace(0, 2 * np.pi, 721); px, py = cx + rr * np.cos(a), cy + rr * np.sin(a)
            m = px ** 2 + py ** 2 <= 1.0001
            ax.plot(px[m], py[m], color='#bbb', lw=0.5)
    ax.plot([-1, 1], [0, 0], color='#bbb', lw=0.5)
    g = S[:, 0, 0]
    ax.scatter(g.real, g.imag, c=fg, cmap='viridis', s=6, zorder=3)
    ax.plot(g.real, g.imag, color='#c8531d', lw=0.8, alpha=0.7)
    for k in np.linspace(0, len(fg) - 1, 6).astype(int):
        ax.annotate('%.1f' % fg[k], (g[k].real, g[k].imag), fontsize=7, xytext=(3, 3), textcoords='offset points')
    ax.set_aspect('equal'); ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05); ax.axis('off')
    ax.set_title('S11 Smith (colour = f, GHz)', fontsize=10)
    fig.tight_layout(); fig.savefig(out_prefix + '_smith.png', dpi=110); plt.close(fig)
    res = dict(s=out_prefix + '_s.png', smith=out_prefix + '_smith.png', f0_GHz=float(fg[0]), f1_GHz=float(fg[-1]))
    if n > 1:
        ph = np.unwrap(np.angle(S[:, 1, 0]))
        gd = -np.gradient(ph, f) / (2 * np.pi) * 1e12
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
        ax[0].plot(fg, np.degrees(ph), color='#1f5fa8', lw=1.2); ax[0].set_title('S21 phase (unwrapped) [°]', fontsize=9); ax[0].grid(color='#ddd')
        ax[1].plot(fg, gd, color='#1f5fa8', lw=1.2); ax[1].set_title('group delay [ps]', fontsize=9); ax[1].grid(color='#ddd')
        lo = np.percentile(gd, 5); hi = np.percentile(gd, 95); ax[1].set_ylim(lo - 0.3 * (hi - lo) - 5, hi + 0.3 * (hi - lo) + 5)
        for a in ax:
            a.set_xlabel('f [GHz]'); a.set_xlim(fg[0], fg[-1]); a.tick_params(labelsize=8)
        fig.tight_layout(); fig.savefig(out_prefix + '_gd.png', dpi=110); plt.close(fig)
        res['gd'] = out_prefix + '_gd.png'
        res['s21_min_dB'] = float(db(S[:, 1, 0]).min())
    res['s11_max_dB'] = float(db(S[:, 0, 0]).max())
    return res


# ---------------------------------------------------------------- engine progress
_PROG = re.compile(r'\[@\s*(?:(\d+)m)?\s*(\d+)s\]\s+Timestep:\s+(\d+)\s+\|\|\s+Speed:\s+([\d.]+)\s+MC/s.*?Energy:\s+~([\d.e+-]+)\s+\(\s*(-?)\s*([\d.]+)dB\)')


def progress(stdout_txt):
    ts, t, sp, en, edb = [], [], [], [], []
    done, maxts, dt = False, None, None
    try:
        for line in open(stdout_txt, errors='replace'):
            m = _PROG.search(line)
            if m:
                ts.append(int(m.group(3))); t.append(int(m.group(1) or 0) * 60 + int(m.group(2)))
                sp.append(float(m.group(4))); en.append(float(m.group(5))); edb.append(-float(m.group(7)) if m.group(6) else float(m.group(7)))
            elif 'Max. number of timesteps:' in line:
                maxts = int(re.search(r'timesteps:\s*(\d+)', line).group(1))
            elif 'FDTD timestep is' in line:
                dt = float(re.search(r'is:\s*([\d.e+-]+)', line).group(1))
            elif line.startswith('Time for') or 'Time for' in line:
                done = True
    except FileNotFoundError:
        pass
    div = False
    if len(en) > 6:
        a = np.array(en)
        mn = a[len(a) // 3:].min()
        div = bool(a[-1] > 100 * mn and a[-1] > 3 * a[-min(8, len(a) - 1)])     # energy rose 20 dB above its floor and still rising
    return dict(ts=ts, t_s=t, speed=sp, energy=en, energy_dB=edb, done=done, max_ts=maxts, dt_s=dt, diverging=div)


def progress_png(stdout_txt, out):
    p = progress(stdout_txt)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    if p['ts']:
        en = np.array(p['energy']); ref = en[:max(3, len(en) // 10)].max() or 1.0
        edb = 10 * np.log10(np.maximum(en, 1e-300) / ref)                # not capped at 0 dB: growth shows as positive dB
        ax[0].plot(p['ts'], edb, color='#1f5fa8', lw=1.2); ax[0].set_title('field energy [dB rel. early max] — must fall', fontsize=9)
        ax[0].axhline(-40, color='#c8531d', lw=0.7, ls='--'); ax[0].set_ylim(min(-60, float(edb.min()) - 5), max(5, float(edb.max()) + 5))
        ax[1].plot(p['ts'], p['speed'], color='#1f5fa8', lw=1.2); ax[1].set_title('speed [M cells/s]', fontsize=9); ax[1].set_ylim(0, max(p['speed']) * 1.2)
        if p['max_ts']:
            ax[0].axvline(p['max_ts'], color='#999', lw=0.7, ls=':')
        if p['diverging']:
            ax[0].text(0.02, 0.9, 'ENERGY GROWING - unstable run', transform=ax[0].transAxes, color='#c8531d', fontsize=11, weight='bold')
    for a in ax:
        a.set_xlabel('timestep'); a.grid(color='#ddd'); a.tick_params(labelsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)
    return p


if __name__ == '__main__':
    a = sys.argv[1:]
    if a and a[0].endswith('.xml'):
        print(geometry(a[0], a[1] if len(a) > 1 else a[0][:-4]))
    elif a and a[0].endswith(('.s2p', '.s1p')):
        print(sparams(a[0], a[1] if len(a) > 1 else a[0][:-4]))
    elif a and a[0].endswith('.txt'):
        print(progress_png(a[0], a[1] if len(a) > 1 else a[0][:-4] + '_progress.png'))
    else:
        print(__doc__)
