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
kicad2view: convert a KiCad board into the data file of the z2labs board viewer.

    python kicad2view.py board.kicad_pcb [board.z2v.json] [--fab JLCPCB:JLC04161H-7628]

The board is read with CSXCAD.KiCadImport (no KiCad needed). Per copper layer the
shapes of one net are merged, the drill holes are cut out and the result is
triangulated, so the viewer only uploads triangles. Dielectric and mask layers
are the board outline minus the holes. Components are boxes over their fab /
courtyard outline; the height comes from the footprint's 3D model (VRML or STEP
bounding box with the model offset, scale and rotation) and, where no model file is
found, is estimated from the footprint name and marked as estimated.
"""

import base64
import json
import math
import os
import re
import sys
import time

import numpy as np
import shapely
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union

FORMAT_VERSION = 1


def _load_kicadimport():
    try:
        from CSXCAD import KiCadImport, FabStackups
        return KiCadImport, FabStackups
    except ImportError:
        pass
    # Fallback: the pure-Python modules of the CSXCAD submodule, without the compiled package
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    base = os.environ.get('Z2_CSXCAD_PY', os.path.join(here, '..', '..', 'CSXCAD', 'python', 'CSXCAD'))
    mods = []
    for name in ('KiCadImport', 'FabStackups'):
        spec = importlib.util.spec_from_file_location('CSXCAD.' + name, os.path.join(base, name + '.py'))
        m = importlib.util.module_from_spec(spec)
        sys.modules['CSXCAD.' + name] = m
        spec.loader.exec_module(m)
        mods.append(m)
    return mods


def b64(a, dtype):
    return base64.b64encode(np.ascontiguousarray(a, dtype=dtype).tobytes()).decode('ascii')


# ---------------------------------------------------------------------------- geometry

def _poly(shape):
    p = Polygon(shape.pts, [h for h in shape.holes if len(h) >= 3]) if len(shape.pts) >= 3 else None
    if p is not None and not p.is_valid:
        p = shapely.make_valid(p)
    return p


def _polygons(g):
    """ All polygons of a geometry (drops lines and points left by boolean operations). """
    if g is None or g.is_empty:
        return []
    if isinstance(g, Polygon):
        return [g]
    if isinstance(g, MultiPolygon):
        return list(g.geoms)
    out = []
    for p in getattr(g, 'geoms', []):
        out += _polygons(p)
    return out


def triangulate(polys):
    """ Triangles (T,3,2) of a list of polygons (constrained Delaunay, holes respected). """
    if not polys:
        return np.zeros((0, 3, 2))
    if hasattr(shapely, 'constrained_delaunay_triangles'):          # shapely >= 2.1
        parts = shapely.get_parts(shapely.constrained_delaunay_triangles(np.array(polys, dtype=object)))
    else:
        # older shapely (e.g. in Pyodide): Delaunay of all vertices, keep the triangles inside
        parts = []
        for p in polys:
            t = shapely.get_parts(shapely.delaunay_triangles(p))
            if len(t):
                parts.append(t[shapely.contains(p, shapely.centroid(t))])
        parts = np.concatenate(parts) if parts else []
    if len(parts) == 0:
        return np.zeros((0, 3, 2))
    c = shapely.get_coordinates(shapely.get_exterior_ring(parts))      # 4 points per triangle, closed
    return c.reshape(-1, 4, 2)[:, :3, :]


def triangulate_each(polys):
    """ triangulate() per polygon; returns triangles and the polygon index of every triangle. """
    out, idx = [], []
    for i, p in enumerate(polys):
        t = triangulate([p])
        out.append(t)
        idx.append(np.full(len(t), i))
    if not out:
        return np.zeros((0, 3, 2)), np.zeros(0, int)
    return np.concatenate(out), np.concatenate(idx)


def rings(polys):
    """ Boundary rings (outer and holes) as a list of (N,2) arrays without the closing point. """
    out = []
    for p in polys:
        for r in [p.exterior] + list(p.interiors):
            c = np.asarray(r.coords)[:-1]
            if len(c) >= 3:
                out.append(c)
    return out


# ---------------------------------------------------------------------------- component heights

# Body height in mm by package name; first match wins. Values are typical datasheet heights;
# they only make the 3D view plausible and are flagged as estimated.
_HEIGHTS = [
    (r'ESP32[-_]?(WROOM|WROVER|MINI)|WROOM|WROVER', 3.2),
    (r'RJ45|RJP|LPJ\d|8P8C|MAGJACK', 13.5),
    (r'USB[-_ ]?A\b|USB_A', 7.0),
    (r'USB', 3.2),
    (r'HDMI', 6.0),
    (r'SD[-_ ]?CARD|MICRO[-_ ]?SD|TF[-_ ]?CARD', 1.9),
    (r'PIN[-_ ]?HEADER|PINHEADER|GBH|HDR|HEADER|2\.54|1X\d|2X\d', 8.5),
    (r'TERMINAL|SCREW|DG\d|KF\d', 10.0),
    (r'CP_ELEC|ELEC|ELECTROLYT|_CAP_\d+X\d+', 6.0),
    (r'HC[-_]?49', 3.6),
    (r'BUTTON|SWITCH|TACT|IT11\d\d|SW_|_SW', 3.5),
    (r'BATTERY|CR20\d\d|CR12\d\d|BAT', 4.0),
    (r'TO[-_]?220', 15.0), (r'TO[-_]?263|D2PAK', 4.5), (r'TO[-_]?252|DPAK', 2.4),
    (r'SOT[-_]?223', 1.7), (r'SOT[-_]?89', 1.5), (r'SOT[-_]?23|SOT[-_]?323|SC[-_]?70|SOT[-_]?353|SOT[-_]?363', 1.1),
    (r'SOD[-_]?123', 1.1), (r'SOD[-_]?323|SOD[-_]?523', 0.9),
    (r'\bSM[ABC]\b|DO[-_]?214|_SM[ABC]\b|\(SM[ABC]\)', 2.3),
    (r'LQFP|TQFP|QFP', 1.5), (r'BGA', 1.3), (r'QFN|DFN|SON|LGA', 0.9),
    (r'TSSOP|MSOP|SSOP|VSSOP', 1.1), (r'SOIC|\bSO[-_]?\d|SOP', 1.6),
    (r'CRYSTAL|XTAL|OSC', 1.0),
    (r'0201', 0.25), (r'0402', 0.35), (r'0603', 0.5), (r'0805', 0.6), (r'1206', 0.65),
    (r'1210', 0.8), (r'1812|2010|2512', 0.7),
]
_NO_BODY = re.compile(r'LOGO|FIDUCIAL|MOUNTINGHOLE|MOUNTING_HOLE|TESTPOINT|TEST_POINT|\bTP\d|JUMPER_SMD|SOLDERJUMPER|'
                      r'SOLDER_JUMPER|\bSJ_|\bSJ\d|NETTIE|NET_TIE|KICAD|OSHW|CE[-_]LOGO|ARROW|HOLE', re.I)


def estimate_height(lib_id, value, area):
    """ Body height in mm: explicit LxWxH in the footprint name, then name patterns, then the value. """
    lib = lib_id.split(':')[-1].upper()
    m = re.search(r'(\d+(?:\.\d+)?)X(\d+(?:\.\d+)?)X(\d+(?:\.\d+)?)(?:MM)?', lib)
    if m:
        return float(m.group(3))
    m = re.search(r'(?:CPOL|CP_|ELEC).*?(\d+(?:\.\d+)?)X(\d+(?:\.\d+)?)', lib)
    if m:                   # radial electrolytic: diameter x height
        return float(m.group(2))
    for text in (lib, (value or '').upper()):
        for pat, h in _HEIGHTS:
            if re.search(pat, text):
                return h
    return float(np.clip(0.6 + 0.12 * math.sqrt(max(area, 0.0)), 0.6, 3.0))


# ---------------------------------------------------------------------------- 3D model heights

_MODEL_DIRS = [r'C:\Program Files\KiCad\{v}\share\kicad\3dmodels', '/usr/share/kicad/3dmodels',
               '/Applications/KiCad/KiCad.app/Contents/SharedSupport/3dmodels']
_model_cache = {}


def _resolve_model(path, board_dir):
    def var(m):
        name = m.group(1)
        if name in os.environ:
            return os.environ[name]
        if name == 'KIPRJMOD':
            return board_dir
        mv = re.match(r'KICAD(\d+)_3DMODEL_DIR', name)
        if mv:
            for d in _MODEL_DIRS:
                d = d.format(v='%s.0' % mv.group(1))
                if os.path.isdir(d):
                    return d
        return m.group(0)
    p = re.sub(r'\$\{([^}]+)\}', var, path)
    if not os.path.isabs(p):
        p = os.path.join(board_dir, p)
    stem = os.path.splitext(p)[0]
    for cand in (p, stem + '.wrl', stem + '.step', stem + '.stp', stem + '.WRL', stem + '.STEP'):
        if os.path.isfile(cand):
            return cand
    return None


def _model_bbox(path):
    """ Bounding box (mm) of a VRML (KiCad units: 0.1 inch) or STEP model, cached. """
    if path in _model_cache:
        return _model_cache[path]
    bb = None
    try:
        txt = open(path, 'r', errors='ignore').read()
        if path.lower().endswith('.wrl'):
            pts = []
            for blk in re.findall(r'point\s*\[([^\]]*)\]', txt):
                v = np.array(re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', blk), float)
                pts.append(v[: len(v) // 3 * 3].reshape(-1, 3))
            if pts:
                a = np.vstack(pts) * 2.54
                bb = (a.min(0), a.max(0))
        else:
            v = re.findall(r'CARTESIAN_POINT\s*\(\s*\'[^\']*\'\s*,\s*\(([^)]*)\)', txt)
            a = np.array([[float(c) for c in t.split(',')[:3]] for t in v if t.count(',') >= 2])
            if len(a):
                bb = (a.min(0), a.max(0))
    except (OSError, ValueError):
        bb = None
    _model_cache[path] = bb
    return bb


# ---------------------------------------------------------------------------- 3D model meshes (VRML)

_VRML_TOK = re.compile(r'"[^"]*"|[{}\[\]]|[^\s{}\[\],"]+')
_VRML_NUM = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$')
_mesh_cache = {}


class _Vrml:
    """ Minimal VRML 2.0 reader: the node tree with DEF/USE, enough for KiCad / VTK / FreeCAD model files. """

    def __init__(self, txt):
        txt = re.sub(r'#[^\n]*', '', txt)
        self.t = _VRML_TOK.findall(txt)
        self.defs = {}
        self.nodes = []
        i, n = 0, len(self.t)
        while i < n:
            if self.t[i] == 'ROUTE':
                i += 4
            elif self._is_node(i):
                node, i = self._node(i)
                self.nodes.append(node)
            else:
                i += 1

    def _is_node(self, i):
        t = self.t
        return i < len(t) and (t[i] in ('DEF', 'USE') or (i + 1 < len(t) and t[i + 1] == '{'))

    def _node(self, i):
        t = self.t
        if t[i] == 'USE':
            return self.defs.get(t[i + 1], {}), i + 2
        name = None
        if t[i] == 'DEF':
            name, i = t[i + 1], i + 2
        node = {'_type': t[i]}
        i += 2
        while i < len(t) and t[i] != '}':
            f = t[i]
            v, i = self._value(i + 1)
            node[f] = v
        if name:
            self.defs[name] = node
        return node, i + 1

    def _value(self, i):
        t = self.t
        if t[i] == '[':
            j = i + 1
            if t[j] == ']':
                return [], j + 1
            if self._is_node(j):
                items = []
                while t[j] != ']':
                    nd, j = self._node(j)
                    items.append(nd)
                return items, j + 1
            k = t.index(']', j)
            return t[j:k], k + 1
        if self._is_node(i):
            return self._node(i)
        j = i
        while j < len(t) and _VRML_NUM.match(t[j]):
            j += 1
        return (t[i], i + 1) if j == i else (t[i:j], j)


def _nums(v, default=()):
    if v is None or isinstance(v, (dict, str)) or not len(v):
        return np.array(default, float)
    return np.array(v, float)


def _axis_angle(r):
    x, y, z, a = r
    n = math.sqrt(x * x + y * y + z * z) or 1.0
    x, y, z = x / n, y / n, z / n
    c, s, C = math.cos(a), math.sin(a), 1 - math.cos(a)
    return np.array([[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
                     [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
                     [z * x * C - y * s, z * y * C + x * s, c + z * z * C]])


def _faces_to_tris(idx):
    """ coordIndex (faces separated by -1) to triangles (fan). """
    if len(idx) % 4 == 0 and (idx[3::4] == -1).all():
        return idx.reshape(-1, 4)[:, :3]
    tri, face = [], []
    for v in np.append(idx, -1).tolist():
        if v < 0:
            tri += [(face[0], face[k], face[k + 1]) for k in range(1, len(face) - 1)]
            face = []
        else:
            face.append(v)
    return np.array(tri, np.int64).reshape(-1, 3)


def _vrml_walk(node, M, out):
    typ = node.get('_type', '')
    if typ == 'Transform':
        T = np.eye(4); T[:3, 3] = _nums(node.get('translation'), (0, 0, 0))[:3]
        Cn = np.eye(4); Cn[:3, 3] = _nums(node.get('center'), (0, 0, 0))[:3]
        R = np.eye(4); R[:3, :3] = _axis_angle(_nums(node.get('rotation'), (0, 0, 1, 0))[:4])
        S = np.eye(4); S[:3, :3] = np.diag(_nums(node.get('scale'), (1, 1, 1))[:3])
        SR = np.eye(4); SR[:3, :3] = _axis_angle(_nums(node.get('scaleOrientation'), (0, 0, 1, 0))[:4])
        M = M @ T @ Cn @ R @ SR @ S @ SR.T @ np.linalg.inv(Cn)
    if typ == 'Shape':
        g = node.get('geometry')
        if not (isinstance(g, dict) and g.get('_type') == 'IndexedFaceSet'):
            return
        coord = g.get('coord')
        pts = _nums(coord.get('point') if isinstance(coord, dict) else None)
        idx = _nums(g.get('coordIndex')).astype(np.int64)
        if len(pts) < 9 or not len(idx):
            return
        pts = pts[: len(pts) // 3 * 3].reshape(-1, 3)
        tri = _faces_to_tris(idx)
        tri = tri[(tri < len(pts)).all(1) & (tri >= 0).all(1)]
        col = (0.8, 0.8, 0.8)
        app = node.get('appearance')
        mat = app.get('material') if isinstance(app, dict) else None
        if isinstance(mat, dict):
            col = tuple(_nums(mat.get('diffuseColor'), col)[:3])
            if _nums(mat.get('transparency'), (0,))[0] > 0.9:
                return
        if len(tri):
            out.append((pts @ M[:3, :3].T + M[:3, 3], tri, col))
        return
    if typ == 'Switch':
        v = node.get('choice') or node.get('children') or []
        w = int(_nums(node.get('whichChoice'), (-1,))[0])
        kids = [v[w]] if isinstance(v, list) and 0 <= w < len(v) else []
    elif typ == 'LOD':
        kids = (node.get('level') or node.get('children') or [])[:1]
    else:
        kids = node.get('children') or []
    if isinstance(kids, dict):
        kids = [kids]
    for c in kids:
        if isinstance(c, dict):
            _vrml_walk(c, M, out)


def _model_mesh(path, max_tri=200000):
    """ Triangles of a VRML model in mm (KiCad VRML unit = 2.54 mm): (pos (N,3), tri (T,3), col (N,3) u8), cached.
        STEP files and meshes above max_tri give None (the component is then drawn as a box). """
    if path in _mesh_cache:
        return _mesh_cache[path]
    res = None
    if path.lower().endswith('.wrl'):
        try:
            doc = _Vrml(open(path, 'r', errors='ignore').read())
            parts = []
            for nd in doc.nodes:
                _vrml_walk(nd, np.eye(4), parts)
            if parts:
                pos, tri, col, off = [], [], [], 0
                for p, t, c in parts:
                    pos.append(p)
                    tri.append(t + off)
                    col.append(np.tile(np.round(np.clip(c, 0, 1) * 255), (len(p), 1)))
                    off += len(p)
                tri = np.concatenate(tri)
                if 0 < len(tri) <= max_tri:
                    pos = np.concatenate(pos) * 2.54
                    used = np.unique(tri)               # drop vertices no triangle uses
                    remap = np.full(len(pos), -1, np.int64)
                    remap[used] = np.arange(len(used))
                    res = (pos[used], remap[tri], np.concatenate(col)[used].astype(np.uint8))
        except (OSError, ValueError, IndexError, KeyError) as e:
            print('  model %s: %s' % (os.path.basename(path), e))
    _mesh_cache[path] = res
    return res


def _model_matrix(F, m):
    """ 4x4 from model coordinates (mm) to board coordinates (mm, y up, z from the board surface of F's side). """
    M = np.eye(4)
    M[:3, :3] = _rot(np.eye(3), m['rotate']).T @ np.diag(m['scale'])
    M[:3, 3] = m['offset']
    if F.side == 'B':                       # KiCad flips back-side models by 180 deg about the Y axis
        M = np.diag([-1.0, 1.0, -1.0, 1.0]) @ M
    a = math.radians(F.angle or 0.0)
    P = np.eye(4)
    P[:2, :2] = [[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]]
    P[:2, 3] = (F.x, F.y)
    return P @ M


def footprint_meshes(F, board_dir, store):
    """ Mesh instances of F's VRML models: list of {model key, 4x4 column-major}; meshes go into store once each. """
    inst = []
    for m in F.models:
        p = _resolve_model(m['path'], board_dir)
        if p and not p.lower().endswith('.wrl'):
            stem = os.path.splitext(p)[0]
            p = next((c for c in (stem + '.wrl', stem + '.WRL') if os.path.isfile(c)), p)
        mesh = _model_mesh(p) if p else None
        if mesh is None:
            continue
        key = os.path.basename(p)
        if key not in store:
            pos, tri, col = mesh
            store[key] = dict(pos=b64(pos.reshape(-1), '<f4'), idx=b64(tri.reshape(-1), '<u4'), col=b64(col.reshape(-1), '<u1'),
                              n_tri=int(len(tri)))
        inst.append(dict(model=key, m=[round(float(v), 6) for v in _model_matrix(F, m).T.reshape(-1)]))
    return inst


# ---------------------------------------------------------------------------- simulation bodies (CSX)

_BOARD_PROP = re.compile(r'^(diel_|plated_holes|mask)|_Cu$|^port_|^Dump|^probe', re.I)


def _box_mesh(p1, p2):
    (x0, y0, z0), (x1, y1, z1) = np.minimum(p1, p2), np.maximum(p1, p2)
    v = np.array([[x, y, z] for z in (z0, z1) for y in (y0, y1) for x in (x0, x1)], float)
    f = [(0, 2, 1), (1, 2, 3), (4, 5, 6), (5, 7, 6), (0, 1, 4), (1, 5, 4), (2, 6, 3), (3, 6, 7), (0, 4, 2), (2, 4, 6), (1, 3, 5), (3, 7, 5)]
    return v, np.array(f)


def _tube_mesh(p1, p2, r_out, r_in=0.0, n=28):
    p1, p2 = np.asarray(p1, float), np.asarray(p2, float)
    ax = p2 - p1
    L = np.linalg.norm(ax)
    if L == 0:
        return np.zeros((0, 3)), np.zeros((0, 3), int)
    ax = ax / L
    u = np.cross(ax, [0, 0, 1] if abs(ax[2]) < 0.9 else [1, 0, 0])
    u /= np.linalg.norm(u)
    w = np.cross(ax, u)
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.cos(a)[:, None] * u + np.sin(a)[:, None] * w
    radii = [r_out] + ([r_in] if r_in > 0 else [])
    V, F = [], []
    for k, r in enumerate(radii):                       # side walls
        b = len(V)
        V += list(p1 + r * ring) + list(p2 + r * ring)
        for i in range(n):
            j = (i + 1) % n
            F += [(b + i, b + j, b + n + i), (b + j, b + n + j, b + n + i)]
    if r_in > 0:                                        # annular end caps
        for e in (0, n):
            for i in range(n):
                j = (i + 1) % n
                F += [(e + i, 2 * n + e + i, e + j), (e + j, 2 * n + e + i, 2 * n + e + j)]
    else:                                               # disc end caps
        for c, e in ((p1, 0), (p2, n)):
            ci = len(V)
            V.append(c)
            F += [(ci, e + i, e + (i + 1) % n) for i in range(n)]
    return np.array(V), np.array(F)


def sim_bodies(xml_path, unit_mm=None):
    """ Non-board bodies of a CSX model file (connector pins, bodies, coax feeds, ...) as triangle meshes in
        board coordinates (mm, the importer's frame): what the simulation really meshes, for the viewer. """
    import xml.etree.ElementTree as ET
    root = ET.parse(xml_path).getroot()
    cs = root if root.tag == 'ContinuousStructure' else root.find('.//ContinuousStructure')
    grid = cs.find('RectilinearGrid')
    scale = (float(grid.get('DeltaUnit', '0.001')) * 1e3) if unit_mm is None and grid is not None else (unit_mm or 1.0)
    pt = lambda e: np.array([float(e.get(k)) for k in 'XYZ']) * scale
    out = []
    for prop in cs.find('Properties'):
        name = prop.get('Name', '')
        if prop.tag not in ('Metal', 'Material', 'ConductingSheet') or _BOARD_PROP.search(name):
            continue
        V, F = [], []
        prims = prop.find('Primitives')
        for q in (prims if prims is not None else []):
            if q.tag == 'Box':
                v, f = _box_mesh(pt(q.find('P1')), pt(q.find('P2')))
            elif q.tag == 'Cylinder':
                v, f = _tube_mesh(pt(q.find('P1')), pt(q.find('P2')), float(q.get('Radius')) * scale)
            elif q.tag == 'CylindricalShell':
                r, sw = float(q.get('Radius')) * scale, float(q.get('ShellWidth')) * scale
                v, f = _tube_mesh(pt(q.find('P1')), pt(q.find('P2')), r + sw / 2, max(r - sw / 2, 1e-4))
            else:
                continue
            F.append(f + sum(len(x) for x in V))
            V.append(v)
        if V:
            V, F = np.concatenate(V), np.concatenate(F)
            kind = 'diel' if prop.tag == 'Material' else 'metal'
            er = None
            if kind == 'diel':
                pr = prop.find('Property')
                er = float(pr.get('Epsilon', '1').split(',')[0]) if pr is not None else None   # anisotropic: x,y,z
            out.append(dict(name=name, kind=kind, er=er, pos=b64(V.reshape(-1), '<f4'), idx=b64(F.reshape(-1), '<u4'), n_tri=int(len(F))))
    return out


def _rot(v, deg):
    rx, ry, rz = (math.radians(-d) for d in deg)      # KiCad rotates the model clockwise for positive angles
    cx, sx, cy, sy, cz, sz = math.cos(rx), math.sin(rx), math.cos(ry), math.sin(ry), math.cos(rz), math.sin(rz)
    R = (np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]]) @ np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
         @ np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]]))
    return v @ R.T


def model_height(F, board_dir):
    """ Body height above the board in mm from the footprint's 3D models, or None. """
    best = None
    for m in F.models:
        p = _resolve_model(m['path'], board_dir)
        bb = _model_bbox(p) if p else None
        if bb is None:
            continue
        lo, hi = bb
        corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])
        c = _rot(corners * np.array(m['scale']), m['rotate']) + np.array(m['offset'])
        h = float(c[:, 2].max())
        if h < 0.05:          # rotation sign conventions differ between model sources: take the extent
            h = float(-c[:, 2].min())
        if 0.05 < h < 60:
            best = max(best or 0.0, h)
    return best


# ---------------------------------------------------------------------------- export

def export(board, name='', simplify=0.002, log=print, board_dir='.', with_models=True, sim_xml=None):
    t0 = time.time()
    cu = board.copper_layers
    cu_index = {l: i for i, l in enumerate(cu)}
    nets = sorted({s.net for s in board.shapes} | {d.net for d in board.drills})
    if '' in nets:
        nets.remove('')
    nets = [''] + nets
    net_index = {n: i for i, n in enumerate(nets)}

    outline = unary_union([Polygon(o) for o in board.outline if len(o) >= 3]) if board.outline else None
    if outline is None or outline.is_empty:
        (x0, x1), (y0, y1) = board.BoundingBox()[:2]
        outline = Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    if len(board.outline) > 1:
        # several closed outlines: the largest is the board, the others are cut-outs inside it
        ps = sorted((Polygon(o) for o in board.outline if len(o) >= 3), key=lambda p: -p.area)
        outline = ps[0]
        for p in ps[1:]:
            if ps[0].contains(p):
                outline = outline.difference(p)
            else:
                outline = outline.union(p)

    # drill holes, per copper layer (a hole cuts only the layers it spans; NPTH cut everything)
    def hole_geom(d):
        if d.slot:
            L, a = d.slot
            dx, dy = 0.5 * max(L - d.d, 0) * math.cos(math.radians(a)), 0.5 * max(L - d.d, 0) * math.sin(math.radians(a))
            return shapely.LineString([(d.x - dx, d.y - dy), (d.x + dx, d.y + dy)]).buffer(d.d / 2, 8) if (dx or dy) else \
                shapely.Point(d.x, d.y).buffer(d.d / 2, 8)
        return shapely.Point(d.x, d.y).buffer(d.d / 2, 8)

    holes_all = []
    holes_per = {l: [] for l in cu}
    for d in board.drills:
        g = hole_geom(d)
        holes_all.append(g)
        if not d.plated:
            for l in cu:
                holes_per[l].append(g)
            continue
        a, b = cu_index.get(d.layers[0], 0), cu_index.get(d.layers[1], len(cu) - 1)
        for l in cu[min(a, b):max(a, b) + 1]:
            holes_per[l].append(g)
    holes_union = unary_union(holes_all) if holes_all else None
    log('  holes: %d (%.1f s)' % (len(holes_all), time.time() - t0))

    layers = {}
    for l in cu:
        groups = {}
        for s in board.ShapesOn(l):
            p = _poly(s)
            if p is not None and not p.is_empty:
                groups.setdefault(s.net, []).append(p)
        hu = unary_union(holes_per[l]) if holes_per[l] else None
        polys, poly_net = [], []
        for net, ps in groups.items():
            g = unary_union(ps)
            if hu is not None:
                g = g.difference(hu)
            g = g.intersection(outline.buffer(0.01))
            if simplify:
                g = g.simplify(simplify, preserve_topology=True)
            for p in _polygons(g):
                if p.area > 1e-6:
                    polys.append(p)
                    poly_net.append(net_index[net])
        tris, tpi = triangulate_each(polys)
        rr, rnet = [], []
        for p, n in zip(polys, poly_net):
            for r in rings([p]):
                rr.append(r)
                rnet.append(n)
        pnet = np.asarray(poly_net, int)
        layers[l] = dict(tris=b64(tris.reshape(-1), '<f4'), tri_net=b64(pnet[tpi] if len(tpi) else [], '<u4'),
                         rings=b64(np.concatenate(rr).reshape(-1) if rr else [], '<f4'),
                         ring_len=b64([len(r) for r in rr], '<u4'), ring_net=b64(rnet, '<u4'),
                         n_tri=int(len(tris)), area=round(sum(p.area for p in polys), 3))
        log('  %-8s %6d polygons %8d triangles (%.1f s)' % (l, len(polys), len(tris), time.time() - t0))

    # dielectric / mask: outline minus all holes
    body = outline.difference(holes_union) if holes_union is not None else outline
    bpolys = _polygons(body.simplify(simplify, preserve_topology=True) if simplify else body)
    btris = triangulate(bpolys)
    brings = rings(bpolys)
    dielectric = dict(tris=b64(btris.reshape(-1), '<f4'), rings=b64(np.concatenate(brings).reshape(-1) if brings else [], '<f4'),
                      ring_len=b64([len(r) for r in brings], '<u4'), n_tri=int(len(btris)))

    # plated barrels (vias and plated pad holes)
    pl = [d for d in board.drills if d.plated]
    drills = dict(n=len(pl), x=b64([d.x for d in pl], '<f4'), y=b64([d.y for d in pl], '<f4'), d=b64([d.d for d in pl], '<f4'),
                  top=b64([min(cu_index.get(d.layers[0], 0), cu_index.get(d.layers[1], len(cu) - 1)) for d in pl], '<u1'),
                  bot=b64([max(cu_index.get(d.layers[0], 0), cu_index.get(d.layers[1], len(cu) - 1)) for d in pl], '<u1'),
                  net=b64([net_index.get(d.net, 0) for d in pl], '<u4'), via=b64([d.kind == 'via' for d in pl], '<u1'))

    # components
    pads_by_ref = {}
    for s in board.shapes:
        if s.ref and s.kind == 'pad':
            pads_by_ref.setdefault(s.ref, []).append(s.pts)
    comps, model_store = [], {}
    for F in board.footprints.values():
        if _NO_BODY.search(F.lib_id) or _NO_BODY.search(F.ref.split('#')[0] if F.ref.startswith('LOGO') else ''):
            continue
        src, pts = None, None
        for k in ('fab', 'courtyard'):
            v = getattr(F, k, None)
            if v is not None and len(v) >= 3:
                src, pts = k, v
                break
        if pts is None:
            v = pads_by_ref.get(F.ref)
            if not v:
                continue
            src, pts = 'pads', np.concatenate(v)
        hull = shapely.MultiPoint(pts).convex_hull
        if src == 'courtyard':
            hull = hull.buffer(-0.15, join_style='mitre')
        if not isinstance(hull, Polygon) or hull.area < 0.02:
            continue
        hull = hull.simplify(0.02)
        h = model_height(F, board_dir)
        est = h is None
        if est:
            h = estimate_height(F.lib_id, F.value, hull.area)
        c = dict(ref=F.ref, value=F.value, lib=F.lib_id, side=F.side, x=round(F.x, 4), y=round(F.y, 4),
                 angle=F.angle, h=round(h, 3), h_est=est, src=src, dnp=bool(F.dnp), models=len(F.models),
                 hull=[[round(a, 4), round(b, 4)] for a, b in np.asarray(hull.exterior.coords)[:-1]])
        if with_models:
            inst = footprint_meshes(F, board_dir, model_store)
            if inst:
                c['meshes'] = inst
        comps.append(c)

    stack = [dict(name=s.name, type=s.type, t=s.thickness, er=s.epsilon_r, tand=s.loss_tangent, material=s.material,
                  z0=s.z0, z1=s.z1) for s in board.stackup]
    (x0, x1), (y0, y1) = outline.bounds[0::2], outline.bounds[1::2]
    data = dict(format='z2v', version=FORMAT_VERSION, name=name, generator='kicad2view',
                stackup_source=board.stackup_source, thickness=board.thickness, bbox=[x0, y0, x1, y1],
                copper=cu, layer_user_names=getattr(board, 'layer_user_names', {}), stack=stack, nets=nets,
                layers=layers, dielectric=dielectric, drills=drills, components=comps, models=model_store,
                warnings=board.warnings)
    data['ports'] = port_candidates(board, board_dir, name)
    if sim_xml is None and name:                        # the simulation model next to the board: <name>_sim.xml
        cand = os.path.join(board_dir, name + '_sim.xml')
        sim_xml = cand if os.path.isfile(cand) else None
    if sim_xml:
        try:
            data['sim_bodies'] = sim_bodies(sim_xml)
            data['sim_source'] = os.path.basename(sim_xml)
            log('  %d simulation bodies from %s' % (len(data['sim_bodies']), os.path.basename(sim_xml)))
        except Exception as e:                           # a broken / foreign XML must not stop the board view
            log('  simulation bodies: %s' % e)
    log('  %d components (%d heights from 3D models, %d drawn from %d VRML models, %dk triangles), %d plated holes, '
        '%d RF port candidates, total %.1f s'
        % (len(comps), sum(not c['h_est'] for c in comps), sum('meshes' in c for c in comps), len(model_store),
           sum(m['n_tri'] for m in model_store.values()) // 1000, len(pl), len(data['ports']), time.time() - t0))
    return data


def port_candidates(board, board_dir='.', name=''):
    """ RF port suggestions from twin/rf_ports (when importable): list of dicts for the viewer. """
    try:
        import rf_ports
        import kicad_sim
    except ImportError:
        return []
    pro = rf_ports.load_project(os.path.join(board_dir, name + '.kicad_pcb')) if name else None
    out = []
    for i, c in enumerate(rf_ports.candidates(board, pro)):
        p = c['pad']
        try:
            port = kicad_sim.make_port(board, p, i + 1)
            geom = dict(start=[round(float(v), 4) for v in port['start']], stop=[round(float(v), 4) for v in port['stop']],
                        axis=port.get('normal_axis'), ref_net=port.get('ref_net'))
        except Exception as e:                       # no reference copper under the pad, etc.
            geom = dict(error=str(e))
        out.append(dict(ref=c['ref'], net=c['net'], layer=c['layer'], score=c['score'], confident=c['confident'], reasons=c['reasons'],
                        lib_id=c['lib_id'], value=c['value'], edge_mm=c['edge_mm'],
                        pad=[round(float(v), 4) for v in (*p.pts.min(0), *p.pts.max(0))], **geom))
    return out


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    ap.add_argument('board')
    ap.add_argument('out', nargs='?')
    ap.add_argument('--fab', help='apply a library stackup, e.g. JLCPCB:JLC04161H-7628')
    ap.add_argument('--simplify', type=float, default=0.002, help='polygon simplification tolerance in mm')
    ap.add_argument('--no-models', action='store_true', help='draw components as boxes, not their VRML 3D models')
    ap.add_argument('--sim', help='CSX model file (openEMS XML) whose connector / extra bodies are shown as the '
                                   'simulation geometry (default: <board>_sim.xml next to the board, if present)')
    a = ap.parse_args(argv)
    KiCadImport, FabStackups = _load_kicadimport()
    t = time.time()
    board = KiCadImport.KiCadBoard(a.board)
    print('parsed %s in %.1f s: %d copper layers, stackup %s' % (a.board, time.time() - t, len(board.copper_layers), board.stackup_source))
    if a.fab:
        fab, nm = a.fab.split(':', 1)
        st = FabStackups.Find(fab=fab, name=nm)
        if not st:
            sys.exit('stackup %s not found' % a.fab)
        FabStackups.Apply(board, st[0])
    name = os.path.splitext(os.path.basename(a.board))[0]
    data = export(board, name, a.simplify, board_dir=os.path.dirname(os.path.abspath(a.board)), with_models=not a.no_models, sim_xml=a.sim)
    out = a.out or os.path.splitext(a.board)[0] + '.z2v.json'
    with open(out, 'w') as fh:
        json.dump(data, fh, separators=(',', ':'))
    print('wrote %s (%.1f MB)' % (out, os.path.getsize(out) / 1e6))


if __name__ == '__main__':
    main(sys.argv[1:])
