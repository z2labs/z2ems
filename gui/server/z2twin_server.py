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
z2twin_server: the local back end of the board viewer.

Serves the viewer and a small JSON API, so the browser can open the KiCad boards of a
folder, place ports, start an openEMS run and show its result. Only files below the
served folder are accessible.

    python z2twin_server.py BOARD_FOLDER [--port 8765] [--threads 4] [--open]
    python z2twin_server.py BOARD_FOLDER --lan          other devices of the network (laptop, phone)

By default only this computer can connect (127.0.0.1). --lan (or --host ADDRESS) listens on
the network and requires an access token: the printed URL carries it once (?token=...), the
browser then keeps it in a cookie. Over Tailscale, keep the default 127.0.0.1 and publish
it with `tailscale serve --bg 8765`: HTTPS and the tailnet's own login, no token needed.

API (all JSON):
    GET  /api/health                      engine, versions, folder
    GET  /api/boards                      .kicad_pcb files below the folder
    POST /api/load    {path, fab?}        viewer data (.z2v), RF port suggestions, stackup suggestions, passives
    POST /api/snap    {path, x, y, fab?}  port at the outer-layer pad nearest to (x, y) (board mm, y up)
    POST /api/field_overlay {result, exc, k, quantity J|V, layer, phase?, lines?, image?} -> heat map + current / iso-voltage lines
    POST /api/run     {path, ports, fmax, fmin?, nf?, fab?, materials?, threads?, coarse?, engine?, converge?}   -> {job}
                      converge: also run the board with 0.6x cells and report whether the mesh is converged
                      engine: auto (GPU build when present) | cpu | cuda
    GET  /api/job/ID                      state, progress, log tail, result when done
    POST /api/job/ID/cancel
    GET  /api/job/ID/touchstone           the .sNp file
    GET  /api/runs?path=                  earlier results of a board
    POST /api/upload?name=X.kicad_pcb     raw board file from the browser -> {path} (stored in <folder>/uploads/<stem>/)
         &project=STEM                    schematic / project / symbol files of that board go next to it
    GET  /api/schematic?path=             schematic of the board's project, flattened per sheet instance, nets
                                          resolved through the board's pads (kisch.py)
    GET  /api/coherence?path=             channel phase-match spec (<board>.coherence.json) and its result
    POST /api/coherence {path, spec}      store the spec and compute it (coherence.py)
    POST /api/pad     {path, x, y}        pad under / nearest to (x, y) (board mm, y up): ref, pad, net
"""

import argparse
import gzip
import hmac
import json
import secrets
import socket
import subprocess
import mimetypes
import os
import re
import sys
import threading
import time
import traceback
import uuid

import numpy as np
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

HERE = os.path.dirname(os.path.abspath(__file__))
GUI = os.path.dirname(HERE)
REPO = os.path.dirname(GUI)
VIEWER = os.path.join(GUI, 'viewer')
for p in (os.path.join(REPO, 'twin'), os.path.join(GUI, 'tools')):
    if p not in sys.path:
        sys.path.insert(0, p)

import kicad_sim as ks          # noqa: E402
ks._setup_python()
import rf_ports                 # noqa: E402
import kisch                    # noqa: E402
import coherence                # noqa: E402
import kicad2view               # noqa: E402
from CSXCAD.KiCadImport import KiCadBoard          # noqa: E402
from CSXCAD import FabStackups                      # noqa: E402

ENERGY = re.compile(r'Timestep:\s*(\d+).*?Speed:\s*([\d.]+)\s*MC/s.*?Energy:\s*~?\s*(\S+)\s*\(\s*-?\s*([\d.]+)dB\)')


class Cancelled(Exception):
    pass


class Job:
    def __init__(self, req, out):
        self.id = uuid.uuid4().hex[:10]
        self.req, self.out = req, out
        self.state, self.error, self.result = 'queued', None, None
        self.log = []
        self.progress = dict(exc=0, n=len(req.get('ports', [])), timestep=0, energy_db=0.0, mcells_s=None, stage='building the model')
        self.cancel = False
        self.t0 = time.time()
        self.touchstone = None

    def view(self, full=False):
        d = dict(id=self.id, state=self.state, error=self.error, progress=self.progress, elapsed_s=round(time.time() - self.t0, 1),
                 log=self.log[-12:], out=self.out)
        if self.result is not None and full:
            d['result'] = self.result
        return d


class Twin:
    """ Boards of one folder, a cache of parsed boards and one simulation at a time. """

    def __init__(self, root, threads=4):
        self.root = os.path.realpath(root)
        self.threads = threads
        self.boards = {}                 # (path, fab) -> (mtime, KiCadBoard)
        self.jobs = {}
        self.lock = threading.Lock()
        self.run_lock = threading.Lock()  # the engine uses all its threads: runs are queued

    # ------------------------------------------------------------ files
    def resolve(self, rel):
        p = os.path.realpath(os.path.join(self.root, rel))
        if not (p == self.root or p.startswith(self.root + os.sep)):
            raise PermissionError('outside the served folder')
        if not os.path.isfile(p):
            raise FileNotFoundError(rel)
        return p

    def list_boards(self):
        out = []
        for d, dirs, files in os.walk(self.root):
            depth = d[len(self.root):].count(os.sep)
            dirs[:] = [x for x in dirs if not x.startswith('.') and not x.endswith('_twin') and depth < 4]
            for f in files:
                if f.endswith('.kicad_pcb'):
                    p = os.path.join(d, f)
                    out.append(dict(path=os.path.relpath(p, self.root).replace(os.sep, '/'), size=os.path.getsize(p),
                                    mtime=os.path.getmtime(p)))
        return sorted(out, key=lambda b: b['path'])

    UPLOAD_EXT = ('.kicad_pcb', '.kicad_sch', '.kicad_pro', '.kicad_sym', '.kicad_dru')

    def upload(self, name, data, project=None):
        """ Store a file sent by the browser in <folder>/uploads/<project>/: the board, and the schematic
        sheets / project file of the same project next to it (the same name replaces the older copy). """
        base = os.path.basename(name or '')
        ext = os.path.splitext(base)[1].lower()
        if ext not in self.UPLOAD_EXT:
            raise ValueError('only KiCad board, schematic, project and symbol files can be uploaded')
        head = data.lstrip()[:16]
        if ext in ('.kicad_pcb', '.kicad_sch', '.kicad_sym') and not head.startswith(b'(kicad_' + ext[7:].encode()):
            if not (ext == '.kicad_sym' and head.startswith(b'(kicad_symbol_lib')):
                raise ValueError('%s does not look like a KiCad %s file' % (base, ext))
        base = re.sub(r'[^\w.,+-]', '_', base)
        proj = re.sub(r'[^\w.+-]', '_', project or os.path.splitext(base)[0])
        d = os.path.join(self.root, 'uploads', proj)
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, base)
        with open(p + '.part', 'wb') as f:
            f.write(data)
        os.replace(p + '.part', p)
        return dict(path=os.path.relpath(p, self.root).replace(os.sep, '/'), size=len(data))

    def schematic(self, rel):
        p = self.resolve(rel)
        root = kisch.find_root(p)
        if not root:
            return dict(found=False, folder=os.path.relpath(os.path.dirname(p), self.root).replace(os.sep, '/'))
        key = ('sch', p)
        mt = max(os.path.getmtime(os.path.join(os.path.dirname(root), f)) for f in os.listdir(os.path.dirname(root))
                 if f.endswith(('.kicad_sch', '.kicad_pcb')))
        with self.lock:
            hit = self.boards.get(key)
            if hit and hit[0] == mt:
                return hit[1]
        sc = kisch.Schematic(root, kisch.pad_nets_from_board(p)).scene()
        sc.update(found=True, root=os.path.basename(root))
        with self.lock:
            self.boards[key] = (mt, sc)
        return sc

    def coherence(self, rel, spec=None):
        p = self.resolve(rel)
        if spec is not None:
            coherence.save_spec(p, spec)
        spec = coherence.load_spec(p)
        return dict(spec=spec, result=coherence.compute(p, spec) if spec.get('groups') else None)

    def board(self, rel, fab=None):
        p = self.resolve(rel)
        key = (p, fab or '')
        mt = os.path.getmtime(p)
        with self.lock:
            hit = self.boards.get(key)
            if hit and hit[0] == mt:
                return p, hit[1]
        b = KiCadBoard(p)
        if fab:
            fn, nm = fab.split(':', 1)
            st = FabStackups.Find(fab=fn, name=nm)
            if not st:
                raise ValueError('stackup %s not found' % fab)
            FabStackups.Apply(b, st[0])
        with self.lock:
            self.boards[key] = (mt, b)
        return p, b

    # ------------------------------------------------------------ api
    def load(self, req):
        p, b = self.board(req['path'], req.get('fab'))
        name = os.path.splitext(os.path.basename(p))[0]
        data = kicad2view.export(b, name, board_dir=os.path.dirname(p), log=lambda *a: None)
        data['path'] = req['path']
        data['stackups'] = [dict(id='%s:%s' % (s.fab, s.name), fab=s.fab, name=s.name, total_mm=round(s.total_mm, 3), score=round(sc, 2),
                                 why=why, generic=s.generic) for sc, s, why in FabStackups.Suggest(b, n=8)]
        data['fab'] = req.get('fab')
        data['passives'] = [dict(ref=c['ref'], kind=c['kind'], value=c['value'], text=c['text'], layer=c['layer'], nets=list(c['nets']))
                            for c in b.Passives()]
        data['warnings'] = list(b.warnings)
        data['runs'] = self.runs(req['path'])
        kp = self.kicad_ports(p, b)
        if kp is not None:                      # the KiCad plugin's ports replace the viewer's suggestions
            data['ports'], data['port_source'] = kp, 'kicad'
        return data

    def field_overlay(self, req):
        """ One field picture for the KiCad plugin: {result: <run>/result.json, exc, k, quantity: 'J'|'V',
            layer, phase: degrees or None (= peak), image, lines}. Returns the RGBA heat map (base64 PNG, row 0
            = max y) over bbox, polylines (current lines: streamlines of the instantaneous surface current;
            voltage: iso-voltage lines) and the legend. """
        import base64
        import fields as fmod
        rp = self.resolve(req['result'])
        run = os.path.dirname(rp)
        idx = (json.load(open(rp, encoding='utf-8')).get('fields') or {})
        if idx.get('version', 1) < 2:
            raise ValueError('this run has no field data of version 2: simulate again')
        e, k, qn = int(req.get('exc', idx['excitations'][0])), int(req.get('k', 0)), req.get('quantity', 'J')
        g = idx['grid']
        names = [l['name'] for l in idx['layers']]
        mask = np.fromfile(os.path.join(run, idx['mask']), np.uint8).reshape(len(names), g['ny'], g['nx'])
        xs = g['x0'] + g['dx'] * np.arange(g['nx'])
        ys = g['y0'] + g['dx'] * np.arange(g['ny'])
        ph = req.get('phase')
        out = dict(bbox=idx['bbox'], f_GHz=idx['f_GHz'][k], exc=e, quantity=qn, power_W=idx.get('power_W', 1.0))
        if qn == 'J':
            li = names.index(req.get('layer') or names[0])
            J = fmod.load(run, idx, 'J', e, k)[li]                     # (ny, nx, 2)
            m = mask[li]
            peak = np.sqrt((np.abs(J) ** 2).sum(-1))
            ph = fmod.peak_phase(J) if ph is None else float(ph)
            vx, vy = fmod.instant(J[..., 0], ph), fmod.instant(J[..., 1], ph)
            mx = float(peak.max()) or 1.0
            rng = float(req.get('db_range', 40.0))
            img = fmod.colormap((20 * np.log10(np.maximum(peak, 1e-30) / mx) + rng) / rng)
            out.update(layer=names[li], phase=ph, legend=dict(kind='dB', unit='A/m', max=mx, range=rng,
                                                               text='|Js| peak %.3g A/m = 0 dB, %g dB range, 1 W incident at P%d' % (mx, rng, e)))
            if req.get('lines', True):
                dsep = float(req.get('dsep_mm', 0)) or max(g['dx'] * 6, (xs[-1] - xs[0]) / 120.0)
                out['dsep'] = dsep
                sl = fmod.streamlines(vx, vy, m, g['x0'], g['y0'], g['dx'], dsep_px=dsep / g['dx'], min_rel=float(req.get('min_rel', 0.05)))
                out['lines'] = [dict(pts=np.asarray(l['pts']).round(4).tolist(), mag=float(np.mean(l['mag'])) / mx) for l in sl]
        else:
            pairs = idx['V']['pairs']
            up = req.get('layer') or pairs[0][0]
            pi = [a for a, _ in pairs].index(up)
            V = fmod.load(run, idx, 'V', e, k)[pi]                     # (ny_e, nx_e) on the E grid
            ex, ey = np.asarray(idx['E']['x']), np.asarray(idx['E']['y'])
            Vg = fmod._resample2(V, ex, ey, xs, ys)                    # -> J grid
            m = mask[names.index(up)] if up in names else np.ones_like(mask[0])
            ph = fmod.peak_phase(Vg * m) if ph is None else float(ph)
            vi = fmod.instant(Vg, ph)
            mx = float(np.abs(Vg * m).max()) or 1.0
            img = fmod.colormap(0.5 + 0.5 * vi / mx, fmod._DIV)
            out.update(layer=up, against=pairs[pi][1], phase=ph,
                       legend=dict(kind='signed', unit='V', max=mx, text='V(%s - %s) at %.0f deg: -%.3g .. +%.3g V, 1 W incident at P%d'
                                   % (up, pairs[pi][1], ph, mx, mx, e)))
            if req.get('lines', True):
                lv = [mx * t for t in (-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)]
                inner = ~fmod._grow(m == 0, 2)                         # well inside the copper: lines across it
                cl = fmod.contour_lines(np.where(inner, vi, np.nan), xs, ys, lv)
                out['lines'] = [dict(pts=c['pts'], mag=c['level'] / mx) for c in cl]
        if req.get('image', True):
            a = np.where(m > 0, 235, 0).astype(np.uint8)[..., None]
            rgba = np.concatenate([img, a], -1)
            pad, px = g['dx'] / 2.0, g['dx']
            if req.get('burn') and out.get('lines'):
                # the lines drawn into a finer picture (white, black edge, chevrons): readable on any
                # map colour, whatever colour the user layer has in KiCad
                sc = int(req.get('scale') or max(1, min(4, 2400 // max(g['nx'], g['ny']))))
                rgba = fmod.burn_lines(rgba, out['lines'], g['x0'], g['y0'], g['dx'], sc, arrows=qn == 'J',
                                       dsep=out.get('dsep'))
                px = g['dx'] / sc
                out['lines_burned'] = True
            if req.get('row0', 'max') == 'max':                       # row 0 = max y (KiCad with the engine y up)
                rgba = rgba[::-1]
            tmp = os.path.join(run, 'fields', '_overlay.png')
            fmod.write_png(tmp, np.ascontiguousarray(rgba))
            out.update(png=base64.b64encode(open(tmp, 'rb').read()).decode('ascii'), width=int(rgba.shape[1]), height=int(rgba.shape[0]),
                       px=px, pad=pad)   # picture: x0 - pad .. x1 + pad, pixel pitch px
        return out

    def kicad_ports(self, p, b):
        """ Ports the KiCad plugin saved for this board (<board>.z2ems/setup.json, port_defs): None when
            there is no such file, else one viewer port per definition, made by port_from_spec so the
            3D markers show exactly what the engine simulates (gap ports: every sub-port). """
        sp = os.path.join(os.path.dirname(p), os.path.splitext(os.path.basename(p))[0] + '.z2ems', 'setup.json')
        if not os.path.exists(sp):
            return None
        try:
            defs = json.load(open(sp, encoding='utf-8')).get('port_defs') or []
        except Exception:
            return None
        out = []
        for d in defs:
            plus = d.get('plus') or {}
            mp, pl = d.get('minus_pads') or [], d.get('minus_plane')
            spec = {'plus': {'ref': plus.get('ref'), 'pad': plus.get('pad'), 'layer': plus.get('layer') or None},
                    'z0': float(d.get('z0') or 50.0), 'name': 'P%d' % int(d.get('n', len(out) + 1))}
            if mp:
                spec['minus'] = {'pads': [{'ref': m.get('ref'), 'pad': m.get('pad')} for m in mp]}
            elif pl:
                spec['minus'] = {'plane': pl}
            refs = {}
            for m in mp:
                refs.setdefault(m.get('ref'), set()).add(str(m.get('pad')))
            minus_txt = ('; '.join('%s pad %s' % (r, ', '.join(sorted(v))) for r, v in refs.items()) if mp
                         else ('plane %s' % pl if pl else 'nearest plane'))
            q = dict(ref=plus.get('ref'), net=plus.get('net', ''), layer=plus.get('layer'), pad_no=plus.get('pad'), z0=spec['z0'],
                     name=spec['name'], spec=spec, kicad=True, on=bool(d.get('excite', True)), confident=True, score=0,
                     minus=minus_txt, reasons=['defined in KiCad (Z2Labs EMS window)'])
            if not plus.get('ref'):
                q['error'] = 'no + terminal'
            else:
                try:
                    g = ks.port_from_spec(b, spec, len(out) + 1)
                    q.update(kind=g.get('kind', 'vertical'), start=g['start'], stop=g['stop'], ref_layer=g.get('ref_layer'),
                             subs=[dict(start=list(sb['start']), stop=list(sb['stop']), axis=sb['axis']) for sb in g.get('subs', [])])
                except Exception as e:
                    q['error'] = str(e)
            out.append(q)
        return out

    def snap(self, req):
        _, b = self.board(req['path'], req.get('fab'))
        pad, port = rf_ports.snap(b, float(req['x']), float(req['y']))
        if pad is None:
            return dict(port=None)
        return dict(port=dict(ref=pad.ref, net=pad.net, layer=pad.layer, start=port['start'], stop=port['stop'], ref_net=port['ref_net'],
                              ref_layer=port['ref_layer'], axis=port['normal_axis'], pad_no=pad.pad,
                              pad=[float(v) for v in (*pad.pts.min(0), *pad.pts.max(0))]))

    def _ports(self, b, specs):
        """ Port definitions from the plugin / viewer: each port is re-made on the server from its pads
            (kicad_sim.port_from_spec: + pad and - pads / plane), so the geometry always matches the board
            the engine sees. """
        return [ks.port_from_spec(b, sp, i + 1) for i, sp in enumerate(specs)]

    def mesh(self, req):
        """ Mesh preview: the lines build() would hand to the engine for these ports and fmax, in board
            millimetres, plus the quality figures. Used by the KiCad plugin's Z2EMS.Mesh overlay. """
        p, b = self.board(req['path'], req.get('fab'))
        ports = self._ports(b, req['ports'])
        if not ports:
            raise ValueError('no ports')
        fmax = float(req.get('fmax', 6e9))
        # same element-end lines as a run with materials.components (build() -> ToCSX -> AddLumpedComponents),
        # so the preview is the simulated mesh also when R/L/C parts are modelled
        b.mesh_hints = None
        if (req.get('materials') or {}).get('components'):
            from CSXCAD import ContinuousStructure
            b.AddLumpedComponents(ContinuousStructure())
        X, Y, Z, res = ks.final_mesh(b, ports, fmax, coarse=float(req.get('coarse', 1.0)))
        X, Y, Z = (sorted(set(float(ks._sig7(v)) for v in a)) for a in (X, Y, Z))
        er = max([l.epsilon_r or 1.0 for l in b.stackup if l.type in ('prepreg', 'core')] + [1.0])
        (bx0, bx1), (by0, by1) = b.BoundingBox()[:2]
        return dict(x=X, y=Y, z=Z, res_mm=res, bbox=[bx0, by0, bx1, by1],
                    layers_z=[dict(layer=c, z=b.SheetZ(c)) for c in b.copper_layers],
                    ports=[dict(number=q['number'], ref=q.get('ref'), net=q.get('net'), kind=q.get('kind', 'vertical'),
                                start=list(map(float, q['start'])), stop=list(map(float, q['stop'])), axis=q['normal_axis'],
                                ref_layer=q.get('ref_layer'), ref_net=q.get('ref_net'),
                                subs=[dict(start=list(map(float, g['start'])), stop=list(map(float, g['stop'])), axis=g['axis'],
                                           minus=g['minus']) for g in q.get('subs', [])]) for q in ports],
                    quality=ks.mesh_quality([X, Y, Z], fmax, er, res))

    def run(self, req):
        p, b = self.board(req['path'], req.get('fab'))
        ports = self._ports(b, req['ports'])
        if not ports:
            raise ValueError('no ports')
        name = os.path.splitext(os.path.basename(p))[0]
        out = os.path.join(os.path.dirname(p), name + '_twin', time.strftime('%Y%m%d-%H%M%S'))
        job = Job(req, out)
        job.ports = ports
        self.jobs[job.id] = job
        threading.Thread(target=self._work, args=(job, p), daemon=True).start()
        return job.view()

    def _work(self, job, path):
        with self.run_lock:
            if job.cancel:
                job.state = 'cancelled'
                return
            job.state = 'running'
            r = job.req
            fmax = float(r.get('fmax', 6e9))
            nf = int(r.get('nf', 400))
            fmin = float(r.get('fmin') or fmax / nf)          # a uniform grid from df: the result also gives a TDR
            mat = dict(r.get('materials') or {})

            def log(m):
                job.log.append(time.strftime('%H:%M:%S ') + m)
                if m.lstrip().startswith('fields:') or 'rc=' in m:
                    job.progress['stage'] = 'post-processing'

            def progress(e, n, line):
                if job.cancel:
                    raise Cancelled()
                if job.progress.get('exc') != e:      # a new excitation starts from zero (the bar showed the last one's 100 %)
                    job.progress.update(timestep=0, energy_db=0.0, frac=0.0, mcells_s=None, max_ts=0, stage='engine setup')
                job.progress['exc'], job.progress['n'] = e, n
                m = ENERGY.search(line)
                if m:
                    job.progress.update(timestep=int(m.group(1)), mcells_s=float(m.group(2)), energy_db=-float(m.group(4)),
                                        stage='time stepping')
                    mx = job.progress.get('max_ts') or 0
                    # done when the energy is down to the end criterion (-40 dB) or the time limit is reached
                    job.progress['frac'] = round(max(min(1.0, job.progress['energy_db'] / -40.0),
                                                     min(1.0, job.progress['timestep'] / mx) if mx else 0.0), 3)
                elif 'Max. number of timesteps' in line:
                    mm = re.search(r'timesteps:\s*(\d+)', line)
                    if mm:
                        job.progress['max_ts'] = int(mm.group(1))

            try:
                os.makedirs(job.out, exist_ok=True)
                json.dump(dict(r, started=time.strftime('%Y-%m-%d %H:%M:%S')), open(os.path.join(job.out, 'request.json'), 'w'), indent=1)
                fl = r.get('fields', True)
                fields = None
                if fl:
                    fq = fl.get('freqs') if isinstance(fl, dict) else None
                    fields = dict(freqs=fq or [fmax * 0.25, fmax * 0.5, fmax * 0.9])
                eng = r.get('engine', 'auto')          # auto | cpu | cuda
                variant = {'cpu': ks.DEFAULT_VARIANT, 'cuda': 'cuda'}.get(eng, 'auto')
                res = ks.simulate(path, None, fmin, fmax, job.out, fab=r.get('fab'), variant=variant, threads=int(r.get('threads') or self_threads()),
                                  coarse=float(r.get('coarse', 1.0)), nf=nf, log=log, materials=mat, ports=job.ports,
                                  progress=progress, fields=fields, converge=bool(r.get('converge')),
                                  excite=r.get('excite') or None)
                if res.get('fields') and res['fields'].get('maps'):          # map paths relative to the served folder
                    res['fields']['base'] = os.path.relpath(job.out, self.root).replace(os.sep, '/')
                    json.dump(res, open(os.path.join(job.out, 'result.json'), 'w'), indent=1)
                job.result = res
                n = len(job.ports)
                tp = [f for f in os.listdir(job.out) if f.endswith('.s%dp' % n)]
                job.touchstone = os.path.join(job.out, tp[0]) if tp else None
                job.state = 'done'
                log('done in %.0f s' % (time.time() - job.t0))
            except Cancelled:
                job.state = 'cancelled'
                log('cancelled')
            except Exception as e:
                job.state, job.error = 'error', '%s: %s' % (type(e).__name__, e)
                log(traceback.format_exc(limit=3))

    def runs(self, rel):
        p = self.resolve(rel)
        d = os.path.join(os.path.dirname(p), os.path.splitext(os.path.basename(p))[0] + '_twin')
        out = []
        if os.path.isdir(d):
            for k in sorted(os.listdir(d), reverse=True)[:20]:
                rj = os.path.join(d, k, 'result.json')
                if os.path.isfile(rj):
                    out.append(dict(run=k, result=os.path.relpath(rj, self.root).replace(os.sep, '/')))
        return out


_TWIN = None


def self_threads():
    return _TWIN.threads if _TWIN else 4


TOKEN = None          # access token when the server listens beyond this computer
MAX_UPLOAD = 64 << 20


def _authorized(handler):
    """ True when no token is required, or the request carries it (cookie, bearer header or ?token=). """
    if not TOKEN:
        return True
    cands = [handler.headers.get('Authorization', '')[7:].strip()]
    for part in (handler.headers.get('Cookie') or '').split(';'):
        k, _, v = part.strip().partition('=')
        if k == 'z2twin':
            cands.append(v)
    cands += parse_qs(urlparse(handler.path).query).get('token', [])
    return any(c and hmac.compare_digest(c, TOKEN) for c in cands)


_CUDA_PROBE = None


def _cuda_state():
    """ ks.cuda_available() runs the engine's --help (seconds on a cold disk); /api/health must answer at once
        (the plugin waits 3 s per port when it looks for a server), so the probe runs in the background and
        health reports None until it is done. """
    global _CUDA_PROBE
    if ks._CUDA_OK is not None:
        return ks._CUDA_OK
    if _CUDA_PROBE is None:
        _CUDA_PROBE = threading.Thread(target=ks.cuda_available, daemon=True)
        _CUDA_PROBE.start()
    return None


class Handler(BaseHTTPRequestHandler):
    server_version = 'z2twin/1'

    def _gate(self):
        """ Token check; a valid ?token= is moved into a cookie and the URL cleaned. Returns True to go on. """
        if not _authorized(self):
            body = (b'<!doctype html><meta name="viewport" content="width=device-width"><title>z2twin</title>'
                    b'<p style="font:16px system-ui;margin:2em">Access token missing or wrong. Open the URL printed by '
                    b'the twin server (it ends with <code>?token=...</code>).</p>')
            self._send(401, body, 'text/html; charset=utf-8')
            return False
        u = urlparse(self.path)
        q = parse_qs(u.query)
        if TOKEN and 'token' in q and self.command == 'GET':
            rest = '&'.join('%s=%s' % (k, v[0]) for k, v in q.items() if k != 'token')
            self._send(303, b'', 'text/plain', {'Location': u.path + ('?' + rest if rest else ''),
                                                'Set-Cookie': 'z2twin=%s; Path=/; HttpOnly; SameSite=Strict; Max-Age=2592000' % TOKEN})
            return False
        return True

    def log_message(self, fmt, *args):            # quiet
        pass

    def _send(self, code, body, ctype='application/json', headers=None):
        data = body if isinstance(body, bytes) else json.dumps(body, separators=(',', ':')).encode('utf-8')
        gz = len(data) > 32768 and 'gzip' in (self.headers.get('Accept-Encoding') or '') and \
            (ctype.startswith('application/json') or ctype.startswith('text/') or ctype == 'application/octet-stream')
        if gz:                                 # board JSON, field arrays (zero off the copper): 5-10x smaller
            data = gzip.compress(data, 5)
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        if gz:
            self.send_header('Content-Encoding', 'gzip')
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Cache-Control', 'no-store')
        for k, v in (headers or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(data)

    def _err(self, code, msg):
        self._send(code, dict(error=msg))

    def _body(self):
        n = int(self.headers.get('Content-Length') or 0)
        return json.loads(self.rfile.read(n) or b'{}') if n else {}

    def do_GET(self):
        if not self._gate():
            return
        u = urlparse(self.path)
        q = parse_qs(u.query)
        try:
            if u.path == '/api/health':
                eng = ks.VARIANTS[ks.DEFAULT_VARIANT]['exe']
                return self._send(200, dict(ok=True, folder=_TWIN.root, engine=eng, engine_found=os.path.exists(eng) or bool(eng and os.sep not in eng),
                                            variant=ks.DEFAULT_VARIANT, threads=_TWIN.threads,
                                            cuda=_cuda_state(), cuda_engine=ks.VARIANTS['cuda']['exe'],
                                            ports_api=2,      # 2: + pad / - pads or plane (kicad_sim.port_from_spec)
                                            fields_api=2,     # 2: Js, E volume, layer voltages, /api/field_overlay
                                            # who this server is: a client (the KiCad plugin) replaces a server of
                                            # another engine folder or version, since a server outlives its client
                                            root=ENGINE_ROOT, version=ENGINE_VERSION, pid=os.getpid(),
                                            busy=sum(1 for j in list(_TWIN.jobs.values()) if j.state in ('queued', 'running'))))
            if u.path == '/api/boards':
                return self._send(200, dict(boards=_TWIN.list_boards(), folder=_TWIN.root))
            if u.path == '/api/runs':
                return self._send(200, dict(runs=_TWIN.runs(q['path'][0])))
            if u.path == '/api/schematic':
                return self._send(200, _TWIN.schematic(q['path'][0]))
            if u.path == '/api/coherence':
                return self._send(200, _TWIN.coherence(q['path'][0]))
            if u.path.startswith('/api/file'):
                p = _TWIN.resolve(q['path'][0])
                return self._send(200, open(p, 'rb').read(), mimetypes.guess_type(p)[0] or 'application/octet-stream')
            m = re.match(r'^/api/job/(\w+)(/touchstone)?$', u.path)
            if m:
                job = _TWIN.jobs.get(m.group(1))
                if not job:
                    return self._err(404, 'no such job')
                if m.group(2):
                    if not job.touchstone:
                        return self._err(404, 'no result yet')
                    return self._send(200, open(job.touchstone, 'rb').read(), 'text/plain',
                                      {'Content-Disposition': 'attachment; filename="%s"' % os.path.basename(job.touchstone)})
                return self._send(200, job.view(full=True))
            return self._static(u.path)
        except (PermissionError, FileNotFoundError) as e:
            return self._err(404, str(e))
        except Exception as e:
            return self._err(500, '%s: %s' % (type(e).__name__, e))

    def do_POST(self):
        if not self._gate():
            return
        u = urlparse(self.path)
        try:
            if u.path == '/api/upload':
                n = int(self.headers.get('Content-Length') or 0)
                if n <= 0 or n > MAX_UPLOAD:
                    return self._err(413, 'board file missing or larger than %d MB' % (MAX_UPLOAD >> 20))
                qq = parse_qs(u.query)
                return self._send(200, _TWIN.upload(qq.get('name', [''])[0], self.rfile.read(n), qq.get('project', [None])[0]))
            if u.path == '/api/shutdown':
                if TOKEN is None and self.client_address[0] not in ('127.0.0.1', '::1'):
                    return self._err(403, 'shutdown only from this computer')
                self._send(200, dict(ok=True, pid=os.getpid()))
                threading.Thread(target=_SRV.shutdown, daemon=True).start()
                return
            req = self._body()
            if u.path == '/api/load':
                return self._send(200, _TWIN.load(req))
            if u.path == '/api/snap':
                return self._send(200, _TWIN.snap(req))
            if u.path == '/api/run':
                return self._send(200, _TWIN.run(req))
            if u.path == '/api/mesh':
                return self._send(200, _TWIN.mesh(req))
            if u.path == '/api/field_overlay':
                return self._send(200, _TWIN.field_overlay(req))
            if u.path == '/api/coherence':
                return self._send(200, _TWIN.coherence(req['path'], req.get('spec')))
            if u.path == '/api/pad':
                p = _TWIN.resolve(req['path'])
                return self._send(200, dict(pad=coherence.pad_near(p, float(req['x']), -float(req['y']))))
            m = re.match(r'^/api/job/(\w+)/cancel$', u.path)
            if m and m.group(1) in _TWIN.jobs:
                _TWIN.jobs[m.group(1)].cancel = True
                return self._send(200, dict(ok=True))
            return self._err(404, 'unknown endpoint')
        except (PermissionError, FileNotFoundError) as e:
            return self._err(404, str(e))
        except (ValueError, KeyError) as e:
            return self._err(400, '%s: %s' % (type(e).__name__, e))
        except Exception as e:
            return self._err(500, '%s: %s' % (type(e).__name__, e))

    def _static(self, path):
        rel = path.lstrip('/') or 'index.html'
        if rel == 'index.html' and os.path.exists(os.path.join(VIEWER, 'index.local.html')) and os.path.isdir(os.path.join(VIEWER, 'three')):
            rel = 'index.local.html'           # offline three.js copy when present
        p = os.path.normpath(os.path.join(VIEWER, rel))      # (symlinks inside the viewer folder are followed)
        if not p.startswith(os.path.normpath(VIEWER) + os.sep) or not os.path.isfile(p):
            return self._err(404, 'not found')
        ctype = mimetypes.guess_type(p)[0] or 'application/octet-stream'
        if p.endswith('.js'):
            ctype = 'text/javascript'
        return self._send(200, open(p, 'rb').read(), ctype)


def network_addresses():
    """ (label, IPv4) of this computer's network interfaces, the Tailscale address marked. """
    out, seen = [], set()
    try:
        ts = subprocess.run(['tailscale', 'ip', '-4'], capture_output=True, text=True, timeout=3).stdout.split()
    except Exception:
        ts = []
    ips = set(ts)
    try:
        ips.update(i[4][0] for i in socket.getaddrinfo(socket.gethostname(), None, socket.AF_INET))
    except OSError:
        pass
    try:                                       # Linux: every interface (LAN, ZeroTier, VPNs)
        ips.update(subprocess.run(['hostname', '-I'], capture_output=True, text=True, timeout=3).stdout.split())
    except Exception:
        pass
    try:                                       # the address of the default route
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(('10.255.255.255', 1))
        ips.add(s.getsockname()[0])
        s.close()
    except OSError:
        pass
    for ip in sorted(ips):
        if ip.startswith('127.') or ':' in ip or ip in seen:
            continue
        seen.add(ip)
        out.append(('tailscale' if ip in ts or ip.startswith('100.') else 'network', ip))
    return out


ENGINE_VERSION = '0.3.2'


class Server(ThreadingHTTPServer):
    # http.server sets SO_REUSEADDR, and on Windows that lets a second server bind a port that is
    # already served: two servers answered on 8765 after a double click (2026-11-11). One port, one server.
    allow_reuse_address = os.name != 'nt'
    daemon_threads = True

    def server_bind(self):
        if os.name == 'nt':
            import socket
            self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        super().server_bind()

    def handle_error(self, request, client_address):
        import sys
        if isinstance(sys.exc_info()[1], (ConnectionError, TimeoutError)):
            return                                  # the browser or the plugin closed the connection: not an error
        super().handle_error(request, client_address)
ENGINE_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
_SRV = None


def main(argv=None):
    global _TWIN
    ap = argparse.ArgumentParser(description='Local back end of the z2labs board viewer (openEMS twin).')
    ap.add_argument('folder', help='folder with .kicad_pcb boards (only files below it are served)')
    ap.add_argument('--port', type=int, default=8765)
    ap.add_argument('--threads', type=int, default=max(1, min(8, (os.cpu_count() or 4) // 2)))
    ap.add_argument('--open', action='store_true', help='open the viewer in the browser')
    ap.add_argument('--host', default='127.0.0.1', help='address to listen on (default: this computer only)')
    ap.add_argument('--lan', action='store_true', help='listen on all networks (same as --host 0.0.0.0), with a token')
    ap.add_argument('--token', default=os.environ.get('Z2TWIN_TOKEN'),
                    help='access token (default: a new random one whenever the server listens beyond 127.0.0.1)')
    a = ap.parse_args(argv)
    global TOKEN
    host = '0.0.0.0' if a.lan else a.host
    local_only = host in ('127.0.0.1', 'localhost', '::1')
    TOKEN = a.token or (None if local_only else secrets.token_urlsafe(18))
    _TWIN = Twin(a.folder, a.threads)
    global _SRV
    try:
        srv = _SRV = Server((host, a.port), Handler)
    except OSError as e:
        print('z2twin: port %d is in use (%s): another server is already there' % (a.port, e), flush=True)
        return 2
    suffix = '?token=' + TOKEN if TOKEN else ''
    url = 'http://127.0.0.1:%d/%s' % (a.port, suffix)
    print('z2twin: serving %s (engine %s%s, %d threads)' % (_TWIN.root, ks.DEFAULT_VARIANT,
          ', GPU build found' if ks.cuda_available() else '', a.threads), flush=True)
    print('  this computer:  %s' % url, flush=True)
    if not local_only:
        for name, ip in network_addresses():
            print('  %-14s  http://%s:%d/%s' % (name + ':', ip, a.port, suffix), flush=True)
        print('  The token gives full access to the boards of this folder; share the URL only with your own devices.', flush=True)
    if a.open:
        webbrowser.open(url)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
