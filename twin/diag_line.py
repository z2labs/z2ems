#!/usr/bin/env python3
""" Diagnostic for the zig-zag / diagonal track divergence: a bare microstrip on the JLC 4-layer stackup (F.Cu trace
    over In1 ground), lumped ports at both ends, PEC copper, PML x 6 — no connector model, no KiCad import.
    Variants:  straight | diag10 (one straight line at 10 deg) | zigzag10 | diag45 | stair10 (10 deg built from
    axis-aligned 0.1 mm steps).  Reports |S11| max and whether the energy grows.

        python diag_line.py <out> [--variants straight,diag10,zigzag10,diag45,stair10] [--fmax 10e9] [--res 0.1]
"""
import json, math, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kicad_sim as ks

H_SUB, ER, W = 0.2104, 4.4, 0.3726
L, MARGIN = 30.0, 5.0


def build(variant, fmax, res_mm=None):
    ks._setup_python()
    from CSXCAD import ContinuousStructure
    from openEMS import openEMS
    FDTD = openEMS(EndCriteria=1e-4, NrTS=300000, MaxTime=40e-9)
    FDTD.SetGaussExcite(0.5 * fmax + 0.15e9, 0.5 * fmax - 0.15e9)
    FDTD.SetBoundaryCond(['PML_8'] * 6)
    CSX = ContinuousStructure(); FDTD.SetCSX(CSX)
    sub = CSX.AddMaterial('sub', epsilon=ER)
    sub.AddBox([-2, -8, 0], [L + 2, 8, H_SUB], priority=1)
    gnd = CSX.AddMetal('gnd'); gnd.AddBox([-2, -8, 0], [L + 2, 8, 0], priority=5)
    tr = CSX.AddMetal('trace')
    y_in, y_out = 0.0, 0.0
    if variant == 'straight':
        pts = [(0, 0), (L, 0)]
    elif variant == 'diag10':
        y_out = L * math.tan(math.radians(10)); pts = [(0, 0), (L, y_out)]
    elif variant == 'diag45':
        y_out = 6.0; pts = [(0, 0), (6.0, 6.0), (L, 6.0)]
    elif variant == 'zigzag10':
        pts = [(0, 0), (L / 2, -(L / 2) * math.tan(math.radians(10))), (L, 0)]
    elif variant == 'stair10':                          # 10 deg from 0.1 mm axis-aligned steps
        pts, x, y = [(0, 0)], 0.0, 0.0
        step = 0.1
        while x < L:
            x2 = min(x + step / math.tan(math.radians(10)), L); pts.append((x2, y)); y += step; pts.append((x2, y)); x = x2
        y_out = y
    else:
        raise SystemExit(variant)
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):        # each segment as a LinPoly strip of width W
        dx, dy = x1 - x0, y1 - y0; n = math.hypot(dx, dy); ux, uy = -dy / n * W / 2, dx / n * W / 2
        px = [x0 + ux, x1 + ux, x1 - ux, x0 - ux]; py = [y0 + uy, y1 + uy, y1 - uy, y0 - uy]
        tr.AddPolygon([px, py], 'z', H_SUB, priority=10)             # points = [x list, y list]
    res = res_mm or ks.C0 / (fmax * math.sqrt(ER)) / 1e-3 / 20
    X = np.concatenate([np.arange(-MARGIN - 2, L + MARGIN + 2.01, res), [0, L, 0 + 0.3, L - 0.3]])
    Y = np.concatenate([np.arange(-8 - MARGIN, 8 + MARGIN + 0.01, res), [-W / 2, W / 2, y_out - W / 2, y_out + W / 2, y_in]])
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        X = np.concatenate([X, [x0, x1]]); Y = np.concatenate([Y, [y0, y1]])
    Z = np.concatenate([np.arange(-MARGIN, H_SUB + MARGIN + 0.01, res), np.linspace(0, H_SUB, 4)])
    g = CSX.GetGrid(); g.SetDeltaUnit(1e-3)
    for ax, lines in zip('xyz', (X, Y, Z)):
        g.SetLines(ax, np.unique(np.round(ks.smooth(lines, res, merge=0.03), 5)))
    p1 = FDTD.AddLumpedPort(1, 50, [0 + 0.3, -W / 2, 0], [0 + 0.3, W / 2, H_SUB], 'z', excite=1, priority=20)
    p2 = FDTD.AddLumpedPort(2, 50, [L - 0.3, y_out - W / 2, 0], [L - 0.3, y_out + W / 2, H_SUB], 'z', excite=0, priority=20)
    return FDTD, CSX, [p1, p2]


def main():
    out = sys.argv[1]
    variants = sys.argv[sys.argv.index('--variants') + 1].split(',') if '--variants' in sys.argv else ['straight', 'diag10', 'zigzag10', 'diag45', 'stair10']
    fmax = float(sys.argv[sys.argv.index('--fmax') + 1]) if '--fmax' in sys.argv else 10e9
    res_mm = float(sys.argv[sys.argv.index('--res') + 1]) if '--res' in sys.argv else None     # cell size (mm) instead of lambda/20
    os.makedirs(out, exist_ok=True)
    rows = []
    for v in variants:
        d = os.path.join(out, v); os.makedirs(d, exist_ok=True)
        FDTD, CSX, ports = build(v, fmax, res_mm)
        xml = os.path.join(d, 'model.xml'); FDTD.Write2XML(xml)
        nl = [len(CSX.GetGrid().GetLines(a)) for a in 'xyz']
        r = ks.run_engine(xml, os.path.join(d, 'sim'), os.environ.get('Z2_VARIANT', 'cuda'), int(os.environ.get('Z2_THREADS', '4')))
        f = np.linspace(0.2e9, fmax, 491)
        for p in ports:
            p.CalcPort(os.path.join(d, 'sim'), f, ref_impedance=50)
        s11 = 20 * np.log10(np.abs(ports[0].uf_ref / ports[0].uf_inc) + 1e-15)
        s21 = 20 * np.log10(np.abs(ports[1].uf_ref / ports[0].uf_inc) + 1e-15)
        import runviz
        pr = runviz.progress(os.path.join(d, 'sim', 'stdout.txt'))
        row = dict(variant=v, res_mm=res_mm, mesh=nl, mcells=float(np.prod(nl)) / 1e6, rc=r['rc'], wall_s=r['wall_s'], ts=r['timesteps'],
                   s11_max_dB=float(s11.max()), s21_min_dB=float(s21.min()), diverging=pr['diverging'], energy_last=pr['energy'][-1] if pr['energy'] else None)
        rows.append(row); print(json.dumps(row), flush=True)
        json.dump(rows, open(os.path.join(out, 'diag.json'), 'w'), indent=1)
    print('DIAG_DONE')


if __name__ == '__main__':
    main()
