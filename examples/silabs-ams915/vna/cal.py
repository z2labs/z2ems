#!/usr/bin/env python3
"""Software SOL (1-port) calibration for the NanoVNA V2 raw data + corrected DUT measurement.

  python cal.py calibrate --port COM6 --name cal_2026-11-03_sma_f --start 100e6 --stop 3e9 --points 1001
      -> prompts OPEN / SHORT / LOAD, stores vna/cal/<name>/{open,short,load}_raw.s1p + cal.json (error terms)
  python cal.py measure --port COM6 --cal vna/cal/<name> --dut WES0111 --out meas/WES0111_standalone.s1p [--repeat 3] [--note "..."]
      -> stores <out> (corrected S11) and <out>.raw.s1p, plus a JSON sidecar with provenance.

Standards are ideal by default (gamma = +1 / -1 / 0). For a kit with known offset delays pass
--open-delay/--short-delay (ps, one-way): gamma = +/-exp(-j*4*pi*f*delay).
Error model (one-port, 3-term): Gm = e00 + e01e10*Ga/(1-e11*Ga)  -> solved from the three standards.
"""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nanovna_v2 import NanoVNAV2, write_touchstone

def read_s1p(path):
    f, s = [], []
    for l in open(path):
        if l.startswith(('!', '#')) or not l.strip(): continue
        v = l.split(); f.append(float(v[0])); s.append(complex(float(v[1]), float(v[2])))
    return np.array(f), np.array(s)

def std_gamma(kind, f, delay_ps=0.0):
    ph = np.exp(-1j * 4 * np.pi * f * delay_ps * 1e-12)
    return {'open': ph, 'short': -ph, 'load': np.zeros_like(f, dtype=complex)}[kind]

def solve_terms(f, m_open, m_short, m_load, g_open, g_short, g_load):
    """Per frequency solve Gm*(1) + Ga*Gm*(-e11)... linear system in (e00, e11, delta) with delta = e00*e11 - e01e10:
       Gm = e00 - delta*Ga + e11*Ga*Gm  ->  [1, -Ga, Ga*Gm] . [e00, delta, e11] = Gm"""
    e00 = np.zeros_like(f, dtype=complex); e11 = e00.copy(); delta = e00.copy()
    for i in range(len(f)):
        A = np.array([[1, -g_open[i], g_open[i] * m_open[i]], [1, -g_short[i], g_short[i] * m_short[i]], [1, -g_load[i], g_load[i] * m_load[i]]])
        x = np.linalg.solve(A, np.array([m_open[i], m_short[i], m_load[i]]))
        e00[i], delta[i], e11[i] = x
    return e00, e11, delta

def correct(gm, e00, e11, delta):
    return (gm - e00) / (gm * e11 - delta)

def sweep(v, start, stop, points, avg, repeat):
    v.set_sweep(start, stop, points, avg)
    runs = [v.sweep() for _ in range(repeat)]
    f = runs[0][0]
    s11 = np.median([r[1] for r in runs], axis=0) if repeat > 1 else runs[0][1]
    return f, s11

def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd')
    c = sub.add_parser('calibrate'); c.add_argument('--port', required=True); c.add_argument('--name', required=True)
    c.add_argument('--start', type=float, default=100e6); c.add_argument('--stop', type=float, default=3e9)
    c.add_argument('--points', type=int, default=1001); c.add_argument('--avg', type=int, default=4); c.add_argument('--repeat', type=int, default=3)
    c.add_argument('--open-delay', type=float, default=0.0); c.add_argument('--short-delay', type=float, default=0.0)
    m = sub.add_parser('measure'); m.add_argument('--port', required=True); m.add_argument('--cal', required=True)
    m.add_argument('--dut', required=True); m.add_argument('--out', required=True); m.add_argument('--repeat', type=int, default=3)
    m.add_argument('--avg', type=int, default=4); m.add_argument('--note', default='')
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    if a.cmd == 'calibrate':
        d = os.path.join(here, 'cal', a.name); os.makedirs(d, exist_ok=True)
        v = NanoVNAV2(a.port); raw = {}
        for k in ('open', 'short', 'load'):
            input('Connect %s to port 1, then Enter: ' % k.upper())
            f, s = sweep(v, a.start, a.stop, a.points, a.avg, a.repeat); raw[k] = s
            write_touchstone(os.path.join(d, k + '_raw.s1p'), f, s, comment='%s standard raw, avg=%d repeat=%d' % (k, a.avg, a.repeat))
            print('  %s: mean |raw| = %.3f' % (k, np.mean(abs(s))))
        e00, e11, delta = solve_terms(f, raw['open'], raw['short'], raw['load'],
                                      std_gamma('open', f, a.open_delay), std_gamma('short', f, a.short_delay), std_gamma('load', f))
        # self-check: corrected load must be ~0, corrected open ~ +1
        chk = {k: float(np.max(abs(correct(raw[k], e00, e11, delta) - std_gamma(k, f, {'open': a.open_delay, 'short': a.short_delay}.get(k, 0))))) for k in raw}
        json.dump(dict(name=a.name, created=time.strftime('%Y-%m-%d %H:%M'), device=v.info(), start=a.start, stop=a.stop, points=a.points,
                       avg=a.avg, repeat=a.repeat, open_delay_ps=a.open_delay, short_delay_ps=a.short_delay, f=f.tolist(),
                       e00=[[z.real, z.imag] for z in e00], e11=[[z.real, z.imag] for z in e11], delta=[[z.real, z.imag] for z in delta],
                       selfcheck_max_residual=chk), open(os.path.join(d, 'cal.json'), 'w'))
        print('saved', d, 'self-check residuals', chk); v.close()
    else:
        cal = json.load(open(os.path.join(a.cal, 'cal.json')))
        f0 = np.array(cal['f']); E = lambda k: np.array([complex(*z) for z in cal[k]])
        v = NanoVNAV2(a.port)
        f, s = sweep(v, cal['start'], cal['stop'], cal['points'], a.avg, a.repeat)
        assert np.allclose(f, f0)
        sc = correct(s, E('e00'), E('e11'), E('delta')); info = v.info(); v.close()
        os.makedirs(os.path.dirname(os.path.abspath(a.out)) or '.', exist_ok=True)
        write_touchstone(a.out, f, sc, comment='%s corrected (SOL %s) %s %s' % (a.dut, cal['name'], time.strftime('%Y-%m-%d %H:%M'), a.note))
        write_touchstone(a.out + '.raw.s1p', f, s, comment='%s raw %s' % (a.dut, a.note))
        i915 = np.argmin(abs(f - 915e6)); z = 50 * (1 + sc[i915]) / (1 - sc[i915])
        json.dump(dict(dut=a.dut, cal=cal['name'], time=time.strftime('%Y-%m-%d %H:%M'), device=info, avg=a.avg, repeat=a.repeat, note=a.note,
                       s11_dB_915=float(20 * np.log10(abs(sc[i915]))), z_915=[z.real, z.imag]), open(a.out + '.json', 'w'), indent=1)
        print('%s: |S11|(915 MHz) = %.1f dB, Z = %.1f %+.1fj ohm -> %s' % (a.dut, 20 * np.log10(abs(sc[i915])), z.real, z.imag, a.out))
if __name__ == '__main__': main()
