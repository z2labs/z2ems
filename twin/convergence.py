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
convergence: is a result mesh-converged? Two runs of the same model at the base resolution and
at a finer one (cells scaled by `refine`, 0.7 = about 3x the cell count); the differences of the
things the twin reports decide:

    f_res      resonance / S11 minimum frequency        tolerance 0.5 %
    tau        S21 delay (phase slope) per path          tolerance 2 ps
    S          complex S11 over the band, linear          tolerance 0.03 (max |dS|)

`compare(base, fine)` takes the result dicts the runners write (f_GHz, s11_re/im or S arrays, tau_ps)
and returns dict(converged, df_pct, dtau_ps, dS, verdict). The verdict is text for the report.
"""
import numpy as np

TOL = dict(f_pct=0.5, tau_ps=2.0, S_lin=0.03)


def _s11(r):
    if 's11_re' in r:
        return np.asarray(r['f_GHz'], float), np.asarray(r['s11_re'], float) + 1j * np.asarray(r['s11_im'], float)
    if 'S_re' in r and 'S11' in r['S_re']:
        return np.asarray(r['f_GHz'], float), np.asarray(r['S_re']['S11'], float) + 1j * np.asarray(r['S_im']['S11'], float)
    return None, None


def _tau(r):
    """ S21 delay of a two-port board result (kicad_sim.simulate), None otherwise """
    if 'S_re' in r and 'S21' in r['S_re'] and r.get('tau_ps') is None:
        f = np.asarray(r['f_GHz'], float)
        s = np.asarray(r['S_re']['S21'], float) + 1j * np.asarray(r['S_im']['S21'], float)
        k = (f >= 0.15 * f.max()) & (f <= 0.6 * f.max())
        return delay_ps(f[k], s[k])
    return r.get('tau_ps')


def compare(base, fine, tol=None):
    t = dict(TOL, **(tol or {}))
    out = dict(converged=True, checks=[])
    fb = base.get('f_res_GHz'); ff = fine.get('f_res_GHz')
    if fb and ff:
        d = 100.0 * (ff - fb) / fb
        out['df_pct'] = round(d, 3); ok = abs(d) <= t['f_pct']
        out['checks'].append(('f_res', '%.4f -> %.4f GHz (%+.2f %%)' % (fb, ff, d), ok)); out['converged'] &= ok
    tb = _tau(base); tf = _tau(fine)
    if tb is not None and tf is not None:
        tb = np.atleast_1d(tb); tf = np.atleast_1d(tf)
        d = float(np.max(np.abs(tf - tb)))
        out['dtau_ps'] = round(d, 2); ok = d <= t['tau_ps']
        out['checks'].append(('delay', 'max |dtau| %.2f ps' % d, ok)); out['converged'] &= ok
    f1, s1 = _s11(base); f2, s2 = _s11(fine)
    if s1 is not None and s2 is not None:
        lo, hi = max(f1.min(), f2.min()), min(f1.max(), f2.max()); f = np.linspace(lo, hi, 400)
        # complex difference of the reflection, linear: 0.03 is -30 dB, a dB comparison would blow up
        # the noise floor of a well matched line
        a = np.interp(f, f1, s1.real) + 1j * np.interp(f, f1, s1.imag); b = np.interp(f, f2, s2.real) + 1j * np.interp(f, f2, s2.imag)
        d = float(np.max(np.abs(a - b)))
        out['dS'] = round(d, 4); ok = d <= t['S_lin']
        out['checks'].append(('S11', 'max |dS11| %.3f (linear)' % d, ok)); out['converged'] &= ok
    if base.get('mcells') and fine.get('mcells'):
        r = fine['mcells'] / base['mcells']
        out['cell_ratio'] = round(r, 2); ok = r >= 1.4
        out['checks'].append(('cells', '%.2f -> %.2f MCells (x%.1f)%s' % (base['mcells'], fine['mcells'], r, '' if ok else ', the finer mesh is not finer: fixed lines dominate'), ok))
        out['converged'] &= ok
    out['verdict'] = ('mesh converged' if out['converged'] else 'NOT converged: refine the mesh') + \
        ' (' + '; '.join('%s %s%s' % (n, m, '' if ok else ' FAIL') for n, m, ok in out['checks']) + ')'
    return out


def delay_ps(f_GHz, s21, lo=None, hi=None):
    """ group delay from the unwrapped phase slope over [lo, hi] GHz """
    f = np.asarray(f_GHz, float); s = np.asarray(s21)
    k = np.ones(len(f), bool)
    if lo is not None: k &= f >= lo
    if hi is not None: k &= f <= hi
    return float(-np.polyfit(2 * np.pi * f[k] * 1e9, np.unwrap(np.angle(s[k])), 1)[0] * 1e12)
