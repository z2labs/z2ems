#!/usr/bin/env python3
"""
mesh_sweep.py: which mesh parameter moves a board's result? (M4b groundwork)

    python mesh_sweep.py board.kicad_pcb --fmax 10e9 --out DIR [--port J1 --port J2] [--variant cuda]

Runs kicad_sim.simulate over a small grid of mesh settings and tabulates Z0, eps_eff, the S21 delay and
|S11| at fmax, plus cells and run time, so the sensitive parameter is visible before the mesher is changed:
  grading ratio (Z2_GRADE_RATIO: off, 2.0, 1.5), lines across a strip/pad (Z2_N_ACROSS: 4, 6, 8),
  resolution (coarse: 1.0, 0.7, 0.5). One parameter varies at a time around the default; the table goes
  to DIR/sweep.json and DIR/sweep.md, a chart to DIR/sweep.png.
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kicad_sim as ks
import convergence

GRID = [('default', {}),
        ('grade_off', {'Z2_GRADE': '0'}), ('grade_2.0', {'Z2_GRADE_RATIO': '2.0'}),
        ('across_6', {'Z2_N_ACROSS': '6'}), ('across_8', {'Z2_N_ACROSS': '8'}),
        ('res_0.7', {'coarse': 0.7}), ('res_0.5', {'coarse': 0.5})]


def one(a, name, env, out):
    saved = {k: os.environ.get(k) for k in ('Z2_GRADE', 'Z2_GRADE_RATIO', 'Z2_N_ACROSS')}
    for k in saved:
        os.environ.pop(k, None)
    coarse = 1.0
    for k, v in env.items():
        if k == 'coarse':
            coarse = v
        else:
            os.environ[k] = v
    t0 = time.time()
    try:
        r = ks.simulate(a.board, a.port, a.fmin, a.fmax, os.path.join(out, name), a.fab, a.variant, a.threads,
                        coarse=coarse, log=lambda m: None)
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    f = np.array(r['f_GHz']) * 1e9
    s11 = np.array(r['S_re']['S11']) + 1j * np.array(r['S_im']['S11'])
    row = dict(name=name, cells=r['mcells'], wall_s=round(time.time() - t0, 1), Z0=r.get('Z0_ohm'), eeff=r.get('eeff'),
               s11_dB_fmax=float(20 * np.log10(abs(s11[-1]) + 1e-15)), s11_dB_max=float(20 * np.log10(np.abs(s11).max() + 1e-15)),
               tau_ps=(round(convergence._tau(r), 2) if 'S21' in r['S_re'] else None),
               quality=r.get('mesh_quality', {}), faults=sum(len(x.get('faults', [])) for x in r['runs']))
    return row


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('board'); ap.add_argument('--port', action='append'); ap.add_argument('--fmin', type=float, default=0.1e9)
    ap.add_argument('--fmax', type=float, default=6e9); ap.add_argument('--fab'); ap.add_argument('--out', required=True)
    ap.add_argument('--variant', default=ks.DEFAULT_VARIANT); ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--only', action='append', help='run only these grid points')
    a = ap.parse_args(argv)
    if not a.port:
        import rf_ports
        ks._setup_python()
        from CSXCAD.KiCadImport import KiCadBoard
        b = KiCadBoard(a.board)
        a.port = [c['ref'] for c in rf_ports.candidates(b, rf_ports.load_project(a.board)) if c['confident']]
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for name, env in GRID:
        if a.only and name not in a.only:
            continue
        row = one(a, name, env, a.out)
        rows.append(row)
        print('%-10s %6.3f MCells %6.0f s  Z0 %6.2f  eeff %6.4f  tau %s ps  |S11|max %6.1f dB  grading %s%s' % (
            name, row['cells'], row['wall_s'], row['Z0'] or 0, row['eeff'] or 0, row['tau_ps'], row['s11_dB_max'],
            row['quality'].get('grading'), '  GPU faults %d' % row['faults'] if row['faults'] else ''), flush=True)
        json.dump(rows, open(os.path.join(a.out, 'sweep.json'), 'w'), indent=1)
    with open(os.path.join(a.out, 'sweep.md'), 'w') as fh:
        fh.write('| setting | MCells | s | Z0 (ohm) | eps_eff | tau (ps) | max S11 (dB) | grading |\n|---|---|---|---|---|---|---|---|\n')
        for r in rows:
            fh.write('| %s | %.3f | %.0f | %.2f | %.4f | %s | %.1f | %s |\n' % (r['name'], r['cells'], r['wall_s'], r['Z0'] or 0, r['eeff'] or 0,
                                                                          r['tau_ps'], r['s11_dB_max'], r['quality'].get('grading')))
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
        names = [r['name'] for r in rows]
        for k, key, lab in ((0, 'Z0', 'Z0 (ohm)'), (1, 'eeff', 'eps_eff'), (2, 'tau_ps', 'delay (ps)')):
            ax[k].bar(names, [r[key] or 0 for r in rows], color='#4c78a8'); ax[k].set_ylabel(lab); ax[k].tick_params(axis='x', rotation=45)
            v = [r[key] for r in rows if r[key]]
            if v:
                ax[k].set_ylim(min(v) - 0.1 * (max(v) - min(v) + 1e-9) - 0.02 * abs(min(v)), max(v) + 0.1 * (max(v) - min(v) + 1e-9) + 0.02 * abs(max(v)))
        fig.suptitle('mesh sweep: %s' % os.path.basename(a.board)); fig.tight_layout()
        fig.savefig(os.path.join(a.out, 'sweep.png'), dpi=120)
    except Exception as e:
        print('no chart:', e)


if __name__ == '__main__':
    main(sys.argv[1:])
