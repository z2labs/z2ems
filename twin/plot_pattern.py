#!/usr/bin/env python3
""" Radiation-pattern plots for an antenna run, in the usual industry form.

    Input:    nf2ff_<f>GHz.h5 written by the antenna evaluation (NF2FF box, theta 0..180 step 2, phi 0..355
              step 5, complex E_theta / E_phi) + result.json (radiation efficiency -> gain)

    Outputs <out>_pattern2d.png (E- and H-plane cuts, polar + rectangular, co-/cross-pol in Ludwig-3),
    <out>_pattern3d.png (3D directivity surface) and <out>_pattern.json (peak, HPBW, SLL, F/B, cross-pol per cut).

        python plot_pattern.py <case_dir> [--out <prefix>] [--f 24.5] [--no3d]
        <case_dir> holds result.json and sim/nf2ff_*.h5 (or is the sim dir itself)
"""
import glob, json, math, os, re, sys

import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import mpl_toolkits
    _mt = os.path.join(os.path.dirname(matplotlib.__file__), '..', 'mpl_toolkits')   # prefer the toolkit shipped with this matplotlib
    if os.path.isdir(_mt) and os.path.abspath(_mt) not in mpl_toolkits.__path__:
        mpl_toolkits.__path__.insert(0, os.path.abspath(_mt))
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    if '3d' not in matplotlib.projections.get_projection_names():
        matplotlib.projections.register_projection(Axes3D)
except ImportError:                      # pragma: no cover
    plt = None

C_Z2, C_Z2X, C_REF, C_REFX = '#1f5fa8', '#8fb3dc', '#c8531d', '#e8a985'   # Z2 EM co/cross, reference co/cross


# ---------------------------------------------------------------- loading
def load_openems(h5path):
    """ theta(deg, n), phi(deg, m), E_theta, E_phi complex (n, m), Dmax (linear) """
    import h5py
    with h5py.File(h5path, 'r') as f:
        th = np.degrees(f['Mesh/theta'][:])
        ph = np.degrees(f['Mesh/phi'][:])
        g = f['nf2ff']
        Et = g['E_theta/FD/f0'][:]
        Ep = g['E_phi/FD/f0'][:]
        Dmax = float(np.ravel(g.attrs['Dmax'])[0])
        fq = float(np.ravel(g.attrs['Frequency'])[0])
    return dict(theta=th, phi=ph, Et=Et, Ep=Ep, Dmax=Dmax, f_GHz=fq / 1e9, src='Z2 EM')



def directivity(p):
    """ Directivity (dBi) of the total field and of the Ludwig-3 co-/cross-polar components on the grid.
        The main polarisation (x or y) is taken from the field at the pattern maximum. """
    Et, Ep, th, ph = p['Et'], p['Ep'], p['theta'], p['phi']
    P = np.abs(Et) ** 2 + np.abs(Ep) ** 2
    scale = p['Dmax'] / np.max(P)                     # |E|^2 * scale = D (linear)
    PH = np.radians(ph)[None, :]
    Ex = Et * np.cos(PH) - Ep * np.sin(PH)            # Ludwig-3, x-polarised reference
    Ey = Et * np.sin(PH) + Ep * np.cos(PH)
    i = np.unravel_index(np.argmax(P), P.shape)
    xpol = abs(Ex[i]) >= abs(Ey[i])
    co, cx = (Ex, Ey) if xpol else (Ey, Ex)
    db = lambda a: 10 * np.log10(np.maximum(a * scale, 1e-12))
    if p.get('circular'):
        D = 10 * np.log10(np.maximum(10 ** (p['D_file'] / 10), 1e-12))
        return dict(D=D, co=D, cx=np.full_like(D, -120.0), pol=None, total_only=True, theta=th, phi=ph,
                    imax=np.unravel_index(np.argmax(D), D.shape))
    return dict(D=db(P), co=db(np.abs(co) ** 2), cx=db(np.abs(cx) ** 2), pol='x' if xpol else 'y',
                theta=th, phi=ph, imax=i)


def cut(D, theta, phi, phi_deg):
    """ Full 360-degree cut through the z axis in the plane phi = phi_deg: returns angle (-180..180, positive
        towards phi_deg) and values; built from the theta 0..180 half-planes at phi_deg and phi_deg + 180. """
    j1 = int(np.argmin(np.abs((phi - phi_deg) % 360)))
    j2 = int(np.argmin(np.abs((phi - phi_deg - 180) % 360)))
    ang = np.concatenate([-theta[::-1][:-1], theta])
    val = np.concatenate([D[::-1, j2][:-1], D[:, j1]])
    return ang, val


def cut_metrics(ang, co, cx=None):
    i = int(np.argmax(co)); pk = float(co[i])
    # HPBW: -3 dB points on both sides of the peak (linear interpolation)
    def edge(step):
        k = i
        while 0 <= k + step < len(co) and co[k + step] > pk - 3:
            k += step
        if not (0 <= k + step < len(co)):
            return None
        a0, a1, v0, v1 = ang[k], ang[k + step], co[k], co[k + step]
        return float(a0 + (a1 - a0) * (v0 - (pk - 3)) / (v0 - v1)) if v0 != v1 else float(a0)
    lo, hi = edge(-1), edge(1)
    hpbw = (hi - lo) if lo is not None and hi is not None else None
    # first side lobe: highest local maximum outside the main lobe (beyond the first minimum on each side)
    sll = None
    for step in (-1, 1):
        k = i
        while 0 <= k + step < len(co):            # walk to the first minimum at least 3 dB below the peak (ripple on a flat top is not a null)
            while 0 <= k + step < len(co) and co[k + step] <= co[k]:
                k += step
            if co[k] <= pk - 3 or not (0 <= k + step < len(co)):
                break
            while 0 <= k + step < len(co) and co[k + step] >= co[k]:
                k += step
        seg = co[:k] if step < 0 else co[k + 1:]
        if len(seg):
            s = float(np.max(seg)) - pk
            sll = s if sll is None else max(sll, s)
    j = int(np.argmin(np.abs(((ang - ang[i] + 180 + 180) % 360) - 180)))   # opposite direction
    fb = pk - float(co[j])
    out = dict(peak_dBi=round(pk, 2), peak_deg=float(ang[i]), hpbw_deg=round(hpbw, 1) if hpbw is not None else None,
               sll_dB=round(sll, 1) if sll is not None else None, fb_dB=round(fb, 1))
    if cx is not None:
        out['xpol_dB'] = round(float(cx[i] - co[i]), 1)      # cross-polar level at the co-polar peak
    return out


def eff_db(p, res=None):
    """ 10 log10 of the radiation efficiency (gain = directivity + this) """
    e = p.get('eff_rad')
    if e is None and res:
        for ff in res.get('farfield', []):
            if abs(ff['f_GHz'] - p['f_GHz']) < 1e-3 and ff.get('eff_rad'):
                e = ff['eff_rad']
    return 10 * math.log10(e) if e else None


# ---------------------------------------------------------------- plots
def plot2d(z2, ref, out, title):
    pz = directivity(z2)
    pr = directivity(ref) if ref else None
    pol = pz['pol']
    planes = [('E-plane', 0 if pol == 'x' else 90), ('H-plane', 90 if pol == 'x' else 0)]
    top = max(pz['D'].max(), pr['D'].max() if pr else -99)
    vmax = 5 * math.ceil(top / 5) + 5
    vmin = vmax - 45
    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1])
    metrics = {}
    reflab = 'reference total' if pr and pr.get('total_only') else 'reference co-pol'
    for k, (name, phd) in enumerate(planes):
        ang, co = cut(pz['co'], pz['theta'], pz['phi'], phd)
        _, cx = cut(pz['cx'], pz['theta'], pz['phi'], phd)
        _, dt = cut(pz['D'], pz['theta'], pz['phi'], phd)
        m = dict(z2=cut_metrics(ang, co, cx))
        if pr:
            angr, cor = cut(pr['co'], pr['theta'], pr['phi'], phd)
            _, cxr = cut(pr['cx'], pr['theta'], pr['phi'], phd)
            m['ref'] = cut_metrics(angr, cor, None if pr.get('total_only') else cxr)
            if pr.get('total_only'):
                co_cmp = dt          # compare total against total
            else:
                co_cmp = co
            # difference statistics on the common grid inside the main lobe (ref co within 10 dB of its peak)
            if len(angr) == len(ang) and np.allclose(angr, ang):
                sel = cor > cor.max() - 10
                m['delta'] = dict(peak_dB=round(m['z2']['peak_dBi'] - m['ref']['peak_dBi'], 2),
                                  rms_mainlobe_dB=round(float(np.sqrt(np.mean((co_cmp[sel] - cor[sel]) ** 2))), 2),
                                  max_mainlobe_dB=round(float(np.max(np.abs(co_cmp[sel] - cor[sel]))), 2),
                                  hpbw_deg=round(m['z2']['hpbw_deg'] - m['ref']['hpbw_deg'], 1)
                                  if m['z2']['hpbw_deg'] is not None and m['ref']['hpbw_deg'] is not None else None)
        metrics[name] = m
        # polar
        ax = fig.add_subplot(gs[0, k], projection='polar')
        ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
        r = np.radians(ang)
        ax.plot(r, np.maximum(co, vmin), color=C_Z2, lw=2, label='Z2 EM co-pol')
        ax.plot(r, np.maximum(cx, vmin), color=C_Z2X, lw=1.2, label='Z2 EM cross-pol')
        if pr:
            ax.plot(np.radians(angr), np.maximum(cor, vmin), color=C_REF, lw=2, ls='--', label=reflab)
            if not pr.get('total_only'):
                ax.plot(np.radians(angr), np.maximum(cxr, vmin), color=C_REFX, lw=1.2, ls='--', label='reference cross-pol')
        ax.set_rlim(vmin, vmax)
        ax.set_rticks(np.arange(vmin + 5, vmax + 1, 10))
        ax.set_rlabel_position(112)
        ax.set_thetagrids(np.arange(0, 360, 30), labels=['%d°' % (a if a <= 180 else a - 360) for a in range(0, 360, 30)])
        ax.grid(color='#cccccc', lw=0.6)
        ax.set_title('%s (φ = %d°), directivity [dBi]' % (name, phd), pad=14, fontsize=11)
        # rectangular
        ax2 = fig.add_subplot(gs[1, k])
        ax2.plot(ang, co, color=C_Z2, lw=2, label='Z2 EM co-pol')
        ax2.plot(ang, cx, color=C_Z2X, lw=1.2, label='Z2 EM cross-pol')
        ax2.plot(ang, dt, color=C_Z2, lw=0.8, ls=':', label='Z2 EM total')
        if pr:
            ax2.plot(angr, cor, color=C_REF, lw=2, ls='--', label=reflab)
            if not pr.get('total_only'):
                ax2.plot(angr, cxr, color=C_REFX, lw=1.2, ls='--', label='reference cross-pol')
        ax2.set_xlim(-180, 180); ax2.set_ylim(vmin, vmax)
        ax2.set_xticks(np.arange(-180, 181, 30))
        ax2.set_xlabel('θ [°]  (0° = boresight, +z)'); ax2.set_ylabel('directivity [dBi]')
        ax2.grid(color='#dddddd', lw=0.6)
        z = m['z2']
        txt = 'Z2 EM: peak %.2f dBi @ %.0f°, HPBW %s°, SLL %s dB, F/B %.1f dB, X-pol %s dB' % (
            z['peak_dBi'], z['peak_deg'], z['hpbw_deg'], z['sll_dB'], z['fb_dB'], z.get('xpol_dB'))
        if pr:
            rr = m['ref']
            txt += '\nref:   peak %.2f dBi @ %.0f°, HPBW %s°, SLL %s dB, F/B %.1f dB, X-pol %s dB' % (
                rr['peak_dBi'], rr['peak_deg'], rr['hpbw_deg'], rr['sll_dB'], rr['fb_dB'], rr.get('xpol_dB'))
            if 'delta' in m:
                d = m['delta']
                txt += '\nΔ peak %+.2f dB, main-lobe RMS %.2f dB (max %.2f dB), Δ HPBW %s°' % (
                    d['peak_dB'], d['rms_mainlobe_dB'], d['max_mainlobe_dB'], d['hpbw_deg'])
        ax2.text(0.01, 0.02, txt, transform=ax2.transAxes, fontsize=7.5, family='monospace', va='bottom',
                 bbox=dict(boxstyle='round', fc='white', ec='#cccccc', alpha=0.9))
        if k == 0:
            ax2.legend(loc='upper right', fontsize=8, ncol=1)
    g = eff_db(z2, z2.get('_res'))
    gr = eff_db(ref) if ref else None
    sub = 'f = %.3f GHz, main polarisation %s (Ludwig-3), Dmax Z2 EM %.2f dBi' % (z2['f_GHz'], pol, pz['D'].max())
    if g is not None:
        sub += ', η_rad %.1f %% → gain %.2f dBi' % (100 * 10 ** (g / 10), pz['D'].max() + g)
    if pr:
        sub += '\nreference: Dmax %.2f dBi' % pr['D'].max()
        if gr is not None:
            sub += ', η_rad %.1f %% → gain %.2f dBi' % (100 * 10 ** (gr / 10), pr['D'].max() + gr)
    fig.suptitle(title + '\n' + sub, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out + '_pattern2d.png', dpi=130)
    plt.close(fig)
    metrics['summary'] = dict(f_GHz=z2['f_GHz'], pol=pol, Dmax_z2_dBi=round(float(pz['D'].max()), 2),
                              gain_z2_dBi=round(float(pz['D'].max() + g), 2) if g is not None else None,
                              Dmax_ref_dBi=round(float(pr['D'].max()), 2) if pr else None,
                              gain_ref_dBi=round(float(pr['D'].max() + gr), 2) if pr and gr is not None else None)
    if pr and pr['D'].shape == pz['D'].shape:
        sel = pr['D'] > pr['D'].max() - 10
        metrics['summary']['rms_mainlobe_3d_dB'] = round(float(np.sqrt(np.mean((pz['D'][sel] - pr['D'][sel]) ** 2))), 2)
        metrics['summary']['theta_max_z2_deg'] = float(pz['theta'][pz['imax'][0]]); metrics['summary']['phi_max_z2_deg'] = float(pz['phi'][pz['imax'][1]])
        metrics['summary']['theta_max_ref_deg'] = float(pr['theta'][pr['imax'][0]]); metrics['summary']['phi_max_ref_deg'] = float(pr['phi'][pr['imax'][1]])
    return metrics


def _surface(ax, D, theta, phi, floor, vmin, vmax, cmap):
    th = np.radians(np.concatenate([theta])); ph = np.radians(np.concatenate([phi, [phi[0] + 360]]))
    Dd = np.concatenate([D, D[:, :1]], axis=1)
    R = np.maximum(Dd - floor, 0)
    T, P = np.meshgrid(th, ph, indexing='ij')
    X = R * np.sin(T) * np.cos(P); Y = R * np.sin(T) * np.sin(P); Z = R * np.cos(T)
    norm = plt.Normalize(vmin, vmax)
    ax.plot_surface(X, Y, Z, facecolors=cmap(norm(Dd)), rstride=1, cstride=1, linewidth=0, antialiased=True, shade=False)
    L = vmax - floor
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.view_init(elev=25, azim=-55)
    return norm


def plot3d(z2, ref, out, title):
    pz = directivity(z2)
    pr = directivity(ref) if ref else None
    top = max(pz['D'].max(), pr['D'].max() if pr else -99)
    vmax = math.ceil(top); floor = vmax - 30                     # 30 dB dynamic range, the usual 3D display
    n = 3 if pr is not None else 1
    fig = plt.figure(figsize=(5.8 * n + 0.8, 6.0))
    cmap = plt.get_cmap('viridis')
    w = 0.92 / n
    pos = lambda k: [0.01 + k * w, 0.0, w - 0.02, 0.9]
    ax = fig.add_axes(pos(0), projection='3d')
    norm = _surface(ax, pz['D'], pz['theta'], pz['phi'], floor, floor, vmax, cmap)
    ax.set_title('Z2 EM, Dmax %.2f dBi' % pz['D'].max(), fontsize=10, y=0.97)
    if pr is not None:
        ax = fig.add_axes(pos(1), projection='3d')
        _surface(ax, pr['D'], pr['theta'], pr['phi'], floor, floor, vmax, cmap)
        ax.set_title('reference, Dmax %.2f dBi' % pr['D'].max(), fontsize=10, y=0.97)
        ax = fig.add_axes(pos(2), projection='3d')
        if pr['D'].shape == pz['D'].shape:
            diff = pz['D'] - pr['D']
            dn = plt.Normalize(-3, 3)
            cm2 = plt.get_cmap('coolwarm')
            # radius from the reference pattern, colour = difference (dB)
            th = np.radians(pr['theta']); ph = np.radians(np.concatenate([pr['phi'], [pr['phi'][0] + 360]]))
            Dd = np.concatenate([pr['D'], pr['D'][:, :1]], axis=1); dd = np.concatenate([diff, diff[:, :1]], axis=1)
            R = np.maximum(Dd - floor, 0); T, P = np.meshgrid(th, ph, indexing='ij')
            ax.plot_surface(R * np.sin(T) * np.cos(P), R * np.sin(T) * np.sin(P), R * np.cos(T), facecolors=cm2(dn(dd)),
                            rstride=1, cstride=1, linewidth=0, shade=False)
            L = vmax - floor
            ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L); ax.set_box_aspect((1, 1, 1))
            ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); ax.view_init(elev=25, azim=-55)
            sel = pr['D'] > pr['D'].max() - 10
            ax.set_title('Z2 EM − reference [dB], main-lobe RMS %.2f dB' % float(np.sqrt(np.mean(diff[sel] ** 2))), fontsize=10, y=0.97)
            cax = fig.add_axes([0.01 + 3 * w - 0.01, 0.2, 0.012, 0.5])
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=dn, cmap=cm2), cax=cax)
            cb.set_label('Δ directivity [dB]')
        else:
            ax.set_title('grids differ, no difference plot', fontsize=10)
    cax = fig.add_axes([0.01 + (2 if pr is not None else 1) * w - 0.03, 0.2, 0.012, 0.5])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label('directivity [dBi]  (radius: %d dB range)' % (vmax - floor))
    fig.suptitle('%s, f = %.3f GHz, 3D directivity (z = boresight)' % (title, z2['f_GHz']), fontsize=11, y=0.985)
    fig.savefig(out + '_pattern3d.png', dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- driver

def make_plots(case_dir, out=None, f_sel=None, do3d=True, title=None):
    sim = case_dir if glob.glob(os.path.join(case_dir, 'nf2ff_*GHz.h5')) else os.path.join(case_dir, 'sim')
    h5s = sorted(glob.glob(os.path.join(sim, 'nf2ff_*GHz.h5')))
    if not h5s:
        raise SystemExit('no nf2ff_*GHz.h5 in ' + sim)
    res = None
    rj = os.path.join(case_dir, 'result.json')
    if os.path.exists(rj):
        res = json.load(open(rj))
    if f_sel is not None:
        h5s = [min(h5s, key=lambda h: abs(float(re.search(r'ff_([\d.]+)GHz', h).group(1)) - f_sel))]
    out = out or os.path.join(case_dir, os.path.basename(os.path.abspath(case_dir)))
    title = title or os.path.basename(os.path.abspath(case_dir))
    allm = {}
    for h in h5s:
        z2 = load_openems(h); z2['_res'] = res
        ref = None
        o = out + ('_%.2fGHz' % z2['f_GHz'] if len(h5s) > 1 else '')
        m = plot2d(z2, ref, o, title)
        if do3d:
            plot3d(z2, ref, o, title)
        m['files'] = dict(z2=h)
        allm['%.4f' % z2['f_GHz']] = m
        print('%s  f=%.3f GHz  Dmax Z2 %.2f dBi%s' % (o, z2['f_GHz'], m['summary']['Dmax_z2_dBi'],
              '  ref %.2f dBi  3D main-lobe RMS %s dB' % (m['summary']['Dmax_ref_dBi'], m['summary'].get('rms_mainlobe_3d_dB')) if ref else ''))
    json.dump(allm, open(out + '_pattern.json', 'w'), indent=1)
    return allm


def main():
    a = sys.argv[1:]
    if not a:
        raise SystemExit(__doc__)
    g = lambda k, d=None: a[a.index(k) + 1] if k in a else d
    make_plots(a[0], out=g('--out'), f_sel=float(g('--f')) if g('--f') else None, do3d='--no3d' not in a,
               title=g('--title'))


if __name__ == '__main__':
    main()
