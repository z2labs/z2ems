"""TH vertical SMA jack (5-pin, 6.35 mm square, legs on 5.08 mm) on the BOTTOM side of WES0115, as 3D CSX parts,
optionally the straight adapter + NanoVNA case below it (board perpendicular to the adapter axis).
Coax: pin d 1.27, PTFE bore d 4.1 (eps 2.08) -> ~49 ohm; port = gap between pin end and the closed barrel bottom.
env: SMA_XY centre (sim frame), SMA_H coax length below the board (9.5 = jack only; ~25 = jack + adapter + VNA jack),
     VNA_BOX="wx,ty,len,dx,dy": PEC case below the coax end (wx along x, ty along y, port offset dx,dy from case centre)."""
import os, numpy as np
from CSXCAD.SmoothMeshLines import SmoothMeshLines
CX, CY = [float(v) for v in os.environ.get('SMA_XY', '24.996178,-25.004406').split(',')]
ZB = float(os.environ.get('SMA_ZB', '0.04429'))     # B.Cu sheet
ZT = float(os.environ.get('SMA_ZT', '1.59429'))     # F.Cu sheet
H = float(os.environ.get('SMA_H', '9.5'))
RP, RD, RB, HB = 0.635, 2.05, 3.175, 1.5
Z_BASE, Z_END = ZB - HB, ZB - H
Z_DIEL, Z_PIN = Z_END + 0.4, Z_END + 0.8           # 0.4 mm closed barrel bottom, 0.4 mm port gap
VB = os.environ.get('VNA_BOX')


def merge(old, new, res, mind=0.03):
    """add lines that are not closer than mind to an existing one, then smooth"""
    old = np.unique(np.round(np.asarray(old, float), 6)); add = []
    for v in sorted(new):
        if np.min(np.abs(old - v)) >= mind and (not add or abs(add[-1] - v) >= mind):
            add.append(v)
    L = np.unique(np.concatenate([old, add]))
    L = SmoothMeshLines(L, res, 1.4)
    L = np.unique(np.round(L, 6)); keep = [L[0]]
    for v in L[1:]:
        if v - keep[-1] >= 0.01:
            keep.append(v)
    return list(keep)


def extra(CSX, X, Y, Z, res):
    pec = CSX.AddMetal('sma_body')
    pec.AddBox([CX - 3.175, CY - 3.175, Z_BASE], [CX + 3.175, CY + 3.175, ZB], priority=40)
    pec.AddCylinder([CX, CY, Z_END], [CX, CY, Z_BASE], RB, priority=40)
    d = CSX.AddMaterial('sma_ptfe', epsilon=2.08)
    d.AddCylinder([CX, CY, Z_DIEL], [CX, CY, ZB], RD, priority=50)
    pin = CSX.AddMetal('sma_pin')
    pin.AddCylinder([CX, CY, Z_PIN], [CX, CY, ZT], RP, priority=60)
    u = list(np.arange(-3.6, 3.6001, 0.2))
    nx = [CX + v for v in u]; ny = [CY + v for v in u]
    nz = [Z_END - 0.3, Z_END, Z_DIEL, 0.5 * (Z_DIEL + Z_PIN), Z_PIN, Z_BASE] + list(np.arange(Z_PIN, ZB - 0.25, 0.5))
    if VB:
        wx, ty, ln, dx, dy = [float(v) for v in VB.split(',')]
        x0, x1 = CX - dx - wx / 2, CX - dx + wx / 2; y0, y1 = CY - dy - ty / 2, CY - dy + ty / 2
        z1, z0 = Z_END, Z_END - ln
        CSX.AddMetal('vna_case').AddBox([x0, y0, z0], [x1, y1, z1], priority=30)
        air = float(os.environ.get('VNA_AIR', '60'))
        nx += [x0 - 0.5, x0 + 0.5, x1 - 0.5, x1 + 0.5, x0 - air, x1 + air]
        ny += [y0 - 0.5, y0 + 0.5, y1 - 0.5, y1 + 0.5, y0 - air, y1 + air]
        nz += [z1 - 0.5, z0 - 0.5, z0 + 0.5, z0 - air]
        print('  extra: NanoVNA case %.0f x %.0f x %.0f mm, x %.1f..%.1f y %.1f..%.1f z %.1f..%.1f' % (wx, ty, ln, x0, x1, y0, y1, z0, z1))
    X = merge(X, nx, res); Y = merge(Y, ny, res); Z = merge(Z, nz, res)
    print('  extra: TH SMA-V at (%.3f, %.3f), coax z %.2f..%.2f, port gap z %.2f..%.2f, min cell x %.3f y %.3f z %.3f' % (
        CX, CY, Z_END, ZB, Z_DIEL, Z_PIN, np.min(np.diff(X)), np.min(np.diff(Y)), np.min(np.diff(Z))))
    return X, Y, Z
