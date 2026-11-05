"""Closed-form microstrip reference: Hammerstad-Jensen (1980) with the finite-thickness correction."""
import math


def msl(w, h, er, t=0.0):
    """ Z0 [Ohm] and eeff of a microstrip, quasi-static (Hammerstad & Jensen, IEEE MTT-S 1980). """
    def z0_air(u):
        f = 6 + (2 * math.pi - 6) * math.exp(-(30.666 / u) ** 0.7528)
        return 60 * math.log(f / u + math.sqrt(1 + (2 / u) ** 2))

    def ee(u):
        a = 1 + math.log((u ** 4 + (u / 52) ** 2) / (u ** 4 + 0.432)) / 49 + math.log(1 + (u / 18.1) ** 3) / 18.7
        b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
        return (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)

    u = w / h
    if t > 0:
        T = t / h
        du1 = T / math.pi * math.log(1 + 4 * math.e / (T * (1 / math.tanh(math.sqrt(6.517 * u))) ** 2))
        dur = 0.5 * (1 + 1 / math.cosh(math.sqrt(er - 1))) * du1
        u1, ur = u + du1, u + dur
    else:
        u1 = ur = u
    e = ee(ur)
    z = z0_air(ur) / math.sqrt(e)
    e *= (z0_air(u1) / z0_air(ur)) ** 2
    return z, e


def width_for(z0, h, er, t=0.0):
    lo, hi = 0.01 * h, 20 * h
    for _ in range(80):
        m = 0.5 * (lo + hi)
        if msl(m, h, er, t)[0] > z0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


if __name__ == '__main__':
    h, er, t = 0.2104, 4.4, 0.035
    w = width_for(50, h, er, t)
    print('w50 = %.4f mm' % w, msl(w, h, er, t), 'at w=0.405:', msl(0.405, h, er, t), 'thin:', msl(0.405, h, er))


def gcpw(w, s, h, er):
    """ Conductor-backed CPW (Simons): Z0, eeff for centre width w, gap s, substrate h, zero copper thickness. """
    from scipy.special import ellipk
    a, b = w / 2, w / 2 + s
    k = a / b
    k3 = math.tanh(math.pi * a / (2 * h)) / math.tanh(math.pi * b / (2 * h))
    r = lambda q: ellipk(q ** 2) / ellipk(1 - q ** 2)          # K(k)/K(k')
    eeff = (1 + er * r(k3) / r(k)) / (1 + r(k3) / r(k))
    return 60 * math.pi / math.sqrt(eeff) / (r(k) + r(k3)), eeff


def gcpw_width_for(z0, s, h, er):
    lo, hi = 0.05, 5.0
    for _ in range(60):
        m = (lo + hi) / 2
        if gcpw(m, s, h, er)[0] > z0:
            lo = m
        else:
            hi = m
    return (lo + hi) / 2
