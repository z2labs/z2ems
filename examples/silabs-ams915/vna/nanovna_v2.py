#!/usr/bin/env python3
"""Minimal NanoVNA V2 / V2 Plus4 driver (NanoRFE binary protocol) -> raw S11/S21 and Touchstone.

Protocol (USB CDC, baud ignored): 1-byte opcodes
  0x00 NOP  0x0d INDICATE (-> '2')  0x10/0x11/0x12 READ 1/2/4 bytes  0x18 READFIFO  0x20/21/22/23 WRITE 1/2/4/8 bytes
Registers: 0x00 sweepStartHz(u64) 0x10 sweepStepHz(u64) 0x20 sweepPoints(u16) 0x22 valuesPerFrequency(u16)
           0x30 valuesFIFO (32-byte records) 0xf0 deviceVariant 0xf1 protocolVersion 0xf2 hwRev 0xf3/0xf4 fw major/minor
FIFO record: fwd0 re/im, rev0 re/im, rev1 re/im (int32 each), freqIndex u16, 10 reserved.  S11 = rev0/fwd0, S21 = rev1/fwd0.
The V2 applies NO calibration over USB: these are raw receiver ratios; SOL(T) error-term correction is done in software (cal.py).

  python nanovna_v2.py --port COM6 info
  python nanovna_v2.py --port COM6 sweep --start 100e6 --stop 3e9 --points 1001 --avg 4 --out dut_raw.s2p
"""
import argparse, json, struct, sys, time
import numpy as np
import serial

class NanoVNAV2:
    REC = struct.Struct('<6iH6s')  # 32 bytes: 6 x int32 + u16 index + 6 reserved
    def __init__(self, port, timeout=2.0):
        self.s = serial.Serial(port, 115200, timeout=timeout)
        self.s.write(b'\x00' * 8 + b'\x0d'); r = self.s.read(1)
        if r != b'2':
            raise IOError('no NanoVNA V2 on %s (indicate -> %r)' % (port, r))
    def _rd(self, reg, n):
        self.s.write(bytes([{1: 0x10, 2: 0x11, 4: 0x12}[n], reg])); b = self.s.read(n)
        if len(b) != n: raise IOError('short read reg 0x%02x' % reg)
        return int.from_bytes(b, 'little')
    def _wr(self, reg, val, n):
        self.s.write(bytes([{1: 0x20, 2: 0x21, 4: 0x22, 8: 0x23}[n], reg]) + int(val).to_bytes(n, 'little'))
    def info(self):
        return dict(variant=self._rd(0xf0, 1), protocol=self._rd(0xf1, 1), hw_rev=self._rd(0xf2, 1),
                    fw='%d.%d' % (self._rd(0xf3, 1), self._rd(0xf4, 1)))
    def set_sweep(self, start, stop, points, avg=1):
        step = (stop - start) / (points - 1)
        self._wr(0x00, round(start), 8); self._wr(0x10, round(step), 8)
        self._wr(0x20, points, 2); self._wr(0x22, avg, 2)
        self.freqs = start + step * np.arange(points); self.points = points; self.avg = avg
    def sweep(self):
        """One full sweep: returns (freqs, S11, S21) averaged over valuesPerFrequency."""
        self._wr(0x30, 0, 1)                      # clear FIFO (any write to 0x30)
        n = self.points * self.avg
        acc11 = np.zeros(self.points, complex); acc21 = np.zeros(self.points, complex); cnt = np.zeros(self.points)
        got = 0
        while got < n:
            k = min(255, n - got)
            self.s.write(bytes([0x18, 0x30, k])); b = self.s.read(32 * k)
            if len(b) != 32 * k: raise IOError('FIFO short read (%d/%d)' % (len(b), 32 * k))
            for i in range(k):
                f0r, f0i, r0r, r0i, r1r, r1i, idx, _ = self.REC.unpack_from(b, 32 * i)
                f0 = complex(f0r, f0i)
                if f0 == 0: continue
                acc11[idx] += complex(r0r, r0i) / f0; acc21[idx] += complex(r1r, r1i) / f0; cnt[idx] += 1
            got += k
        cnt[cnt == 0] = 1
        return self.freqs, acc11 / cnt, acc21 / cnt
    def close(self): self.s.close()

def write_touchstone(path, f, s11, s21=None, comment=''):
    with open(path, 'w') as fh:
        fh.write('! %s\n# Hz S RI R 50\n' % comment)
        for i in range(len(f)):
            if s21 is None: fh.write('%.0f %.9f %.9f\n' % (f[i], s11[i].real, s11[i].imag))
            else: fh.write('%.0f %.9f %.9f %.9f %.9f 0 0 0 0\n' % (f[i], s11[i].real, s11[i].imag, s21[i].real, s21[i].imag))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--port', required=True)
    sub = ap.add_subparsers(dest='cmd'); sub.add_parser('info')
    sw = sub.add_parser('sweep'); sw.add_argument('--start', type=float, default=100e6); sw.add_argument('--stop', type=float, default=3e9)
    sw.add_argument('--points', type=int, default=1001); sw.add_argument('--avg', type=int, default=4); sw.add_argument('--repeat', type=int, default=1)
    sw.add_argument('--out', required=True); sw.add_argument('--note', default='')
    a = ap.parse_args()
    v = NanoVNAV2(a.port)
    if a.cmd == 'info':
        print(json.dumps(v.info())); return
    v.set_sweep(a.start, a.stop, a.points, a.avg)
    runs = []
    for r in range(a.repeat):
        t = time.time(); f, s11, s21 = v.sweep(); runs.append((s11, s21)); print('sweep %d: %.1f s' % (r, time.time() - t))
    s11 = np.median([x[0] for x in runs], axis=0) if a.repeat > 1 else runs[0][0]
    s21 = np.median([x[1] for x in runs], axis=0) if a.repeat > 1 else runs[0][1]
    info = v.info(); v.close()
    write_touchstone(a.out, f, s11, s21, 'NanoVNA V2 raw (uncalibrated) %s fw %s %s avg=%d repeat=%d %s' % (
        json.dumps(info), info['fw'], time.strftime('%Y-%m-%d %H:%M'), a.avg, a.repeat, a.note))
    print('wrote', a.out, '|S11| dB at 915 MHz ~', round(20 * np.log10(abs(np.interp(915e6, f, abs(s11)))), 1))
if __name__ == '__main__': main()
