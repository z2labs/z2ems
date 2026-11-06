# twin: KiCad board to S-parameters

`kicad_sim.py` runs openEMS on a KiCad board in one step: board and stackup from
`CSXCAD.KiCadImport` (optionally a fab stackup from `CSXCAD.FabStackups`),
lumped 50 Ohm ports at pads, automatic mesh, one run per excited port, Touchstone
output and, for a two-port line, Z0 and effective permittivity.

    python kicad_sim.py board.kicad_pcb --port J1 --port J2 --fmax 6e9
    python kicad_sim.py board.kicad_pcb --port J1 --port J2 --fab JLCPCB:JLC04161H-7628 --copper solid

| option | meaning |
|---|---|
| `--port SPEC` | a port, repeat per port: `REF` (that footprint's signal pad, vertical port to the nearest reference copper below), `REF:PAD/GPAD,GPAD` (gap port from those ground pads to the pad, e.g. `J1:1/2,3,4,5`, `U3:7/6,8`; `GREF:GPAD` for another footprint's pad), `REF:PAD@LAYER` (vertical port to that layer, e.g. `ANT1:1@In1.Cu`) |
| `--fmin/--fmax` | frequency range in Hz |
| `--copper sheet` | zero-thickness copper on the dielectric face (default, fast) |
| `--copper solid` | copper with its thickness (slower: thin cells) |
| `--mask` | include the solder mask |
| `--variant` | `cpu` (default), `cuda` (GPU engine, needs an NVIDIA GPU and a CUDA build) or `auto` (GPU when available) |
| `--components` | two-pad R/L/C footprints as lumped elements (the KiCad plugin turns this on by default) |
| `--converge` | repeat the run on a finer mesh and report the change |

Output in `<board>_sim/` (or `--out`): `<board>.s2p`, `result.json`, one folder
per excitation with the openEMS model and probes.

## How the model is built

- **Ports**: a port is a + pad and a - reference.
  - *Gap port* (- = ground pads on the same layer: an edge connector's GND pads,
    an IC's neighbouring GND pins): one lumped sub-port in the copper plane per
    ground pad, across the real gap (from the + pad edge to the nearest copper
    of the ground net in that direction, so a ground pour between the pads
    counts), each with n·Z0, driven together; the port voltage is their mean,
    the port current their sum. On the BWSMA edge-launch coupon the gap ports
    give |S21| -0.41 / -1.88 dB at 1 / 6 GHz against -0.40 / -1.77 dB with
    vertical ports.
  - *Vertical port* (- = a plane, or the nearest reference copper): at the pad
    edge that faces the nearest board edge, inset 0.1 mm into the pad, from
    the sheet of the reference layer below up to the pad sheet.
- **Copper as sheets** lies on the dielectric face of each layer
  (`KiCadBoard.SheetZ`): the trace-to-plane distance of an outer microstrip is
  exactly the prepreg thickness. Dielectrics reach from sheet to sheet, so there
  are no air gaps at inner layers.
- **Mesh**: λ/20 in the densest dielectric; thirds rule on the edges of tracks
  and rectangular pads of the port nets (no line on a metal edge: such a line
  widens the strip by half a cell, -8 % Z0 on the coupon); every dielectric at
  least four cells; a line through every plated hole of the reference net.
  Lines and primitive coordinates are rounded to 1e-5 mm.
- **Time**: energy end criterion 1e-4, but at most 25 board transit times:
  lossless plane-via loops keep a static field from the DC part of the
  excitation, which never decays.

## Engine

`kicad_sim.py` finds openEMS in this order: `Z2_OPENEMS` (path of the executable), `bin/openEMS[.exe]` of an
engine package next to `twin/`, `openEMS` on PATH. The GPU engine is the same executable with `--engine=cuda`
(a `-DWITH_CUDA=ON` build); `Z2_OPENEMS_CUDA` points to a separate GPU build. The CUDA engine checks its own
memory (`--gpu-guard`, a checksum of the coefficients at every energy check) and a run that dies of a GPU fault is
repeated (`Z2_GPU_RETRIES`, default 2), so old and worn GPUs give a clear error instead of a wrong result.

Useful environment switches: `Z2_OPENEMS_ARGS` (extra engine options; one given here replaces the default),
`Z2_CONFORMAL=0` (no conformal PEC), `Z2_GRADE=0` / `Z2_GRADE_RATIO` (mesh grading), `Z2_AIR_MARGIN` (air around
the board in mm, antennas need about a quarter wavelength), `Z2_REFINE_BOX="x0,y0,x1,y1,r;..."` (uniform mesh
lines every r mm in a box), `Z2_PORT_TRACK_EXIT=1` (port across the track that leaves a THT pad),
`Z2_PORT_FORCE=<ports.json>` (reuse stored port definitions), `Z2_EXTRA_PY=<file>` (add 3D parts the board file
does not hold, e.g. a through-hole SMA body or an instrument case; see `examples/silabs-ams915/sma_v.py`),
`Z2_LOSS_F` (frequency at which the dielectric loss tangent is matched).

## Fields

Runs from the viewer or the KiCad plugin also store the fields (`fields.py`), for every driven port, at three
frequencies (0.25, 0.5 and 0.9 x fmax). Everything is scaled to **1 W incident at the driven port** (peak
phasors, phase 0 = the incident wave at the port), so the numbers are physical:

- the **surface current density Js** (A/m) on every copper layer, from a rot(H) dump on the sheet times the dual
  cell height: the signal current and its return current in the planes;
- the **E-field** in a box around the whole stack at reduced resolution, so a cut plane can be chosen after the run;
- the **voltage** of each copper layer against the next one down, V = -integral of Ez between the two sheets.

Check on the 50 ohm coupon (1 W, 1.5 GHz): 0.183 A on the strip and 0.195 A in the In1 return path for 0.2 A
expected, 9.8-10.1 V for 10 V expected. The viewer shows them as colour maps, arrows and current / E-field lines,
animated in phase; the KiCad plugin draws the maps and the lines on the board. The files are float16 arrays next to
`result.json` (`fields/`, a few MB per port and frequency); `result.json["fields"]` indexes them (version 2).

## Validation

| example | result with the default settings | check |
|---|---|---|
| `examples/coupon_msl50.kicad_pcb` (50 Ω microstrip, 40 mm, JLC04161H-7628) | Z0 49.50 Ω, εeff 3.353, \|S11\| < −28 dB 0.1–6 GHz (CPU and CUDA engine agree to 0.005 Ω) | Hammerstad–Jensen with 35 µm copper: 49.35 Ω / 3.337 |
| `examples/series_c1p.kicad_pcb` (1 pF series capacitor on the same line, `--components`) | \|S21\| −5.9 dB at 1 GHz, −2.5 dB at 2 GHz | ideal 1 pF between 50 Ω ports: −5.5 / −2.1 dB (the rest: pad and gap parasitics) |

`make_coupon.py` writes the microstrip coupon: 0.3726 mm track on F.Cu over In1 (0.2104 mm 7628 prepreg,
εr 4.4), GND planes on In1 and B.Cu stitched with vias, SMD launch pads J1/J2.

## Edge-SMA launch (`launch_study.py`, `bwsma_model.py`)

`launch_study.py` simulates an edge-mounted SMA launch with a connector model (centre pin on the pad, 50 Ω PTFE coax
stub, metal flange, lumped port at the coax end). `bwsma_model.py` is the model of the low-cost BAT Wireless
BWSMA-KE-P001 edge SMA jack (LCSC C496550); `models/BWSMA-KE-P001.step/.wrl` are its 3D models
(regenerating them needs `cadquery`; the simulation itself does not).
`examples/bwsma_air_launch_coupon.kicad_pcb` is a 50 × 12 mm GCPW coupon with an "air" launch for that connector: the
board is routed away under the first 2.8 mm of the centre pin, and a 1.55 mm pad at the pin tip compensates the
series inductance. With the default mesh it gives max |S11| −19.5 dB up to 3 GHz, −17.8 dB up to 6 GHz and
|S21| better than −0.8 dB up to 6 GHz (simulated). The naive footprint for the same connector (coupon 29 of
`make_panel.py`) gives max |S11| −0.9 dB and |S21| down to −8.9 dB in the same band. VNA measurements of the
manufactured coupons will follow. Simulate with `python make_panel.py out --sim 1,29`.

## Tools

| script | purpose |
|---|---|
| `mesh_drc.py` | "DRC on the mesh": labels the staircase copper islands per layer, reports SHORT (two nets on one island) and OPEN (a signal net split on a layer) before a run |
| `mesh_sweep.py` | one mesh parameter at a time, table and plot |
| `convergence.py` | the `--converge` check |
| `diag_line.py` | straight / 10° / zig-zag / 45° microstrip lines: checks that diagonal copper stays connected in the staircase mesh |
| `plot_pattern.py` | 2D and 3D radiation-pattern plots of an antenna run |
| `z2mon.py`, `runviz.py` | run monitor in the browser: geometry, progress, S-parameters, currents, far field |
| `crop.py` | cut a region of a large board into a small simulation model |
