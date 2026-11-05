# Silicon Labs MSC-AMS915 antenna kit: measurement and simulation tools

The Silicon Labs MSC-AMS915 "915 MHz antenna matrix" kit holds nine SMA antenna boards (printed ILA/IFA, meander,
BIFA, chip and wire antennas) with measurement reports (application notes AN847, AN848, AN853). It is a good
measured-truth corpus for sub-GHz antenna simulation.

**The Silicon Labs files are not included here** (design files, Gerbers, PADS layout, application notes). Download
the kit documentation from silabs.com and put it into one folder; the scripts read it from `AMS915_DIR`
(default: the current directory).

## Contents

| folder | content |
|---|---|
| `vna/nanovna_v2.py`, `vna/cal.py` | NanoVNA V2 (V2 Plus4) over USB serial without third-party GUI tools: sweep to Touchstone, SOL calibration |
| `meas/` | own raw NanoVNA measurements (WES0115 panic-button IFA through an SMA adapter, uncalibrated; open standard) |
| `kicad/` | PADS ASCII → KiCad 10 pipeline: `fix_panel.py` (PADS copper polygons imported on In1.Cu → F.Cu), `add_round_copper.py` (round-board GND pours and arc antenna pieces the importer drops), `split_panel.py` (panel → 9 boards), `set_matching.py` (matching-network values from AN848), `refill_zones.py`, `net_antenna_copper.py`, `fix_ifa_short.py`, `conn_check.py` (feed-to-antenna copper connectivity), `gerber2geom.py` / `geom2kicad.py` (Gerber reference path), `export_nopanel.py` |
| `tools/gerb_raster.py` | rasterises Gerbers of the original and of the KiCad result for a visual diff |
| `sim/` | Z2 EMS S11 results (Touchstone) of the antenna boards, 0.1–3 GHz |
| `sma_v.py` | `Z2_EXTRA_PY` example: a through-hole vertical SMA body, coax chain and the NanoVNA metal case as extra 3D parts |

Run the KiCad scripts with the KiCad 10 python (they use `pcbnew`); `conn_check.py`, `geom2kicad.py`,
`gerber2geom.py` and `fix_ifa_short.py` also need `shapely`.

## Finding (WES0115 panic-button IFA)

Simulated alone, the board resonates near 985 MHz; the NanoVNA measurement shows 918 MHz. Adding the TH SMA, the
adapter and the NanoVNA's metal case as extra 3D parts (`sma_v.py`) moves the simulated resonance to 926 MHz: the
instrument case acts as the ground counterpoise of this small antenna. The kit's own reports were measured on a
motherboard with a hand on it, which is part of the device under test in the same way.
