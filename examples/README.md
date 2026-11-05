# Examples

| file | what it is | run |
|---|---|---|
| `coupon_msl50.kicad_pcb` | 50 Ω microstrip, 40 mm, JLCPCB JLC04161H-7628 4-layer stackup, SMD launch pads J1/J2 | `python ../twin/kicad_sim.py coupon_msl50.kicad_pcb --port J1 --port J2 --fmax 6e9` → Z0 49.5 Ω, εeff 3.35 |
| `series_c1p.kicad_pcb` | the same line with a 1 pF series capacitor (C1) | add `--components` → \|S21\| −5.9 dB at 1 GHz (ideal 1 pF: −5.5 dB) |
| `dbg.kicad_pcb` | small 4-layer board used by the KiCad plugin's capability check | |
| `bwsma_air_launch_coupon.kicad_pcb` | 50 × 12 mm GCPW coupon with the "air" launch for the BAT Wireless BWSMA-KE-P001 edge SMA (LCSC C496550) | `python ../twin/make_panel.py out --sim 1` (connector model included) |
| `silabs-ams915/` | tools and results for the Silicon Labs MSC-AMS915 915 MHz antenna kit (NanoVNA measurement, PADS → KiCad import, Z2 EMS results) | see its README |

All boards open in KiCad 10. With the KiCad plugin installed, open a board and press the Z2Labs EMS buttons:
stackup preset, mesh preview, simulate, field overlay.
