"""Write a 50 Ohm microstrip test coupon as a KiCad board (JLC04161H-7628 stackup).

    python make_coupon.py out.kicad_pcb [--width 0.3726] [--length 40]

F.Cu: one straight track between two SMD launch pads at the board edges (refs J1, J2).
In1.Cu and B.Cu: filled GND planes, tied by two rows of stitching vias 2.5 mm off the track. The stackup section is the JLCPCB 4-layer
JLC04161H-7628 stackup (1 oz outer, 0.5 oz inner).
"""
import argparse
import uuid

STACKUP = '''    (stackup
      (layer "F.SilkS" (type "Top Silk Screen"))
      (layer "F.Mask" (type "Top Solder Mask") (thickness 0.01524) (epsilon_r 3.8) (loss_tangent 0))
      (layer "F.Cu" (type "copper") (thickness 0.035))
      (layer "dielectric 1" (type "prepreg") (thickness 0.2104) (material "7628 NP-155F") (epsilon_r 4.4) (loss_tangent 0.02))
      (layer "In1.Cu" (type "copper") (thickness 0.0152))
      (layer "dielectric 2" (type "core") (thickness 1.065) (material "Core NP-155F") (epsilon_r 4.43) (loss_tangent 0.02))
      (layer "In2.Cu" (type "copper") (thickness 0.0152))
      (layer "dielectric 3" (type "prepreg") (thickness 0.2104) (material "7628 NP-155F") (epsilon_r 4.4) (loss_tangent 0.02))
      (layer "B.Cu" (type "copper") (thickness 0.035))
      (layer "B.Mask" (type "Bottom Solder Mask") (thickness 0.01524) (epsilon_r 3.8) (loss_tangent 0))
      (layer "B.SilkS" (type "Bottom Silk Screen"))
      (copper_finish "ENIG")
      (dielectric_constraints no))'''


def U():
    return str(uuid.uuid4())


def board(w, L, H=14.0, pad_len=1.0, via_pitch=2.5, via_offset=2.5):
    y = H / 2
    vias = ''
    x = via_pitch
    while x < L - via_pitch / 2:      # GND stitching vias tie In1 and B.Cu together (no parallel-plate cavity)
        for yy in (y - via_offset, y + via_offset):
            vias += '  (via (at %g %g) (size 0.6) (drill 0.3) (layers "F.Cu" "B.Cu") (net 1) (uuid "%s"))\n' % (x, yy, U())
        x += via_pitch
    rect = lambda layer: '(xy 0 0) (xy %g 0) (xy %g %g) (xy 0 %g)' % (L, L, H, H)
    zones = ''
    for layer in ('In1.Cu', 'B.Cu'):
        zones += '''  (zone (net 1) (net_name "GND") (layer "%s") (uuid "%s") (hatch edge 0.5)
    (connect_pads (clearance 0.2)) (min_thickness 0.2) (filled_areas_thickness no)
    (fill yes (thermal_gap 0.5) (thermal_bridge_width 0.5))
    (polygon (pts %s))
    (filled_polygon (layer "%s") (pts %s)))
''' % (layer, U(), rect(layer), layer, rect(layer))
    fps = ''
    for ref, x in (('J1', pad_len / 2), ('J2', L - pad_len / 2)):
        fps += '''  (footprint "z2labs:RF_Launch_Pad" (layer "F.Cu") (uuid "%s") (at %g %g)
    (property "Reference" "%s" (at 0 -1.5 0) (layer "F.SilkS") (uuid "%s") (effects (font (size 0.8 0.8) (thickness 0.12))))
    (property "Value" "SMA edge launch" (at 0 1.5 0) (layer "F.Fab") (uuid "%s") (effects (font (size 0.8 0.8) (thickness 0.12))))
    (attr smd)
    (pad "1" smd rect (at 0 0) (size %g %g) (layers "F.Cu" "F.Mask") (net 2 "RF") (uuid "%s"))
  )
''' % (U(), x, y, ref, U(), U(), pad_len, w, U())
    return '''(kicad_pcb (version 20241229) (generator "pcbnew") (generator_version "9.0")
  (general (thickness 1.6) (legacy_teardrops no))
  (paper "A4")
  (title_block (title "50 Ohm microstrip coupon, JLC04161H-7628") (company "Z2Labs"))
  (layers
    (0 "F.Cu" signal) (4 "In1.Cu" signal "In1_GND") (6 "In2.Cu" signal) (2 "B.Cu" signal)
    (9 "F.Adhes" user "F.Adhesive") (11 "B.Adhes" user "B.Adhesive") (13 "F.Paste" user) (15 "B.Paste" user)
    (5 "F.SilkS" user "F.Silkscreen") (7 "B.SilkS" user "B.Silkscreen") (1 "F.Mask" user) (3 "B.Mask" user)
    (25 "Edge.Cuts" user) (27 "Margin" user) (31 "F.CrtYd" user "F.Courtyard") (29 "B.CrtYd" user "B.Courtyard")
    (35 "F.Fab" user) (33 "B.Fab" user))
  (setup
%s
    (pad_to_mask_clearance 0)
    (pcbplotparams (layerselection 0x00010fc_ffffffff) (outputdirectory "")))
  (net 0 "") (net 1 "GND") (net 2 "RF")
%s  (segment (start %g %g) (end %g %g) (width %g) (layer "F.Cu") (net 2) (uuid "%s"))
  (gr_rect (start 0 0) (end %g %g) (stroke (width 0.05) (type default)) (fill no) (layer "Edge.Cuts") (uuid "%s"))
%s%s)
''' % (STACKUP, fps, pad_len / 2, y, L - pad_len / 2, y, w, U(), L, H, U(), vias, zones)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--width', type=float, default=0.3726)
    ap.add_argument('--length', type=float, default=40.0)
    a = ap.parse_args()
    open(a.out, 'w').write(board(a.width, a.length))
    print('wrote', a.out)
