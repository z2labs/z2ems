#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2026 Zoltan Doczi -- GPL-3.0-or-later
"""
bwsma_model: 3D model (STEP + VRML) of the BAT Wireless BWSMA-KE-P001 (LCSC C496550), a cheap
straight edge-mount SMA jack with an offset pin, from its datasheet drawing (July-2-2025-A01):

    flange 6.5 x 6.5 x 1.5, barrel 8 mm (1/4-36 UNEF, 6.25), total length 13.5
    centre pin d 0.95, four 0.9 x 0.9 legs 4 mm long; two legs in the plane of the pin
    (y = +-2.55, they sit on the top copper), two on the bottom side, 1.7 mm between the
    two leg rows (1.6 mm board)

Model frame (= footprint frame of make_panel.Board._bwsma): origin on the board top surface at the
board edge, on the pin axis line; +x into the board, the barrel points to -x, z up (board bottom
at z = -1.6). Needs cadquery.

    python bwsma_model.py out_dir      # writes BWSMA-KE-P001.step and .wrl
"""
import math
import os
import sys

import cadquery as cq

PIN_D = 0.95
LEG = 0.9               # leg cross-section (square)
LEG_Y = 2.55            # leg centres, +-y
LEG_L = 4.0             # pin and leg length beyond the flange face
GAP = 1.7               # between the top and the bottom leg row (board slot)
FLANGE = 6.5
FLANGE_T = 1.5
BARREL_L = 8.0
THREAD_D = 6.25         # 1/4-36 UNEF
THREAD_P = 25.4 / 36
THREAD_H = 0.38          # thread depth (0.54 P, UN external)
NECK_D, NECK_L = 5.0, 1.0
BORE_D = 4.6            # SMA jack outer-conductor bore at the mating face
PTFE_D = 4.1            # PTFE dielectric through the flange (back face: PTFE ring around the pin)
SOCKET_D = 1.27         # socket contact
NAME = 'BWSMA-KE-P001'

ZC = PIN_D / 2          # pin axis above the top copper


def _cyl_x(d, x0, x1, y=0.0, z=ZC):
    """ Cylinder of diameter d along x from x0 to x1 (x1 > x0). """
    return cq.Workplane('YZ', origin=(x0, y, z)).circle(d / 2).extrude(x1 - x0)


def _thread(x0, x1):
    """ 1/4-36 UNEF external thread from x0 (mating face) to x1: a minor-diameter core plus a helical 60 deg
        ridge swept along the barrel, trimmed by a chamfered major-diameter envelope (0.3 mm lead-in at the face,
        thread runs out into the neck). """
    L = x1 - x0
    h = THREAD_H                                      # radial depth of the ridge
    r_min = THREAD_D / 2 - h
    half = math.tan(math.radians(30)) * h             # half width of the ridge at the root grows with depth
    flat = THREAD_P / 8                               # crest flat
    root = flat / 2 + half
    # profile in the XZ plane at the root radius (helix along +Z), slightly sunk into the core for a clean union
    prof = (cq.Workplane('XZ').polyline([(r_min - 0.05, -root - 0.02), (r_min, -root), (r_min + h - 0.01, -flat / 2),
                                          (r_min + h - 0.01, flat / 2), (r_min, root), (r_min - 0.05, root + 0.02)]).close())
    helix = cq.Wire.makeHelix(THREAD_P, L + 2 * THREAD_P, r_min - 0.05)
    ridge = prof.sweep(cq.Workplane().add(helix), isFrenet=True).translate((0, 0, -THREAD_P))
    core = cq.Workplane('XY').circle(r_min + 0.01).extrude(L)
    env = cq.Workplane('XY').circle(THREAD_D / 2).extrude(L).faces('<Z').chamfer(0.35)
    t = cq.Workplane('XY').newObject([core.val().fuse(ridge.val().intersect(env.val())).clean()])   # crest 0.01 inside the envelope: no coincident faces
    # helix axis Z -> barrel axis X, start at the mating face
    return t.rotate((0, 0, 0), (0, 1, 0), 90).translate((x0, 0, ZC))


def build():
    x_f = -FLANGE_T                                   # flange back face (barrel side)
    x_tip = x_f - BARREL_L                            # front face of the jack
    # body: flange + neck + threaded barrel with 36 TPI grooves, front chamfer, bore
    flange = cq.Workplane('XY').box(FLANGE_T, FLANGE, FLANGE).translate((-FLANGE_T / 2, 0, ZC)).edges('|X').fillet(0.3)
    neck = _cyl_x(NECK_D, x_f - NECK_L, x_f)
    thread = _thread(x_tip, x_f - NECK_L)
    body = flange.union(neck).union(thread)
    # coax through the whole body: outer-conductor bore at the mating face, PTFE bore through the flange
    body = body.cut(_cyl_x(BORE_D, x_tip - 0.1, x_tip + 3.0)).cut(_cyl_x(PTFE_D, x_tip, 0.1))
    # back face: shallow square pocket leaving a 0.6 mm rim (as on the part), the PTFE disc sits in it
    body = body.cut(cq.Workplane('XY').box(0.2, FLANGE - 1.2, FLANGE - 1.2).translate((-0.1 + 0.001, 0, ZC)))
    # PTFE: fills the bore up to 0.3 mm behind the mating face and runs through the flange, flush with the pocket floor
    ptfe = _cyl_x(BORE_D, x_tip + 0.3, x_tip + 3.0).union(_cyl_x(PTFE_D, x_tip + 3.0, -0.2))
    ptfe = ptfe.cut(_cyl_x(SOCKET_D, x_tip, x_tip + 4.0)).cut(_cyl_x(PIN_D, x_tip + 3.9, 0.1))
    # socket contact: slotted tube at the face, its hole visible in the PTFE; it carries the pin
    socket = _cyl_x(SOCKET_D, x_tip + 0.5, x_tip + 4.0).cut(_cyl_x(0.92, x_tip + 0.4, x_tip + 3.0))
    for a in (0, 90):                                 # four slots
        slot = cq.Workplane('XY').box(2.0, 0.18, 1.6).translate((x_tip + 1.4, 0, ZC)).rotate((0, 0, ZC), (1, 0, ZC), a)
        socket = socket.cut(slot)
    # centre pin on the top copper, rounded tip
    pin = _cyl_x(PIN_D, x_tip + 3.9, LEG_L - PIN_D / 2).union(cq.Workplane('XY').sphere(PIN_D / 2).translate((LEG_L - PIN_D / 2, 0, ZC)))
    # legs: top row in the pin plane (bottom face on the top copper), bottom row below the board slot
    legs = None
    for zb in (0.0, -GAP - LEG):
        for sy in (-1, 1):
            leg = cq.Workplane('XY').box(LEG_L + 0.2, LEG, LEG).translate((LEG_L / 2 - 0.1, sy * LEG_Y, zb + LEG / 2))   # rooted in the back pocket.faces('>X').edges().chamfer(0.15)
            legs = leg if legs is None else legs.union(leg)
    body = body.union(legs)
    gold = cq.Color(0.86, 0.68, 0.24)
    assy = cq.Assembly(name=NAME)
    assy.add(body, name='body', color=gold)
    assy.add(ptfe, name='ptfe', color=cq.Color(0.95, 0.95, 0.93))
    assy.add(socket, name='socket', color=cq.Color(0.93, 0.75, 0.30))
    assy.add(pin, name='pin', color=cq.Color(0.93, 0.75, 0.30))
    return assy


def write(out):
    os.makedirs(out, exist_ok=True)
    assy = build()
    step = os.path.join(out, NAME + '.step')
    assy.save(step, 'STEP')
    wrl = os.path.join(out, NAME + '.wrl')
    # KiCad VRML units are 0.1 inch: scale the tessellated assembly by 1/2.54
    from cadquery.occ_impl.exporters.assembly import exportVRML
    scaled = cq.Assembly(name=NAME)
    for ch in assy.children:
        scaled.add(ch.obj.val().scale(1 / 2.54) if hasattr(ch.obj, 'val') else ch.obj, name=ch.name, color=ch.color)
    exportVRML(scaled, wrl, tolerance=0.005, angularTolerance=0.2)
    return step, wrl


if __name__ == '__main__':
    print('\n'.join(write(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models'))))
