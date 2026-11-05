"""Fix KiCad PADS-importer artefacts in the WES0110-01 panel.
PADS level-1 COPPER (COPCLS) polygons - the printed antenna shapes, feed patches and
small GND patches - land on In1.Cu as zones. They are solid copper shapes, so convert
them to filled PCB_SHAPE polygons on F.Cu with their net kept. Output: *-panel-fixed.kicad_pcb
"""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import pcbnew
SRC=os.path.join(ROOT, "kicad/WES0110-01-panel.kicad_pcb")
DST=os.path.join(ROOT, "kicad/WES0110-01-panel-fixed.kicad_pcb")
b=pcbnew.LoadBoard(SRC); in1=b.GetLayerID("In1.Cu"); n=0
for z in list(b.Zones()):
    if z.GetIsRuleArea() or z.GetFirstLayer()!=in1: continue
    o=z.Outline()
    for oi in range(o.OutlineCount()):
        s=pcbnew.PCB_SHAPE(b); s.SetShape(pcbnew.SHAPE_T_POLY); s.SetLayer(pcbnew.F_Cu)
        s.SetFilled(True); s.SetWidth(0)
        pts=pcbnew.VECTOR_VECTOR2I([o.CVertex(i,oi,-1) for i in range(o.Outline(oi).PointCount())])
        s.SetPolyPoints(pts); s.SetNetCode(z.GetNetCode()); b.Add(s); n+=1
    b.Delete(z)
print("converted In1.Cu copper polygons -> F.Cu filled shapes:", n)
pcbnew.SaveBoard(DST,b)
