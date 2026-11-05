import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import pcbnew
SRC=os.path.join(ROOT, "kicad/WES0110-01-panel-fixed2.kicad_pcb")
DST=os.path.join(ROOT, "_extract/panel_noframe.kicad_pcb")
b=pcbnew.LoadBoard(SRC)
for z in list(b.Zones()):
    if z.GetNetname()=="$$$10": b.Delete(z)
pcbnew.SaveBoard(DST,b); print("ok")
