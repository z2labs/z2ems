"""WES0115 IFA: the PADS GND shorting stub (3.0 x 0.46 mm at the arc's radial leg) got the antenna net when the
round-board copper was restored (GND net lookup returned 0, then the no-net->antenna-net step renamed it).
Set it back to GND and refill, so the GND pour connects to it."""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import pcbnew
f=os.path.join(ROOT, "kicad/boards/WES0115-01-APF915P_panic-IFA.kicad_pcb")
b=pcbnew.LoadBoard(f); gnd=b.GetNetInfo().GetNetItem('GND'); print('GND code', gnd.GetNetCode())
for d in b.GetDrawings():
    if d.GetClass()=='PCB_SHAPE' and d.GetLayer()==pcbnew.F_Cu:
        bb=d.GetBoundingBox()
        if abs(bb.GetLeft()/1e6-34.35)<0.1 and abs(bb.GetTop()/1e6-20.31)<0.1:
            d.SetNetCode(gnd.GetNetCode()); print('stub re-netted to', d.GetNetname(), 'bbox %.2f..%.2f x %.2f..%.2f'%(bb.GetLeft()/1e6,bb.GetRight()/1e6,bb.GetTop()/1e6,bb.GetBottom()/1e6))
zs=pcbnew.ZONES(); [zs.append(z) for z in b.Zones() if not z.GetIsRuleArea()]
pcbnew.ZONE_FILLER(b).Fill(zs)
import sys  # needs shapely in the KiCad python (pip install --target <dir> shapely; PYTHONPATH=<dir>)
from shapely.geometry import Polygon; from shapely.ops import unary_union
def polys(sps):
    out=[]
    for i in range(sps.OutlineCount()):
        o=sps.Outline(i); pts=[(o.CPoint(k).x/1e6,o.CPoint(k).y/1e6) for k in range(o.PointCount())]
        if len(pts)>=3: out.append(Polygon(pts).buffer(0))
    return out
G=unary_union([p for z in b.Zones() if z.GetNetname()=='GND' and z.IsOnLayer(pcbnew.F_Cu) for p in polys(z.GetFilledPolysList(pcbnew.F_Cu))])
for d in b.GetDrawings():
    if d.GetClass()=='PCB_SHAPE' and d.GetLayer()==pcbnew.F_Cu:
        sps=pcbnew.SHAPE_POLY_SET(); d.TransformShapeToPolygon(sps,pcbnew.F_Cu,0,int(0.005e6),pcbnew.ERROR_INSIDE)
        P=unary_union(polys(sps)); print(' poly net %-7s dist to GND fill %.3f mm'%(d.GetNetname(), P.distance(G)))
pcbnew.SaveBoard(f+'.tmp',b); import os; os.replace(f+'.tmp',f); print('saved')
