"""Split the imported WES0110-01 PADS panel into 9 single-antenna KiCad boards.
Run with KiCad 10 python:  python split_panel.py
- PADS layer 24 (documentation) lines landed on In1.Cu: they are the routed board
  outlines with mouse-bite tab gaps. We cluster them, take the convex hull of the
  cluster vertices as the real Edge.Cuts outline, and keep every item whose bbox
  lies inside that hull bbox.
"""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import pcbnew, os, math
from collections import defaultdict

SRC = os.path.join(ROOT, "kicad/WES0110-01-panel-fixed2.kicad_pcb")
OUT = os.path.join(ROOT, "kicad/boards")
NAMES = {1:"WES0111-01-APL915M_medium-ILA", 2:"WES0112-01-ACM915D_ceramic",
         3:"WES0113-01-AWH915S_small-helical", 4:"WES0114-01-AWH915M_medium-helical",
         5:"WES0115-01-APF915P_panic-IFA", 6:"WES0116-01-APL915P_panic-ILA",
         7:"WES0117-01-APN915D_meander", 8:"WES0118-01-APL915S_small-ILA",
         9:"WES0119-01-APB915D_BIFA"}
MM = 1e6
GAP = 4.0*MM     # max tab gap to bridge when clustering outline segments

def hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3: return pts
    def cross(o,a,b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lo=[]
    for p in pts:
        while len(lo)>=2 and cross(lo[-2],lo[-1],p)<=0: lo.pop()
        lo.append(p)
    up=[]
    for p in reversed(pts):
        while len(up)>=2 and cross(up[-2],up[-1],p)<=0: up.pop()
        up.append(p)
    return lo[:-1]+up[:-1]

def main():
    b = pcbnew.LoadBoard(SRC)
    in1 = b.GetLayerID("In1.Cu")
    lines = [d for d in b.GetDrawings() if d.GetLayer()==in1 and d.GetShapeStr()=="Line"]
    # cluster by endpoint proximity
    segs = [((d.GetStart().x,d.GetStart().y),(d.GetEnd().x,d.GetEnd().y)) for d in lines]
    parent = list(range(len(segs)))
    def find(i):
        while parent[i]!=i: parent[i]=parent[parent[i]]; i=parent[i]
        return i
    def near(p,q): return math.hypot(p[0]-q[0],p[1]-q[1]) <= GAP
    for i in range(len(segs)):
        for j in range(i+1,len(segs)):
            if any(near(p,q) for p in segs[i] for q in segs[j]):
                parent[find(i)] = find(j)
    clusters = defaultdict(list)
    for i in range(len(segs)): clusters[find(i)].append(i)
    clusters = [c for c in clusters.values() if len(c) >= 4]
    # (kind, params in mm, KiCad panel coords): 0113 = 40x40 square placed like 0114 (SMA 20 mm from left, 30 mm from top);
    # 0115/0116 = Ø30 round boards centred on their SMA (PADS COPCUT CIRCLE r=15 mm, origin == SMA)
    EXTRA = {3:("rect",(-151.9,133.51,-111.9,173.51)), 5:("circle",(-86.90,168.51,15.0)), 6:("circle",(-86.95,208.52,15.0))}
    print("outline clusters:", len(clusters))
    smas = {f.GetReference(): f.GetPosition() for f in b.GetFootprints() if f.GetReference().startswith("SMA")}

    report = []
    jobs = [("cluster",idx) for idx in clusters] + [("extra",(n,k,v)) for n,(k,v) in EXTRA.items()]
    for ci, (jt, idx) in enumerate(jobs):
      if jt=="extra":
        n,kind,v = idx
        if kind=="rect": x0,y0,x1,y1=[int(a*MM) for a in v]; hp=[(x0,y0),(x1,y0),(x1,y1),(x0,y1)]; circ=None
        else: cx,cy,r=v; circ=(int(cx*MM),int(cy*MM),int(r*MM)); x0,y0,x1,y1=circ[0]-circ[2],circ[1]-circ[2],circ[0]+circ[2],circ[1]+circ[2]; hp=None
        xs=[x0,x1]; ys=[y0,y1]
      else:
        circ=None
        pts = [p for i in idx for p in segs[i]]
        hp = hull(pts)
        xs=[p[0] for p in hp]; ys=[p[1] for p in hp]
        I=int(1.0*MM)  # routing path is 1 mm outside the real edge
        x0,x1,y0,y1=min(xs)+I,max(xs)-I,min(ys)+I,max(ys)-I
        hp=[(x0,y0),(x1,y0),(x1,y1),(x0,y1)]; xs=[x0,x1]; ys=[y0,y1]
      box = pcbnew.BOX2I(pcbnew.VECTOR2I(min(xs),min(ys)), pcbnew.VECTOR2I(max(xs)-min(xs), max(ys)-min(ys)))
      box.Inflate(int(0.3*MM))
      sma = [r for r,p in smas.items() if box.Contains(p)]
      if len(sma)!=1:
          print("cluster",ci,"SMA match",sma,"-> skipped"); continue
      n = int(sma[0][3:]); name = NAMES[n]
      nb = pcbnew.LoadBoard(SRC)
      edge = pcbnew.Edge_Cuts
      # delete everything outside box
      for f in list(nb.GetFootprints()):
          if not box.Contains(f.GetPosition()): nb.Delete(f)
      for t in list(nb.GetTracks()):
          if not (box.Contains(t.GetStart()) and box.Contains(t.GetEnd())): nb.Delete(t)
      for z in list(nb.Zones()):
          zb = z.GetBoundingBox()
          if not (box.Contains(zb.GetOrigin()) and box.Contains(zb.GetEnd())): nb.Delete(z)
      for d in list(nb.GetDrawings()):
          db = d.GetBoundingBox()
          keep = box.Contains(db.GetOrigin()) and box.Contains(db.GetEnd())
          if not keep or d.GetLayer()==in1: nb.Delete(d)
      # closed outline
      if circ:
          c = pcbnew.PCB_SHAPE(nb); c.SetShape(pcbnew.SHAPE_T_CIRCLE); c.SetLayer(edge)
          c.SetCenter(pcbnew.VECTOR2I(circ[0],circ[1])); c.SetEnd(pcbnew.VECTOR2I(circ[0]+circ[2],circ[1])); c.SetWidth(int(0.1*MM)); nb.Add(c)
      for k in range(len(hp) if hp else 0):
          s = pcbnew.PCB_SHAPE(nb); s.SetShape(pcbnew.SHAPE_T_SEGMENT); s.SetLayer(edge)
          s.SetStart(pcbnew.VECTOR2I(*hp[k])); s.SetEnd(pcbnew.VECTOR2I(*hp[(k+1)%len(hp)]))
          s.SetWidth(int(0.1*MM)); nb.Add(s)
      # move so board origin = bbox top-left
      off = pcbnew.VECTOR2I(-min(xs)+int(10*MM), -min(ys)+int(10*MM))
      for it in list(nb.GetFootprints())+list(nb.GetTracks())+list(nb.Zones())+list(nb.GetDrawings()):
          it.Move(off)
      # stackup: 2 layer, 1.55 mm (AN848)
      ds = nb.GetDesignSettings(); ds.SetCopperLayerCount(2); ds.SetBoardThickness(int(1.55*MM))

      zs=[z for z in nb.Zones() if z.GetZoneName().startswith("GND_round")]
      if zs:
          zv=pcbnew.ZONES(); [zv.append(z) for z in zs]
          pcbnew.ZONE_FILLER(nb).Fill(zv)
          print("  round GND fill:", [round(z.GetFilledArea()/1e12,1) for z in zs], "mm2")
      fn = os.path.join(OUT, name + ".kicad_pcb")
      pcbnew.SaveBoard(fn, nb)
      w=(max(xs)-min(xs))/MM; h=(max(ys)-min(ys))/MM
      report.append((n,name,round(w,2),round(h,2),len(nb.GetFootprints()),len(nb.GetTracks()),len(nb.Zones()),len(hp) if hp else 1))
      print("saved", name, round(w,2),"x",round(h,2),"mm  fp",len(nb.GetFootprints()),"tracks",len(nb.GetTracks()),"zones",len(nb.Zones()), "circle" if circ else "rect")
    with open(os.path.join(OUT,"split_report.txt"),"w") as f:
        for r in report: f.write("\t".join(map(str,r))+"\n")

if __name__=="__main__": main()
