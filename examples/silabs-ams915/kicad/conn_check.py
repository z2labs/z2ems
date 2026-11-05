"""Connectivity check: for each board, build a union of copper per layer (pads, tracks, filled polygons, zones)
and report which copper islands the SMA signal pad reaches, and where the antenna polygons sit relative to it."""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import sys  # needs shapely in the KiCad python (pip install --target <dir> shapely; PYTHONPATH=<dir>)
import pcbnew, glob, os, sys
from shapely.geometry import Polygon, MultiPolygon, Point
from shapely.ops import unary_union
MM=1e6
def poly_of(sps):
    out=[]
    for i in range(sps.OutlineCount()):
        o=sps.Outline(i); pts=[(o.CPoint(k).x/MM, o.CPoint(k).y/MM) for k in range(o.PointCount())]
        if len(pts)>=3: out.append(Polygon(pts).buffer(0))
    return out
for f in sorted(glob.glob(os.path.join(ROOT, "kicad/boards/*.kicad_pcb"))):
    b=pcbnew.LoadBoard(f); name=os.path.basename(f)[:7]
    if name not in ('WES0111','WES0115','WES0116','WES0117','WES0118'): continue
    sma=[fp for fp in b.GetFootprints() if fp.GetReference().startswith('SMA')][0]
    sig=[p for p in sma.Pads() if p.GetNetname()!='GND'][0]
    L=pcbnew.F_Cu
    items=[]   # (label, polygon)
    for fp in b.GetFootprints():
        for p in fp.Pads():
            if p.IsOnLayer(L):
                items.append(('pad %s.%s[%s]'%(fp.GetReference(),p.GetNumber(),p.GetNetname()), unary_union(poly_of(p.GetEffectivePolygon(L)))))
    for t in b.GetTracks():
        if t.GetClass()=='PCB_TRACK' and t.GetLayer()==L:
            sps=pcbnew.SHAPE_POLY_SET(); t.TransformShapeToPolygon(sps,L,0,int(0.005*MM),pcbnew.ERROR_INSIDE)
            items.append(('track[%s]'%t.GetNetname(), unary_union(poly_of(sps))))
    for d in b.GetDrawings():
        if d.GetClass()=='PCB_SHAPE' and d.GetLayer()==L:
            sps=pcbnew.SHAPE_POLY_SET(); d.TransformShapeToPolygon(sps,L,0,int(0.005*MM),pcbnew.ERROR_INSIDE)
            bb=d.GetBoundingBox(); items.append(('poly[%s] %.0fx%.0f'%(d.GetNetname(),bb.GetWidth()/MM,bb.GetHeight()/MM), unary_union(poly_of(sps))))
    for z in b.Zones():
        if not z.GetIsRuleArea() and z.IsOnLayer(L):
            items.append(('zone[%s]'%z.GetNetname(), unary_union(poly_of(z.GetFilledPolysList(L)))))
    items=[(l,g) for l,g in items if g and not g.is_empty]
    # union-find by touching (tolerance 1 um)
    n=len(items); par=list(range(n))
    def fnd(i):
        while par[i]!=i: par[i]=par[par[i]]; i=par[i]
        return i
    for i in range(n):
        for j in range(i+1,n):
            if items[i][1].distance(items[j][1])<0.001: par[fnd(i)]=fnd(j)
    ser={'WES0111':'RS1','WES0115':'RS5','WES0116':'RS6','WES0117':'LS7','WES0118':'LS8'}[name]
    sfp=[fp for fp in b.GetFootprints() if fp.GetReference()==ser][0]
    apad=[p for p in sfp.Pads() if p.GetNetname()!=sig.GetNetname()][0]
    print('== antenna side of',ser,'pad',apad.GetNumber(),'net',apad.GetNetname())
    start=next(i for i,(l,g) in enumerate(items) if l.startswith('pad %s.%s'%(ser,apad.GetNumber())))
    root=fnd(start); grp=[items[i] for i in range(n) if fnd(i)==root]
    print('==',name,'signal pad',sig.GetNetname(),'| island members:',len(grp))
    kinds={}
    for l,g in grp: kinds[l.split('[')[0].split(' ')[0]]=kinds.get(l.split('[')[0].split(' ')[0],0)+1
    print('   island kinds:',kinds, '| pads in island:', sorted(set(l.split()[1] for l,g in grp if l.startswith('pad')))[:12])
    island=unary_union([g for l,g in grp])
    for l,g in items:
        if l.startswith('poly') and fnd(items.index((l,g)))!=root:
            print('   NOT connected: %s  gap to island %.3f mm'%(l, g.distance(island)))
    big=[(l,g) for l,g in grp if l.startswith('poly')]
    print('   antenna polys in island:', [l for l,g in big])
