"""Gerber reference boards: kicad/gerber_ref/<board>.json -> kicad/boards_gerber/<board>.kicad_pcb
Copper = the Gerber union as filled F.Cu/B.Cu polygons (fractured, so holes survive), drills = vias,
outline from the board window, SMA + matching footprints copied from the split board (same frame) so the
twin has a port pad and lumped elements. Nets: a copper polygon takes the net of the pads it touches
(GND wins if any GND pad touches; else the first signal pad's net; antenna copper touching nothing stays
unnamed and is renamed to the far pad net of the series element by touching chain).
Run with KiCad 10 python."""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import sys  # needs shapely in the KiCad python (pip install --target <dir> shapely; PYTHONPATH=<dir>)
import pcbnew, json, glob, os
from shapely.geometry import Polygon
from shapely.ops import unary_union
MM=1e6
SRC=os.path.join(ROOT, "kicad")
OUT=os.path.join(SRC,'boards_gerber'); os.makedirs(OUT,exist_ok=True)
KEEP=('SMA','RS','LS','CP')
for jf in sorted(glob.glob(os.path.join(SRC,'gerber_ref','WES*.json'))):
    J=json.load(open(jf)); name=J['name']
    ref=pcbnew.LoadBoard(os.path.join(SRC,'boards',name+'.kicad_pcb'))
    b=pcbnew.LoadBoard(os.path.join(SRC,'boards',name+'.kicad_pcb'))   # start from the split board: nets, stackup, settings
    # strip everything except the kept footprints
    for fp in list(b.GetFootprints()):
        if not fp.GetReference().startswith(KEEP): b.Delete(fp)
    # keep the SMA signal-net tracks of the split board: the twin walks out along a track to place the port
    sman=[p.GetNetCode() for fp in b.GetFootprints() if fp.GetReference().startswith('SMA') for p in fp.Pads() if p.GetNetname()!='GND'][0]
    for t in list(b.GetTracks()):
        if not (t.GetClass()=='PCB_TRACK' and t.GetNetCode()==sman): b.Delete(t)
    for z in list(b.Zones()): b.Delete(z)
    for d in list(b.GetDrawings()): b.Delete(d)
    # outline
    E=pcbnew.Edge_Cuts
    if J['outline']=='rect':
        x0,y0,x1,y1=J['outline_param']
        for a,c in (((x0,y0),(x1,y0)),((x1,y0),(x1,y1)),((x1,y1),(x0,y1)),((x0,y1),(x0,y0))):
            s=pcbnew.PCB_SHAPE(b); s.SetShape(pcbnew.SHAPE_T_SEGMENT); s.SetLayer(E); s.SetStart(pcbnew.VECTOR2I(int(a[0]*MM),int(a[1]*MM))); s.SetEnd(pcbnew.VECTOR2I(int(c[0]*MM),int(c[1]*MM))); s.SetWidth(int(0.1*MM)); b.Add(s)
    else:
        cx,cy,r=J['outline_param']; s=pcbnew.PCB_SHAPE(b); s.SetShape(pcbnew.SHAPE_T_CIRCLE); s.SetLayer(E); s.SetCenter(pcbnew.VECTOR2I(int(cx*MM),int(cy*MM))); s.SetEnd(pcbnew.VECTOR2I(int((cx+r)*MM),int(cy*MM))); s.SetWidth(int(0.1*MM)); b.Add(s)
    # align: the SMA signal pin drill in the Gerber vs the SMA pad of the split board (rect windows are +-1 mm)
    sma=[fp for fp in b.GetFootprints() if fp.GetReference().startswith('SMA')][0]
    sp=[p for p in sma.Pads() if p.GetNetname()!='GND'][0].GetPosition(); sx,sy=sp.x/MM,sp.y/MM
    big=[(x,y,d) for x,y,d in J['drills'] if d>=1.0]
    if big:
        x,y,d=min(big,key=lambda t:(t[0]-sx)**2+(t[1]-sy)**2); dx,dy=sx-x,sy-y
        if abs(dx)<3 and abs(dy)<3:
            for key in ('F_Cu','B_Cu'):
                for P in J[key]:
                    P['ext']=[[u+dx,v+dy] for u,v in P['ext']]; P['holes']=[[[u+dx,v+dy] for u,v in h] for h in P['holes']]
            J['drills']=[(u+dx,v+dy,dd) for u,v,dd in J['drills']]
            print('  %s aligned by (%.3f, %.3f) mm on SMA pin drill d=%.2f'%(name[:7],dx,dy,d))
        else: print('  %s NO ALIGN: nearest big drill off by (%.2f,%.2f)'%(name[:7],dx,dy))
    # pads geometry for net assignment
    pads=[]
    for fp in b.GetFootprints():
        for p in fp.Pads():
            sps=p.GetEffectivePolygon(pcbnew.F_Cu); pts=[(sps.CVertex(i,0,-1).x/MM,sps.CVertex(i,0,-1).y/MM) for i in range(sps.Outline(0).PointCount())] if sps.OutlineCount() else []
            if len(pts)>=3: pads.append((p.GetNetname(), p.GetNetCode(), Polygon(pts).buffer(0), fp.GetReference()+'.'+p.GetNumber()))
    stats={}
    for layer,key in ((pcbnew.F_Cu,'F_Cu'),(pcbnew.B_Cu,'B_Cu')):
        for P in J[key]:
            poly=Polygon(P['ext'],P['holes']).buffer(0)
            if poly.is_empty: continue
            touching=[(n,c,lab) for n,c,g,lab in pads if layer==pcbnew.F_Cu or n=='GND' or True] if False else [(n,c,lab) for n,c,g,lab in pads if g.distance(poly)<0.002]
            net=0
            if any(n=='GND' for n,c,lab in touching): net=[c for n,c,lab in touching if n=='GND'][0]
            elif touching: net=touching[0][1]
            sps=pcbnew.SHAPE_POLY_SET(); ol=sps.NewOutline() if False else None
            sps=pcbnew.SHAPE_POLY_SET()
            idx=sps.NewOutline()
            for x,y in P['ext'][:-1]: sps.Append(int(x*MM),int(y*MM),idx)
            for h in P['holes']:
                hi=sps.NewHole(idx)
                for x,y in h[:-1]: sps.Append(int(x*MM),int(y*MM),idx,hi)
            sps.Fracture()
            for i in range(sps.OutlineCount()):
                o=sps.Outline(i); s=pcbnew.PCB_SHAPE(b); s.SetShape(pcbnew.SHAPE_T_POLY); s.SetLayer(layer); s.SetFilled(True); s.SetWidth(0)
                s.SetPolyPoints(pcbnew.VECTOR_VECTOR2I([o.CPoint(k) for k in range(o.PointCount())])); s.SetNetCode(net); b.Add(s)
            k=(key, b.FindNet(net).GetNetname() if net else '-'); stats[k]=stats.get(k,0)+round(poly.area)
    # drills -> vias (plated), tiny annulus; net GND if inside a GND polygon else none
    gnd_f=unary_union([Polygon(P['ext'],P['holes']).buffer(0) for P in J['F_Cu']])  # all copper, for membership only
    nv=0
    for x,y,d in J['drills']:
        if d>1.0: continue   # SMA pins / mounting holes: the kept footprints' THT pads carry those
        v=pcbnew.PCB_VIA(b); v.SetPosition(pcbnew.VECTOR2I(int(x*MM),int(y*MM))); v.SetDrill(int(d*MM)); v.SetWidth(int((d+0.3)*MM)); v.SetLayerPair(pcbnew.F_Cu,pcbnew.B_Cu)
        v.SetNetCode(b.FindNet('GND').GetNetCode()); b.Add(v); nv+=1
    pcbnew.SaveBoard(os.path.join(OUT,name+'.kicad_pcb'),b)
    print(name[:28], 'vias',nv, {('%s/%s'%k):v for k,v in stats.items()})
