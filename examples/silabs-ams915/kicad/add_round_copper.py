"""Add the copper the KiCad PADS importer dropped on the two round (panic-button) boards:
GND pour circles (PADS POUROUT CIRCLE) and the level-1 COPCLS polygons (arc antenna, feed stub)
read straight from the PADS ASCII. Applied to the fixed panel before splitting."""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import re, pcbnew
ASC=os.path.join(ROOT, "MSC-AMS915/ASCII/WES0110-01_ASCII.asc")
PCB=os.path.join(ROOT, "kicad/WES0110-01-panel-fixed.kicad_pcb")
OUT=os.path.join(ROOT, "kicad/WES0110-01-panel-fixed2.kicad_pcb")
U=1.5e6; OX=201.903822; OY=223.514406   # PADS mm -> KiCad: x-OX, OY-y
CENTERS=[(115.0,55.0),(115.0,15.0)]
def near(x,y): return any(abs(x-cx)<16 and abs(y-cy)<16 for cx,cy in CENTERS)
L=open(ASC,encoding='cp1250',errors='replace').read().split('\n')
b=pcbnew.LoadBoard(PCB); MM=1e6
def P(x,y): return pcbnew.VECTOR2I(int((x-OX)*MM), int((OY-y)*MM))
def netcode(name):
    n=b.GetNetInfo().GetNetItem(name) if name else None
    return n.GetNetCode() if n else 0
added=0; newzones=[]
# 1) COPCLS copper polygons in *LINES*
i=L.index(next(l for l in L if l.startswith('*LINES*'))); j=L.index(next(l for l in L if l.startswith('*VIA*')))
hdr=re.compile(r'^(\S+)\s+(COPPER)\s+(-?\d+)\s+(-?\d+)\s+(\d+)\s*(\S*)')
pc=re.compile(r'^(COPCLS)\s+(\d+)\s+(\d+)\s+(\d+)')
k=i; cur=None
while k<j:
    m=hdr.match(L[k])
    if m: cur=(int(m.group(3))/U,int(m.group(4))/U,m.group(6)); k+=1; continue
    m=pc.match(L[k])
    if m and cur:
        n=int(m.group(2)); lvl=int(m.group(4)); pts=[]
        for q in range(n):
            k+=1; f=L[k].split(); pts.append((cur[0]+int(f[0])/U, cur[1]+int(f[1])/U))
        if near(cur[0],cur[1]) and lvl in (1,2):
            s=pcbnew.PCB_SHAPE(b); s.SetShape(pcbnew.SHAPE_T_POLY); s.SetLayer(pcbnew.F_Cu if lvl==1 else pcbnew.B_Cu)
            s.SetFilled(True); s.SetWidth(0); s.SetPolyPoints(pcbnew.VECTOR_VECTOR2I([P(x,y) for x,y in pts[:-1] if True]))
            s.SetNetCode(netcode(cur[2] if cur[2] and not cur[2].startswith('$$$') else '')); b.Add(s); added+=1
            print('poly',lvl,repr(cur[2]),n,'pts at %.1f,%.1f'%cur[:2])
    k+=1
# 2) circular GND pours in *POUR*
i=L.index(next(l for l in L if l.startswith('*POUR*')))
for k in range(i,len(L)-3):
    m=re.match(r'^(\S+)\s+POUROUT\s+(-?\d+)\s+(-?\d+)\s+\d+\s+\d+\s+\S+\s+(\S+)',L[k])
    if m and L[k+1].startswith('CIRCLE'):
        cx,cy=int(m.group(2))/U,int(m.group(3))/U; lvl=int(L[k+1].split()[4]); r=abs(int(L[k+2].split()[0]))/U
        if not near(cx,cy): continue
        import math
        z=pcbnew.ZONE(b); z.SetLayer(pcbnew.F_Cu if lvl==1 else pcbnew.B_Cu); z.SetNetCode(netcode(m.group(4)))
        z.SetZoneName("GND_round_%d"%lvl)
        ol=z.Outline(); ol.NewOutline()
        for t in range(72):
            a=2*math.pi*t/72; p=P(cx+r*math.cos(a), cy+r*math.sin(a)); ol.Append(p.x,p.y)
        z.SetLocalClearance(int(0.3*MM)); z.SetMinThickness(int(0.25*MM)); z.SetPadConnection(pcbnew.ZONE_CONNECTION_FULL)
        z.SetIsFilled(False); z.SetNeedRefill(True)
        b.Add(z); newzones.append(z); added+=1; print('zone circle',lvl,m.group(4),'c=%.2f,%.2f r=%.3f'%(cx,cy,r))
print('added',added)
zv=pcbnew.ZONES(); [zv.append(z) for z in newzones]
filler=pcbnew.ZONE_FILLER(b); print('fill ok', filler.Fill(zv))
for z in newzones: print(' filled area', round(z.GetFilledArea()/1e12,1),'mm2')
try:
    ok=pcbnew.SaveBoard(OUT,b); print('save',ok)
except Exception as e: print('SAVE ERROR',repr(e))
import os; print('size',os.path.getsize(OUT))
