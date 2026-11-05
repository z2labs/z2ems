import re, sys, math
from PIL import Image, ImageDraw
def render(path, out, px_per_mm=10, crop=None):
    txt=open(path,errors='replace').read()
    fmt=re.search(r'%FSLAX(\d)(\d)Y(\d)(\d)',txt); xd=int(fmt.group(2))
    unit=25.4 if '%MOIN' in txt else 1.0
    sc=unit/10**xd
    aps={}
    for m in re.finditer(r'%ADD(\d+)([CRO]),([\d.]+)(?:X([\d.]+))?',txt):
        aps[int(m.group(1))]=(m.group(2),float(m.group(3))*unit,float(m.group(4))*unit if m.group(4) else None)
    for m in re.finditer(r'%ADD(\d+)(\w+),([\d.]+)',txt):
        if int(m.group(1)) not in aps: aps[int(m.group(1))]=('C',float(m.group(3))*unit,None)
    body=re.sub(r'%[^%]*%','',txt).replace('\r','').replace('\n','')
    cmds=[c for c in body.split('*') if c]
    x=y=0.0; ap=None; mode=1; pts=[]; prims=[]; quad=False; region=None; lastop=2
    xmin=ymin=1e9; xmax=ymax=-1e9
    for c in cmds:
        if c.startswith('G04'): continue
        m=re.match(r'^G54D(\d+)$',c) or re.match(r'^D(\d+)$',c)
        if m and int(m.group(1))>=10: ap=int(m.group(1)); continue
        if c=='G75': quad=True; continue
        if c.startswith('G36'): region=[]; continue
        if c.startswith('G37'):
            if region and len(region)>2: prims.append(('P',region))
            region=None; continue
        g=re.search(r'G0([123])',c)
        if g: mode=int(g.group(1))
        mx=re.search(r'X(-?\d+)',c); my=re.search(r'Y(-?\d+)',c); mi=re.search(r'I(-?\d+)',c); mj=re.search(r'J(-?\d+)',c)
        d=re.search(r'D0([123])$',c)
        if not d and (mx or my) and lastop==1: d=re.match(r'(1)','1')
        if d: lastop=int(d.group(1))
        nx=float(mx.group(1))*sc if mx else x; ny=float(my.group(1))*sc if my else y
        if d and region is not None:
            op=int(d.group(1))
            if op==2 and region: prims.append(('P',region)); region=[]
            region.append((nx,ny)); x,y=nx,ny; continue
        if d:
            op=int(d.group(1))
            if op==1 and ap in aps:
                if mode==1 or not (mi or mj):
                    prims.append(('L',x,y,nx,ny,aps[ap]))
                else:
                    i=float(mi.group(1))*sc if mi else 0; j=float(mj.group(1))*sc if mj else 0
                    cx,cy=x+i,y+j; r=math.hypot(i,j)
                    a0=math.atan2(y-cy,x-cx); a1=math.atan2(ny-cy,nx-cx)
                    if mode==2:   # CW
                        while a1>=a0: a1-=2*math.pi
                    else:
                        while a1<=a0: a1+=2*math.pi
                    n=max(4,int(abs(a1-a0)*r*2)+2); px,py=x,y
                    for k in range(1,n+1):
                        a=a0+(a1-a0)*k/n; qx,qy=cx+r*math.cos(a),cy+r*math.sin(a)
                        prims.append(('L',px,py,qx,qy,aps[ap])); px,py=qx,qy
            elif op==3 and ap in aps: prims.append(('F',nx,ny,aps[ap]))
            for px_,py_ in ((nx,ny),(x,y)): xmin=min(xmin,px_);xmax=max(xmax,px_);ymin=min(ymin,py_);ymax=max(ymax,py_)
        x,y=nx,ny
    if crop: xmin,ymin,xmax,ymax=crop
    W=int((xmax-xmin)*px_per_mm)+20; H=int((ymax-ymin)*px_per_mm)+20
    im=Image.new('L',(W,H),255); dr=ImageDraw.Draw(im)
    T=lambda X,Y:((X-xmin)*px_per_mm+10, H-((Y-ymin)*px_per_mm+10))
    for p in prims:
        if p[0]=='P':
            dr.polygon([T(X,Y) for X,Y in p[1]],fill=0); continue
        if p[0]=='L':
            _,x0,y0,x1,y1,a=p; w=max(1,int(a[1]*px_per_mm)); dr.line([T(x0,y0),T(x1,y1)],fill=0,width=w)
            r=w/2
            for X,Y in ((x0,y0),(x1,y1)):
                u,v=T(X,Y); dr.ellipse([u-r,v-r,u+r,v+r],fill=0)
        else:
            _,X,Y,a=p; u,v=T(X,Y); r=a[1]*px_per_mm/2
            if a[0]=='R': rh=(a[2] or a[1])*px_per_mm/2; dr.rectangle([u-r,v-rh,u+r,v+rh],fill=0)
            else: dr.ellipse([u-r,v-r,u+r,v+r],fill=0)
    im.save(out); print(out, 'extent mm: %.1f..%.1f x %.1f..%.1f'%(xmin,xmax,ymin,ymax), len(prims),'prims')
if __name__=='__main__':
    render(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv)>3 else 10,
           tuple(map(float,sys.argv[4].split(','))) if len(sys.argv)>4 else None)
