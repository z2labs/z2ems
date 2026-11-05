#!/usr/bin/env python3
""" z2mon: watch the Z2 EM runs of a folder in the browser, while they run and after.

        python z2mon.py <runs root> [--port 8766] [--host 0.0.0.0] [--open]       (--open starts the default browser)

    Every run that make_panel.sim_coupon / calib_sim / antenna_cases produce (a folder with model.xml, sim/stdout.txt,
    run.json) shows up in the Runs tab (the evolution: id, label, mesh, |S11| max, time, status). Selecting a run gives
    Geometry (top view, launch zoom with cuts, 3D), Progress (energy / speed while it runs, with an unstable-run
    warning), S-parameters (dB, Smith, group delay), Fields (current-density maps per copper layer and frequency,
    when the run was made with Z2_DUMPS=1) and Far field (the pattern plots of antenna runs). Pictures are made by
    runviz.py / fields.py / plot_pattern.py and cached next to the run. No internet needed (works on the offline
    laptop); the same pictures feed the KiCad plugin's result view.
"""
import glob, json, os, re, sys, threading, time, webbrowser
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs, quote

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import runviz

ROOT = None
LOCK = threading.Lock()


_CACHE = {}


def find_runs(root):
    """ run folders: anything with model.xml and sim/ (the twin layout <out>/sim/cNN/excE/) or antenna cases (<id>/sim/);
        entries are cached by the mtime of model.xml / stdout.txt / run.json, so thousands of old runs cost nothing """
    runs = []
    for xml in glob.glob(os.path.join(root, '**', 'model.xml'), recursive=True):
        d = os.path.dirname(xml)
        if os.path.basename(d) == 'sim':                 # antenna_cases: <case>/sim/model.xml (the twin copies model.xml into sim/ too)
            d = os.path.dirname(d)
            if os.path.exists(os.path.join(d, 'model.xml')):
                continue
            if not os.path.exists(os.path.join(d, 'sim', 'stdout.txt')) and not os.path.exists(os.path.join(d, 'result.json')):
                continue
        sim = os.path.join(d, 'sim')
        rj = os.path.join(d, 'run.json')
        stdout_p = os.path.join(sim, 'stdout.txt')
        key = (os.path.getmtime(xml), os.path.getmtime(stdout_p) if os.path.exists(stdout_p) else 0, os.path.getmtime(rj) if os.path.exists(rj) else 0)
        c = _CACHE.get(d)
        if c and c[0] == key:
            runs.append(c[1])
            continue
        info = {}
        if os.path.exists(rj):
            try:
                info = json.load(open(rj))
            except Exception:
                info = {}
        if not info:                                     # older twin runs: <out>/sim/cNN.json next to the run folder
            m = re.match(r'c(\d+)$', os.path.basename(os.path.dirname(d)))
            cj = os.path.join(os.path.dirname(os.path.dirname(d)), os.path.basename(os.path.dirname(d)) + '.json') if m else os.path.join(d, 'result.json')
            if os.path.exists(cj):
                try:
                    c = json.load(open(cj))
                    info = dict(id=c.get('id'), label=c.get('label') or c.get('id'), mcells=c.get('mcells'), mesh=c.get('mesh'), wall_s=c.get('wall_s'),
                                s11_max_3G=c.get('s11_max_3G'), s11_max_6G=c.get('s11_max_6G'), s21_min_6G=c.get('s21_min_6G'), status='done')
                except Exception:
                    pass
        rel = os.path.relpath(d, root)
        st = info.get('status')
        stdout = os.path.join(sim, 'stdout.txt')
        if not st:
            st = 'done' if os.path.exists(stdout) and 'Time for' in open(stdout, errors='replace').read()[-3000:] else ('running' if os.path.exists(stdout) else 'queued')
        mt = max(os.path.getmtime(xml), os.path.getmtime(stdout) if os.path.exists(stdout) else 0)
        t0 = None                                        # start: run.json "started", else when model.xml was written
        for fmt_ in ('%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H:%M', '%Y-%m-%dT%H:%M:%S'):
            try:
                t0 = time.mktime(time.strptime(str(info.get('started') or '')[:19], fmt_))
                break
            except ValueError:
                pass
        if t0 is None:
            t0 = os.path.getmtime(xml)
        t1 = None                                        # end: start + wall time, else the last write of the output
        if st not in ('running', 'queued'):
            try:
                t1 = t0 + float(info['wall_s']) if info.get('wall_s') else mt
            except (TypeError, ValueError):
                t1 = mt
        entry = dict(dir=rel, name=info.get('label') or rel, id=info.get('id'), status=st, mtime=mt, t_start=t0, t_end=t1,
                         mcells=info.get('mcells'), mesh=info.get('mesh'), wall_s=info.get('wall_s'),
                         s11_3=info.get('s11_max_3G'), s11_6=info.get('s11_max_6G'), s21_6=info.get('s21_min_6G'),
                         started=info.get('started'), batch=os.path.relpath(os.path.dirname(os.path.dirname(d)) if os.path.basename(os.path.dirname(d)).startswith('c') else os.path.dirname(d), root))
        _CACHE[d] = (key, entry)
        runs.append(entry)
    runs.sort(key=lambda r: (r.get('t_start') or r['mtime'], r['mtime']))   # by start time (the evolution chart); the table shows it reversed
    return runs


_3D_CLS = ['cu', 'pin', 'body', 'via']


def geom3d(root, rel):
    """ metal faces of model.xml for the small 3D viewer: faces = [[class, x0,y0,z0, x1,y1,z1, ...], ...] (mm),
        lo/hi of the metal, the dielectric box; cached in <run>/viz/geom3d.json by the mtime of model.xml """
    d = os.path.join(root, rel)
    if not os.path.abspath(d).startswith(os.path.abspath(root)):
        return {'error': 'outside root'}
    xml = os.path.join(d, 'model.xml')
    if not os.path.exists(xml):
        xml = os.path.join(d, 'sim', 'model.xml')
    if not os.path.exists(xml):
        return {'error': 'no model.xml yet'}
    cache = os.path.join(d, 'viz')
    os.makedirs(cache, exist_ok=True)
    cj = os.path.join(cache, 'geom3d.json')
    if os.path.exists(cj) and os.path.getmtime(cj) >= os.path.getmtime(xml):
        try:
            return json.load(open(cj))
        except Exception:
            pass
    np = runviz.np
    prims, _mesh = runviz.load_model(xml)
    faces = []
    lo = np.array([np.inf] * 3); hi = -lo
    dlo = lo.copy(); dhi = hi.copy()

    def add(ci, pts):
        pts = np.asarray(pts, dtype=float)
        if len(pts) < 3:
            return
        a = np.cross(pts[1] - pts[0], pts[-1] - pts[0])
        if len(pts) == 4 and float(np.dot(a, a)) < 1e-14:
            return                                              # degenerate face of a zero-thickness box
        nonlocal lo, hi
        lo = np.minimum(lo, pts.min(0)); hi = np.maximum(hi, pts.max(0))
        faces.append([ci] + [round(float(v), 4) for v in pts.ravel()])

    for p in sorted(prims, key=lambda q: q['pri']):
        col, al, cls = runviz.style(p)
        if cls in ('diel', 'ptfe'):
            if p['kind'] == 'box':
                dlo = np.minimum(dlo, p['lo']); dhi = np.maximum(dhi, p['hi'])
            elif p['kind'] == 'poly' and p['axis'] == 2:
                dlo = np.minimum(dlo, [p['pts'][:, 0].min(), p['pts'][:, 1].min(), min(p['z0'], p['z1'])])
                dhi = np.maximum(dhi, [p['pts'][:, 0].max(), p['pts'][:, 1].max(), max(p['z0'], p['z1'])])
            continue
        if cls not in _3D_CLS:
            continue
        ci = _3D_CLS.index(cls)
        if p['kind'] == 'box':
            (x0, y0, z0), (x1, y1, z1) = p['lo'], p['hi']
            v = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0], [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
            for f in ([0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4]):
                add(ci, v[f])
        elif p['kind'] == 'poly' and p['axis'] == 2:
            for z in ((p['z0'], p['z1']) if abs(p['z1'] - p['z0']) > 1e-9 else (p['z0'],)):
                add(ci, np.column_stack([p['pts'], np.full(len(p['pts']), z)]))
        elif p['kind'] in ('cyl', 'shell'):
            n = 12; a = np.linspace(0, 2 * np.pi, n, endpoint=False)
            s, e, r = p['s'], p['e'], p['r']
            if abs(s[2] - e[2]) > 1e-9 and abs(s[0] - e[0]) < 1e-9 and abs(s[1] - e[1]) < 1e-9:      # vertical
                r0 = np.column_stack([s[0] + r * np.cos(a), s[1] + r * np.sin(a), np.full(n, s[2])]); r1 = r0.copy(); r1[:, 2] = e[2]
            elif abs(s[1] - e[1]) < 1e-9 and abs(s[2] - e[2]) < 1e-9:                                   # along x
                r0 = np.column_stack([np.full(n, s[0]), s[1] + r * np.cos(a), s[2] + r * np.sin(a)]); r1 = r0.copy(); r1[:, 0] = e[0]
            else:                                                                                        # along y
                r0 = np.column_stack([s[0] + r * np.cos(a), np.full(n, s[1]), s[2] + r * np.sin(a)]); r1 = r0.copy(); r1[:, 1] = e[1]
            for i in range(n):
                add(ci, [r0[i], r0[(i + 1) % n], r1[(i + 1) % n], r1[i]])
            add(ci, r0); add(ci, r1)
    if not faces:
        return {'error': 'no metal in model.xml'}
    out = dict(faces=faces, cls=_3D_CLS, lo=[float(x) for x in lo], hi=[float(x) for x in hi], zs=2.0,
               diel=[[float(x) for x in dlo], [float(x) for x in dhi]] if np.isfinite(dlo).all() else None)
    try:
        json.dump(out, open(cj, 'w'), separators=(',', ':'))
    except OSError:
        pass
    return out


def run_detail(root, rel):
    d = os.path.join(root, rel)
    if not os.path.abspath(d).startswith(os.path.abspath(root)):
        return {'error': 'outside root'}
    sim = os.path.join(d, 'sim')
    xml = os.path.join(d, 'model.xml')
    cache = os.path.join(d, 'viz')
    os.makedirs(cache, exist_ok=True)
    out = dict(dir=rel)
    try:
        out['run'] = json.load(open(os.path.join(d, 'run.json')))
    except Exception:
        out['run'] = {}
    # geometry (cached by model.xml mtime)
    gj = os.path.join(cache, 'geom.json')
    if os.path.exists(xml) and (not os.path.exists(gj) or os.path.getmtime(gj) < os.path.getmtime(xml)):
        with LOCK:
            try:
                g = runviz.geometry(xml, os.path.join(cache, 'g'))
                json.dump(g, open(gj, 'w'))
            except Exception as e:
                json.dump({'error': str(e)}, open(gj, 'w'))
    if os.path.exists(gj):
        out['geometry'] = json.load(open(gj))
    # progress (always fresh)
    stdout = os.path.join(sim, 'stdout.txt')
    if os.path.exists(stdout):
        with LOCK:
            out['progress'] = runviz.progress_png(stdout, os.path.join(cache, 'progress.png'))
        out['progress']['png'] = os.path.join(cache, 'progress.png')
        out['progress']['tail'] = open(stdout, errors='replace').read()[-1500:]
    # S-parameters: the touchstone of the coupon (sibling of the run folder) or the case
    s2p = out['run'].get('s2p')
    if not s2p:
        cand = glob.glob(os.path.join(os.path.dirname(d), '*.s2p')) + glob.glob(os.path.join(os.path.dirname(os.path.dirname(d)), '*.s2p')) + \
            glob.glob(os.path.join(d, '*.s2p')) + glob.glob(os.path.join(d, '*.s1p'))
        m = re.match(r'(c\d+)', os.path.basename(os.path.dirname(d)) if os.path.basename(d).startswith('exc') else os.path.basename(d))
        cand = [c for c in cand if not m or m.group(1) in os.path.basename(c)] or cand
        s2p = cand[0] if cand else None
    if s2p and os.path.exists(s2p):
        sj = os.path.join(cache, 'sparams.json')
        if not os.path.exists(sj) or os.path.getmtime(sj) < os.path.getmtime(s2p):
            with LOCK:
                try:
                    json.dump(runviz.sparams(s2p, os.path.join(cache, 's'), title=out['run'].get('label', '')), open(sj, 'w'))
                except Exception as e:
                    json.dump({'error': str(e)}, open(sj, 'w'))
        out['sparams'] = json.load(open(sj))
    # fields (made by make_panel when Z2_DUMPS=1)
    if out['run'].get('fields'):
        out['fields'] = out['run']['fields']
        out['fields']['base'] = os.path.join(d, 'fields')
    # far field (antenna_cases + plot_pattern)
    pats = sorted(glob.glob(os.path.join(d, '*_pattern2d.png')) + glob.glob(os.path.join(d, '*_pattern3d.png')))
    if pats:
        out['farfield'] = pats
        pj = glob.glob(os.path.join(d, '*_pattern.json'))
        if pj:
            try:
                out['farfield_metrics'] = json.load(open(pj[0]))
            except Exception:
                pass
    rj = os.path.join(d, 'result.json')
    if os.path.exists(rj):
        try:
            r = json.load(open(rj))
            out['result'] = {k: v for k, v in r.items() if not isinstance(v, list)}
        except Exception:
            pass
    return out


PAGE = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Z2 EM run monitor</title>
<style>
body{font-family:system-ui,Segoe UI,Arial,sans-serif;margin:0;background:#f4f5f7;color:#222}
header{background:#1f2d3d;color:#fff;padding:8px 14px;display:flex;gap:16px;align-items:center}
header b{font-size:15px} header span{opacity:.8;font-size:12px}
nav{display:flex;gap:4px;padding:6px 10px;background:#e4e7ec;border-bottom:1px solid #ccd}
nav button{border:1px solid #bbc;background:#fff;padding:6px 14px;border-radius:5px;cursor:pointer;font-size:13px}
nav button.on{background:#1f5fa8;color:#fff;border-color:#1f5fa8}
main{padding:10px 14px} table{border-collapse:collapse;font-size:12.5px;background:#fff} th,td{border:1px solid #dde;padding:3px 8px;text-align:right}
th{background:#eef1f5} td:nth-child(2),td:nth-child(3){text-align:left} tr.sel{background:#dbe7f7} tr:hover{background:#eef4fb;cursor:pointer}
.ok{color:#1a7f3c} .run{color:#b36b00;font-weight:600} .bad{color:#c8231d;font-weight:600}
img{max-width:100%;border:1px solid #ddd;background:#fff;margin:4px 0} .row{display:flex;flex-wrap:wrap;gap:10px}
.row>div{flex:1 1 480px} pre{background:#111;color:#cfe;font-size:11px;padding:8px;max-height:260px;overflow:auto}
select,label{font-size:13px} .menu{margin:6px 0;display:flex;gap:12px;align-items:center;flex-wrap:wrap}
canvas{background:#fff;border:1px solid #ddd} .kv{font-size:12.5px;line-height:1.5} .warn{background:#fde8e6;color:#8a1c14;padding:6px 10px;border-radius:5px;font-weight:600}
.top{display:flex;gap:10px;align-items:stretch;margin-bottom:8px} #evo{flex:1 1 auto;min-width:0;width:100%;height:250px;cursor:pointer}
#v3w{flex:0 0 360px;position:relative;background:#fff;border:1px solid #ddd} #v3d{display:block;width:100%;height:250px;border:0;cursor:grab;touch-action:none}
#v3cap{position:absolute;left:7px;top:5px;right:7px;font-size:11px;line-height:1.35;color:#333;pointer-events:none}
#v3hint{position:absolute;right:7px;bottom:4px;font-size:10px;color:#999;pointer-events:none}
td.t{white-space:nowrap;font-variant-numeric:tabular-nums}
@media(max-width:900px){.top{flex-direction:column} #v3w{flex-basis:auto}}
</style></head><body>
<header><b>Z2 EM run monitor</b><span id="root"></span><span id="clock"></span>
<label style="margin-left:auto;font-size:12px"><input type="checkbox" id="follow" checked> follow the latest run</label>
<label style="font-size:12px"><input type="checkbox" id="auto" checked> auto-refresh 5 s</label></header>
<nav><button data-t="runs" class="on">Runs / evolution</button><button data-t="geom">Geometry</button><button data-t="prog">Progress</button>
<button data-t="spar">S-parameters</button><button data-t="field">Fields / currents</button><button data-t="ff">Far field</button></nav>
<main>
<section id="t-runs"><div class="menu"><label>batch <select id="batch"></select></label><span id="nrun"></span></div>
<div class="top"><canvas id="evo" width="1100" height="250"></canvas>
<div id="v3w"><canvas id="v3d" width="360" height="250"></canvas><div id="v3cap"></div><div id="v3hint">drag: rotate · wheel: zoom · double-click: reset</div></div></div>
<div id="runs"></div></section>
<section id="t-geom" hidden><div class="kv" id="gmesh"></div><div class="row"><div><img id="g_launch"></div><div><img id="g_3d"></div></div><img id="g_top"></section>
<section id="t-prog" hidden><div id="pwarn"></div><div class="kv" id="pkv"></div><img id="p_png"><pre id="ptail"></pre></section>
<section id="t-spar" hidden><div class="kv" id="skv"></div><img id="s_s"><div class="row"><div><img id="s_smith"></div><div><img id="s_gd"></div></div></section>
<section id="t-field" hidden><div class="menu"><label>layer <select id="flayer"></select></label><label>f <select id="ffreq"></select></label><span id="fpeak"></span></div><img id="f_img"><div class="kv" id="fnote"></div></section>
<section id="t-ff" hidden><div id="ffimgs"></div><pre id="ffm"></pre></section>
</main>
<script>
const $=id=>document.getElementById(id); let RUNS=[], SEL=null, DETAIL=null, TAB='runs', BATCH='';
const F=p=>'/file?path='+encodeURIComponent(p);
document.querySelectorAll('nav button').forEach(b=>b.onclick=()=>{TAB=b.dataset.t;document.querySelectorAll('nav button').forEach(x=>x.classList.toggle('on',x===b));
 ['runs','geom','prog','spar','field','ff'].forEach(t=>$('t-'+t).hidden=(t!==TAB));});
async function refresh(){
 const r=await (await fetch('/api/runs')).json(); RUNS=r.runs; $('root').textContent=r.root; $('clock').textContent=new Date().toLocaleTimeString();
 const batches=[...new Set(RUNS.map(x=>x.batch))]; const sel=$('batch'); const cur=sel.value||BATCH;
 sel.innerHTML='<option value="">all</option>'+batches.map(b=>`<option ${b===cur?'selected':''}>${b}</option>`).join(''); BATCH=sel.value;
 const rows=RUNS.filter(x=>!BATCH||x.batch===BATCH);
 $('nrun').textContent=rows.length+' runs, '+rows.filter(x=>x.status==='running').length+' running';
 if($('follow').checked){const last=rows.filter(x=>x.status==='running').slice(-1)[0]||rows.slice(-1)[0]; if(last&&last.dir!==SEL){SEL=last.dir;}}
 $('runs').innerHTML='<table><tr><th>#</th><th>run (newest on top)</th><th>batch</th><th>status</th><th>cells [M]</th><th>mesh</th><th>|S11|max&lt;3G</th><th>&lt;6G</th><th>S21min&lt;6G</th><th>time [s]</th><th>start</th><th>end</th></tr>'+
  rows.slice().reverse().map((x,j)=>{const i=rows.length-1-j;return `<tr class="${x.dir===SEL?'sel':''}" data-d="${x.dir}"><td>${i+1}</td><td>${x.id!=null?x.id+' ':''}${x.name}</td><td>${x.batch}</td><td class="${x.status==='done'?'ok':x.status==='running'?'run':'bad'}">${x.status}</td><td>${x.mcells?x.mcells.toFixed(2):''}</td><td>${x.mesh?x.mesh.join('×'):''}</td><td>${fmt(x.s11_3)}</td><td>${fmt(x.s11_6)}</td><td>${fmt(x.s21_6,2)}</td><td>${x.wall_s?x.wall_s.toFixed(0):''}</td><td class="t">${fT(x.t_start)}</td><td class="t">${x.status==='running'?'<span class="run">running</span>':fT(x.t_end,x.t_start)}</td></tr>`}).join('')+'</table>';
 document.querySelectorAll('#runs tr[data-d]').forEach(tr=>tr.onclick=()=>{SEL=tr.dataset.d;$('follow').checked=false;load();});
 evo(rows); if(SEL) load();
}
function fmt(v,d){return (v==null)?'':(v>0?'<span class="bad">+'+v.toFixed(d||1)+' ⚠</span>':v.toFixed(d||1));}
const p2=n=>String(n).padStart(2,'0');
function fT(t,ref){if(!t)return'';const d=new Date(t*1000),hm=p2(d.getHours())+':'+p2(d.getMinutes());
 if(ref&&new Date(ref*1000).toDateString()===d.toDateString())return hm;return d.getFullYear()+'-'+p2(d.getMonth()+1)+'-'+p2(d.getDate())+' '+hm;}
const esc=s=>String(s==null?'':s).replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
let EVO=[]; const EV={L:48,R:14,T:20,B:42};
function evo(rows){const c=$('evo');const cw=c.clientWidth||1100;if(c.width!==cw)c.width=cw;const g=c.getContext('2d');g.clearRect(0,0,c.width,c.height);EVO=rows;
 if(!rows.length){g.fillStyle='#888';g.font='12px sans-serif';g.fillText('no runs',EV.L,40);return;}
 const v=rows.map(x=>x.s11_6).filter(x=>x!=null),lo=Math.min(-30,...v),hi=Math.max(0,...v);
 const W=c.width-EV.L-EV.R,H=c.height-EV.T-EV.B,X=i=>EV.L+(rows.length>1?i*W/(rows.length-1):W/2),Y=y=>EV.T+H-(y-lo)/(hi-lo||1)*H;
 g.font='11px sans-serif';g.lineWidth=1;
 for(let y=Math.ceil(lo/10)*10;y<=hi;y+=10){g.strokeStyle='#e4e4e4';g.beginPath();g.moveTo(EV.L,Y(y));g.lineTo(c.width-EV.R,Y(y));g.stroke();g.fillStyle='#555';g.fillText(y+' dB',4,Y(y)+4);}
 // time axis: dashed separator + date label where the day changes, then evenly spaced date/time ticks
 const D=t=>new Date((t||0)*1000);let lastDay=null,lastLab=-1e9;g.textAlign='left';
 rows.forEach((x,i)=>{const d=D(x.t_start),k=d.toDateString();if(lastDay!==null&&k!==lastDay){const xs=(X(i-1)+X(i))/2;g.strokeStyle='#9aa';g.setLineDash([4,3]);g.beginPath();g.moveTo(xs,EV.T);g.lineTo(xs,EV.T+H);g.stroke();g.setLineDash([]);
  if(xs-lastLab>78&&xs<c.width-EV.R-60){g.fillStyle='#678';g.fillText(p2(d.getMonth()+1)+'-'+p2(d.getDate())+' →',xs+3,EV.T+11);lastLab=xs;}}lastDay=k;});
 g.strokeStyle='#888';g.beginPath();g.moveTo(EV.L,EV.T+H);g.lineTo(c.width-EV.R,EV.T+H);g.stroke();
 g.textAlign='center';g.fillStyle='#444';const nT=Math.max(2,Math.floor(W/105));let px=-1e9;
 for(let k=0;k<nT;k++){const i=Math.round(k*(rows.length-1)/Math.max(nT-1,1)),xp=X(i);if(xp-px<80)continue;px=xp;const d=D(rows[i].t_start);
  g.strokeStyle='#888';g.beginPath();g.moveTo(xp,EV.T+H);g.lineTo(xp,EV.T+H+4);g.stroke();
  g.fillText(p2(d.getHours())+':'+p2(d.getMinutes()),xp,EV.T+H+16);g.fillText(p2(d.getMonth()+1)+'-'+p2(d.getDate()),xp,EV.T+H+29);}
 g.textAlign='left';
 g.strokeStyle='#1f5fa8';g.lineWidth=1.5;g.beginPath();let first=true;rows.forEach((x,i)=>{if(x.s11_6==null)return;const a=X(i),b=Y(Math.min(x.s11_6,hi));first?g.moveTo(a,b):g.lineTo(a,b);first=false;});g.stroke();
 rows.forEach((x,i)=>{const sel=x.dir===SEL;if(x.s11_6==null){if(x.status==='running'){g.fillStyle='#b36b00';g.fillRect(X(i)-3,EV.T+H-7,6,6);}return;}
  g.fillStyle=x.s11_6>0?'#c8231d':(sel?'#c8531d':'#1f5fa8');g.beginPath();g.arc(X(i),Y(Math.min(x.s11_6,hi)),sel?5:3,0,7);g.fill();});
 const si=rows.findIndex(x=>x.dir===SEL);if(si>=0){g.strokeStyle='rgba(200,83,29,.35)';g.beginPath();g.moveTo(X(si),EV.T);g.lineTo(X(si),EV.T+H);g.stroke();}
 g.fillStyle='#333';g.fillText('max |S11| 0–6 GHz per run, in start order (x: start date/time) — red = unstable (>0 dB), orange square = running; click a point to open it',EV.L+2,12);}
function evoIdx(e){const c=$('evo'),r=c.getBoundingClientRect(),x=(e.clientX-r.left)*c.width/r.width,W=c.width-EV.L-EV.R;
 return EVO.length?Math.max(0,Math.min(EVO.length-1,Math.round((x-EV.L)/(W/Math.max(EVO.length-1,1))))):-1;}
$('evo').onclick=e=>{const i=evoIdx(e);if(i<0)return;SEL=EVO[i].dir;$('follow').checked=false;evo(EVO);load();};
$('evo').onmousemove=e=>{const i=evoIdx(e);if(i<0)return;const x=EVO[i];$('evo').title=`${x.id!=null?x.id+' ':''}${x.name}\n${fT(x.t_start)} → ${x.status==='running'?'running':fT(x.t_end,x.t_start)}\n|S11|max<6G: ${x.s11_6!=null?x.s11_6.toFixed(1)+' dB':'–'}`;};

// ---- small 3D viewer of the selected (followed) run: metal of model.xml, software-rendered on a 2D canvas (no libraries, works offline)
const V3={d:null,dir:null,yaw:-0.75,pitch:0.62,zoom:1,drag:null,raf:0};
const V3COL=[[230,168,23],[214,39,40],[138,138,138],[70,70,70]];
async function load3d(){if(!SEL)return;if(V3.dir===SEL&&V3.d)return;V3.dir=SEL;V3.d=null;draw3d();const want=SEL;
 let d;try{d=await (await fetch('/api/geom3d?dir='+encodeURIComponent(want))).json();}catch(e){d={error:String(e)};}
 if(V3.dir===want){V3.d=d;draw3d();}}
function req3d(){if(!V3.raf)V3.raf=requestAnimationFrame(()=>{V3.raf=0;draw3d();});}
function draw3d(){const c=$('v3d'),w=c.clientWidth||360;if(c.width!==w)c.width=w;const g=c.getContext('2d');g.fillStyle='#fff';g.fillRect(0,0,c.width,c.height);const d=V3.d;
 if(!d||!d.faces){g.fillStyle='#888';g.font='12px sans-serif';g.textAlign='center';g.fillText(d?(d.error||'no geometry'):(SEL?'loading geometry…':'no run selected'),c.width/2,c.height/2+10);g.textAlign='left';return;}
 const lo=d.lo,hi=d.hi,zs=d.zs||2,cx=(lo[0]+hi[0])/2,cy=(lo[1]+hi[1])/2,cz=(lo[2]+hi[2])/2;
 const TOP=44,BOT=18,ca=Math.cos(V3.yaw),sa=Math.sin(V3.yaw),cp=Math.cos(V3.pitch),sp=Math.sin(V3.pitch);
 const R=(x,y,z,a,b,p,q)=>{x-=cx;y-=cy;z=(z-cz)*zs;const x1=x*a-y*b,y1=x*b+y*a;return[x1,-(y1*q+z*p),y1*p-z*q];};
 if(!d._s0||d._w!==c.width){let a0=1e9,a1=-1e9,b0=1e9,b1=-1e9;const A=Math.cos(-0.75),B=Math.sin(-0.75),Pc=Math.cos(0.62),Ps=Math.sin(0.62);   // fit once, at the default view
  for(const X of[lo[0],hi[0]])for(const Y of[lo[1],hi[1]])for(const Z of[lo[2],hi[2]]){const r=R(X,Y,Z,A,B,Pc,Ps);a0=Math.min(a0,r[0]);a1=Math.max(a1,r[0]);b0=Math.min(b0,r[1]);b1=Math.max(b1,r[1]);}
  d._s0=0.92*Math.min((c.width-16)/Math.max(2*Math.max(-a0,a1),1e-6),(c.height-TOP-BOT)/Math.max(2*Math.max(-b0,b1),1e-6));d._w=c.width;}
 const S=V3.zoom*d._s0,ox=c.width/2,oy=TOP+(c.height-TOP-BOT)/2;
 const P=(x,y,z)=>{const r=R(x,y,z,ca,sa,cp,sp);return[ox+r[0]*S,oy+r[1]*S,r[2]];};
 const L=[];for(const f of d.faces){const n=(f.length-1)/3,q=new Array(n);let dep=0;for(let k=0;k<n;k++){q[k]=P(f[1+3*k],f[2+3*k],f[3+3*k]);dep+=q[k][2];}L.push([dep/n,f[0],q]);}
 L.sort((a,b)=>b[0]-a[0]);const edge=L.length<6000;
 // dielectric outline (behind)
 const box=(a,b,col)=>{const v=[[a[0],a[1],a[2]],[b[0],a[1],a[2]],[b[0],b[1],a[2]],[a[0],b[1],a[2]],[a[0],a[1],b[2]],[b[0],a[1],b[2]],[b[0],b[1],b[2]],[a[0],b[1],b[2]]].map(p=>P(...p));
  g.strokeStyle=col;g.lineWidth=1;g.beginPath();for(const[i,j]of[[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]]){g.moveTo(v[i][0],v[i][1]);g.lineTo(v[j][0],v[j][1]);}g.stroke();};
 if(d.diel)box(d.diel[0],d.diel[1],'rgba(44,160,44,.45)');
 g.lineWidth=0.3;g.strokeStyle='rgba(0,0,0,.28)';
 for(const[dep,ci,q]of L){const a=q[0],b=q[1],m=q[q.length>3?q.length-1:2];
  const ux=b[0]-a[0],uy=b[1]-a[1],uz=(b[2]-a[2])*S,vx=m[0]-a[0],vy=m[1]-a[1],vz=(m[2]-a[2])*S;
  const nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx,nl=Math.hypot(nx,ny,nz)||1;
  const sh=0.45+0.55*Math.abs((0.35*nx-0.55*ny+0.75*nz)/nl),k=V3COL[ci]||V3COL[0];
  g.fillStyle=`rgb(${k[0]*sh|0},${k[1]*sh|0},${k[2]*sh|0})`;g.beginPath();g.moveTo(q[0][0],q[0][1]);for(let i=1;i<q.length;i++)g.lineTo(q[i][0],q[i][1]);g.closePath();g.fill();if(edge)g.stroke();}
 // axis triad + scale
 const t0=[30,c.height-26],ax=[[1,0,0,'x','#c8231d'],[0,1,0,'y','#1a7f3c'],[0,0,1,'z','#1f5fa8']];g.lineWidth=1.5;g.font='10px sans-serif';
 for(const[x,y,z,nm,col]of ax){const x1=x*ca-y*sa,y1=x*sa+y*ca,px=x1*16,py=-(y1*sp+z*cp)*16;g.strokeStyle=col;g.fillStyle=col;g.beginPath();g.moveTo(t0[0],t0[1]);g.lineTo(t0[0]+px,t0[1]+py);g.stroke();g.fillText(nm,t0[0]+px*1.25-3,t0[1]+py*1.25+3);}
 $('v3hint').textContent=`${(hi[0]-lo[0]).toFixed(1)} × ${(hi[1]-lo[1]).toFixed(1)} mm, z ×${zs} · drag: rotate · wheel: zoom · dbl-click: reset`;}
(()=>{const c=$('v3d'),pt=e=>e.touches?[e.touches[0].clientX,e.touches[0].clientY]:[e.clientX,e.clientY];
 const down=e=>{V3.drag=pt(e);c.style.cursor='grabbing';},up=()=>{V3.drag=null;c.style.cursor='grab';};
 const move=e=>{if(!V3.drag)return;const p=pt(e);V3.yaw+=(p[0]-V3.drag[0])*0.01;V3.pitch=Math.max(-1.5,Math.min(1.5,V3.pitch+(p[1]-V3.drag[1])*0.01));V3.drag=p;if(e.touches)e.preventDefault();req3d();};
 c.addEventListener('mousedown',down);window.addEventListener('mouseup',up);window.addEventListener('mousemove',move);
 c.addEventListener('touchstart',down,{passive:true});c.addEventListener('touchend',up);c.addEventListener('touchmove',move,{passive:false});
 c.addEventListener('wheel',e=>{e.preventDefault();V3.zoom=Math.max(0.2,Math.min(40,V3.zoom*(e.deltaY<0?1.15:1/1.15)));req3d();},{passive:false});
 c.addEventListener('dblclick',()=>{V3.yaw=-0.75;V3.pitch=0.62;V3.zoom=1;req3d();});
 window.addEventListener('resize',()=>{evo(EVO);req3d();});})();
function cap3d(){const x=RUNS.find(r=>r.dir===SEL);if(!x){$('v3cap').innerHTML='';return;}const p=(DETAIL&&DETAIL.dir===SEL&&DETAIL.progress)||{};
 const last=p.ts&&p.ts.length?p.ts[p.ts.length-1]:null,pct=(last!=null&&p.max_ts)?Math.min(100,Math.round(100*last/p.max_ts)):null;
 const st=x.status==='running'?`<span class="run">running</span>${pct!=null?' '+pct+' %':''}${p.energy_dB&&p.energy_dB.length?' · energy '+p.energy_dB[p.energy_dB.length-1]+' dB':''}`:`<span class="${x.status==='done'?'ok':'bad'}">${esc(x.status)}</span>`;
 $('v3cap').innerHTML=`<b>${x.id!=null?esc(x.id)+' ':''}${esc(x.name)}</b><br>${esc(x.batch)} · ${fT(x.t_start)} · ${st}${p.diverging?' · <span class="bad">UNSTABLE</span>':''}`;}
async function load(){ if(!SEL)return; load3d(); cap3d(); DETAIL=await (await fetch('/api/run?dir='+encodeURIComponent(SEL))).json(); const d=DETAIL,r=d.run||{}; cap3d();
 const g=d.geometry||{}; $('gmesh').innerHTML=`<b>${r.id!=null?r.id+' ':''}${r.label||SEL}</b> — ${g.mesh?('mesh '+g.mesh.n.join('×')+' = '+g.mesh.cells.toFixed(2)+' M cells, smallest cell '+g.mesh.min_mm.map(x=>x?x.toFixed(3):'?').join(' / ')+' mm'):''} ${g.error||''}`;
 $('g_launch').src=g.launch?F(g.launch):''; $('g_3d').src=g.iso?F(g.iso):''; $('g_top').src=g.top?F(g.top):'';
 const p=d.progress||{}; $('pwarn').innerHTML=p.diverging?'<div class="warn">ENERGY IS GROWING — this run is numerically unstable (check the boundary / margins / geometry)</div>':'';
 $('pkv').innerHTML=p.ts?`timestep ${p.ts.slice(-1)[0]}${p.max_ts?' / '+p.max_ts:''}, ${p.t_s.slice(-1)[0]} s, ${p.speed.slice(-1)[0]} M cells/s, energy ${p.energy_dB.slice(-1)[0]} dB (dt ${p.dt_s?(p.dt_s*1e15).toFixed(2)+' fs':'?'}) — ${p.done?'finished':'running'}`:'no engine output yet';
 $('p_png').src=p.png?F(p.png)+'&t='+Date.now():''; $('ptail').textContent=p.tail||'';
 const s=d.sparams||{}; $('skv').innerHTML=s.s?`|S11| max ${s.s11_max_dB.toFixed(1)} dB${s.s21_min_dB!=null?', |S21| min '+s.s21_min_dB.toFixed(2)+' dB':''} (${s.f0_GHz}–${s.f1_GHz} GHz)`:'no S-parameters yet (the run is not finished)';
 $('s_s').src=s.s?F(s.s):''; $('s_smith').src=s.smith?F(s.smith):''; $('s_gd').src=s.gd?F(s.gd):'';
 const f=d.fields; const L=$('flayer'),Q=$('ffreq'); if(f&&f.maps){const layers=Object.keys(f.maps); if(L.options.length!==layers.length){L.innerHTML=layers.map(l=>`<option>${l}</option>`).join('');} if(Q.options.length!==f.f_GHz.length){Q.innerHTML=f.f_GHz.map((x,i)=>`<option value="${i}">${x} GHz</option>`).join('');}
  const li=L.value||layers[0],qi=+Q.value||0; $('f_img').src=F(f.base+'/'+f.maps[li][qi].split('/').pop()); $('fpeak').textContent='peak '+f.peak_db[li][qi]+' dB (0 dB = strongest current of any layer at this f, '+f.db_range+' dB range)'; $('fnote').textContent='current density |J| (rot H on the copper plane) — frequency-domain dump; run with Z2_DUMPS=1 Z2_DUMP_F=1,3,6';}
 else{$('f_img').src='';$('fnote').textContent='no field dumps for this run (set Z2_DUMPS=1 for the next one)';}
 $('ffimgs').innerHTML=(d.farfield||[]).map(x=>`<img src="${F(x)}">`).join('')||'<i>no far-field data (antenna runs only)</i>'; $('ffm').textContent=d.farfield_metrics?JSON.stringify(d.farfield_metrics,null,1).slice(0,3000):'';
 document.querySelectorAll('#runs tr[data-d]').forEach(tr=>tr.classList.toggle('sel',tr.dataset.d===SEL));
}
$('flayer').onchange=$('ffreq').onchange=()=>load(); $('batch').onchange=()=>{BATCH=$('batch').value;refresh();};
refresh(); setInterval(()=>{if($('auto').checked)refresh();},5000);
</script></body></html>'''


class H(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _json(self, obj, code=200):
        b = json.dumps(obj).encode()
        self.send_response(code); self.send_header('Content-Type', 'application/json'); self.send_header('Content-Length', str(len(b))); self.end_headers(); self.wfile.write(b)

    def do_GET(self):
        u = urlparse(self.path); q = parse_qs(u.query)
        try:
            if u.path == '/' or u.path == '/index.html':
                b = PAGE.encode()
                self.send_response(200); self.send_header('Content-Type', 'text/html; charset=utf-8'); self.send_header('Content-Length', str(len(b))); self.end_headers(); self.wfile.write(b)
            elif u.path == '/api/runs':
                self._json(dict(root=ROOT, runs=find_runs(ROOT)))
            elif u.path == '/api/run':
                self._json(run_detail(ROOT, q.get('dir', [''])[0]))
            elif u.path == '/api/geom3d':
                self._json(geom3d(ROOT, q.get('dir', [''])[0]))
            elif u.path == '/file':
                p = os.path.abspath(q.get('path', [''])[0])
                if not p.startswith(os.path.abspath(ROOT)) or not os.path.isfile(p):
                    self.send_error(404); return
                b = open(p, 'rb').read()
                ct = 'image/png' if p.endswith('.png') else 'text/plain'
                self.send_response(200); self.send_header('Content-Type', ct); self.send_header('Content-Length', str(len(b))); self.send_header('Cache-Control', 'no-cache'); self.end_headers(); self.wfile.write(b)
            else:
                self.send_error(404)
        except Exception as e:
            self._json({'error': str(e)}, 500)


def main():
    global ROOT
    a = sys.argv[1:]
    if not a:
        raise SystemExit(__doc__)
    ROOT = os.path.abspath(a[0])
    port = int(a[a.index('--port') + 1]) if '--port' in a else 8766
    host = a[a.index('--host') + 1] if '--host' in a else '127.0.0.1'       # --host 0.0.0.0 for the ZeroTier / LAN address (read-only pages)
    srv = ThreadingHTTPServer((host, port), H)
    url = 'http://127.0.0.1:%d/' % port
    print('z2mon: %s  ->  %s' % (ROOT, url), flush=True)
    if '--open' in a:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    srv.serve_forever()


if __name__ == '__main__':
    main()
