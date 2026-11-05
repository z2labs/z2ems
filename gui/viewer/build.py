"""Build the standalone viewer page (index.html) from viewer.body.html.

The body file is the page content without the document skeleton, so it can also be
embedded in another page; index.html adds the skeleton for local use.
"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
body = open(os.path.join(here, 'viewer.body.html'), encoding='utf-8').read()
local = '--local' in sys.argv          # test build: three.js from ./three instead of the CDN
if local:
    body = body.replace('https://cdn.jsdelivr.net/npm/three@0.170.0/', './three/')
page = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
        '</head>\n<body>\n' + body + '\n</body>\n</html>\n')
out = os.path.join(here, 'index.local.html' if local else 'index.html')
# Python sources the in-browser KiCad import runs in Pyodide
import shutil
os.makedirs(os.path.join(here, 'py'), exist_ok=True)
csx = os.environ.get('Z2_CSXCAD_PY', os.path.join(here, '..', '..', 'CSXCAD', 'python', 'CSXCAD'))
for src in (os.path.join(csx, 'KiCadImport.py'),
            os.path.join(here, '..', 'tools', 'kicad2view.py'),
            os.path.join(here, '..', '..', 'twin', 'rf_ports.py'),
            os.path.join(here, '..', '..', 'twin', 'kicad_sim.py')):
    if os.path.exists(src):
        shutil.copy(src, os.path.join(here, 'py', os.path.basename(src)))
    else:
        print('missing', src)
open(out, 'w', encoding='utf-8').write(page)
print('wrote', out)
