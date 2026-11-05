"""Give the no-net copper polygons (printed antennas, feed patches) the SMA signal net so the twin meshes them
(it builds mesh lines only from port-net copper) and treats them as connected copper."""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import pcbnew, glob, os
for f in sorted(glob.glob(os.path.join(ROOT, "kicad/boards/*.kicad_pcb"))):
    b=pcbnew.LoadBoard(f)
    sma=[fp for fp in b.GetFootprints() if fp.GetReference().startswith('SMA')][0]
    sig=[p for p in sma.Pads() if p.GetNetname()!='GND'][0]; nc=sig.GetNetCode()
    n=0
    for d in b.GetDrawings():
        if d.GetClass()=='PCB_SHAPE' and pcbnew.IsCopperLayer(d.GetLayer()) and d.GetNetCode()==0:
            d.SetNetCode(nc); n+=1
    pcbnew.SaveBoard(f+'.tmp',b); os.replace(f+'.tmp',f)
    print(os.path.basename(f)[:28], 'net', sig.GetNetname(), 'shapes re-netted', n)
