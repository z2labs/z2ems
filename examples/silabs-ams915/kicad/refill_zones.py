"""Refill all zones with KiCad (the PADS-imported stored fills are the pour outlines WITHOUT pad/track
clearances -> the GND fill covered the SMA centre pin and shorted the EM port)."""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import pcbnew, glob, os
for f in sorted(glob.glob(os.path.join(ROOT, "kicad/boards/*.kicad_pcb"))):
    b=pcbnew.LoadBoard(f)
    ds=b.GetDesignSettings()
    zs=pcbnew.ZONES()
    for z in b.Zones():
        if not z.GetIsRuleArea():
            if z.GetLocalClearance()<int(0.2e6): z.SetLocalClearance(int(0.2e6))
            z.SetMinThickness(int(0.15e6)); z.SetPadConnection(pcbnew.ZONE_CONNECTION_FULL)
            zs.append(z)
    pcbnew.ZONE_FILLER(b).Fill(zs)
    info=[]
    for z in b.Zones():
        if z.GetIsRuleArea(): continue
        L=z.GetFirstLayer(); fp=z.GetFilledPolysList(L)
        info.append('%s/%s %.0fmm2 holes=%d'%(b.GetLayerName(L), z.GetNetname() or '-', fp.Area()/1e12, sum(fp.HoleCount(i) for i in range(fp.OutlineCount()))))
    pcbnew.SaveBoard(f+'.tmp',b); os.replace(f+'.tmp',f)
    print(os.path.basename(f)[:28], '|', '; '.join(info))
