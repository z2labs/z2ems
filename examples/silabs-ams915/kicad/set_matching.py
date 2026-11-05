"""Set the matching-network values (AN848 figures 2/17/32/47/77/92/107/122, text for WES0115) as footprint Values
so the twin (--components) places them as lumped elements. N.F. = DNP (open pads)."""
import os
# folder holding the Silicon Labs MSC-AMS915 files (download them from silabs.com) and the kicad/ work tree
ROOT = os.environ.get("AMS915_DIR", os.getcwd())
import pcbnew, glob, os
VAL = {'RS1':'0R','CP12':'1.5pF','CP11':'DNP',
       'LS2':'10nH','CP22':'1.2pF','CP21':'DNP',
       'LS3':'22nH','CP32':'0.75pF','CP31':'DNP',
       'LS4':'10nH','CP42':'1pF','CP41':'DNP',
       'RS5':'0R','CP52':'DNP','CP51':'DNP',
       'RS6':'0R','CP62':'4.7pF','CP61':'DNP',
       'LS7':'22nH','CP72':'0.75pF','CP71':'DNP',
       'LS8':'33nH','CP82':'DNP','CP81':'2.7pF',
       'CP91':'3.3pF'}
for f in sorted(glob.glob(os.path.join(ROOT, "kicad/boards/*.kicad_pcb"))):
    b=pcbnew.LoadBoard(f); done=[]
    for fp in b.GetFootprints():
        r=fp.GetReference()
        if r in VAL: fp.SetValue(VAL[r]); done.append('%s=%s'%(r,VAL[r]))
    pcbnew.SaveBoard(f+'.tmp',b); os.replace(f+'.tmp',f)
    print(os.path.basename(f), ' '.join(done))
