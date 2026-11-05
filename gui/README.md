# Z2 EMS viewer and server

A browser 3D view of a KiCad board for the visual check before simulation:
exploded copper, dielectric and mask layers with their true thickness ratios,
plated hole barrels, and component bodies that can be hidden per side or one
by one.

- `tools/kicad2view.py`: KiCad board to viewer data (`.z2v.json`). Uses
  `CSXCAD.KiCadImport` (no KiCad needed) and `shapely`. Per copper layer the
  shapes of each net are merged, holes cut out and the result triangulated.
- `viewer/viewer.body.html`: the viewer (three.js WebGPURenderer, WebGL 2
  fallback). `viewer/build.py` wraps it into the standalone `viewer/index.html`.

## Use

    python tools/kicad2view.py path/to/board.kicad_pcb viewer/board.z2v.json
    python tools/kicad2view.py board.kicad_pcb out.z2v.json --fab JLCPCB:JLC04161H-7628
    cd viewer && python -m http.server 8000      # open http://localhost:8000/

`index.html` loads `board.z2v.json` next to it, `?board=<file>` loads another
file, and any `.z2v.json` can be dropped onto the view. Over `file://` the
browser blocks the fetch; use the local server.

## Local twin (simulate from the viewer)

    python gui/server/z2twin_server.py <board folder> --open

serves the viewer and a JSON API on 127.0.0.1 (only files below the given
folder). In this mode the viewer

- lists the `.kicad_pcb` files of the folder and imports them with the same
  Python importer (no Pyodide);
- suggests the RF ports (connector footprints, net names, net classes) and
  snaps hand-placed ports to the nearest outer-layer pad; ports are numbered in
  list order, the arrow points from the reference copper to the signal pad, and
  each port has its own reference impedance Z0;
- offers the fab stackups that fit the board (selecting one reloads the board
  with it), laminate dispersion (Djordjevic-Sarkar via Debye poles), copper
  loss with roughness, and the fitted R/L/C footprints as lumped elements;
- runs openEMS (one excitation per port, queued, cancellable) with live
  progress, and shows |S| in dB, the Smith chart, the TDR impedance profile
  and the phase; the Touchstone file downloads from the result;
- maps the current density of every copper layer (port-1 excitation, three
  frequencies) onto the layers of the 3D view: hide layers to see the return
  current on the inner planes.

Results are stored next to the board in `<board>_twin/<date-time>/`
(`result.json`, `.sNp`, `fields/*.png`) and listed under “Earlier runs”.
Without the server (e.g. the published page) the viewer is a viewer: a
`result.json` can still be dropped onto it to see the plots.

The engine: `Z2_OPENEMS`, else `bin/openEMS` of the engine package, else `openEMS` on PATH (see
`twin/README.md`).

## Remote use: a strong PC runs, a laptop or phone looks

The engine runs where the server runs; the browser only shows the board and the results, so a
small laptop or a phone is enough on the viewing side.

- **Same network:** `python gui/server/z2twin_server.py <board folder> --lan --threads 16`
  prints one URL per network address, each ending in `?token=...`. Open it once on the laptop or
  phone (scan it from the terminal, or send it to yourself); the token then lives in a cookie for
  30 days. Without the token the server answers 401. `--token` (or `Z2TWIN_TOKEN`) keeps the same
  token across restarts.
- **Tailscale (recommended away from home):** keep the server on this computer only (no `--lan`)
  and publish it in the tailnet: `tailscale serve --bg 8765`. Every device logged in to the tailnet
  opens `https://<machine>.<tailnet>.ts.net/` with HTTPS and no token; nothing is exposed to the
  internet. Do not use `tailscale funnel` for client boards: that makes the page public.
- A GPU build of openEMS (`--engine=cuda`, see openEMS `-DWITH_CUDA=ON`) is picked automatically
  (`Engine: Auto`) when an NVIDIA driver is installed. Tested on Linux with Pascal GPUs (GTX 1080 / 1080 Ti);
  the Windows CUDA build of the release package is experimental.

## In the viewer

- Separation and Z scale sliders: Z scale multiplies every layer thickness by
  the same factor, so the ratios stay true. Components are drawn at true height.
- Click a copper shape or a barrel: its net is highlighted on all layers
  (Esc clears). Click a dielectric or mask: thickness, εr, tan δ, material.
- Click a component, then H or Delete hides it; "Show hidden" brings them back.
- Double-click a layer in the stackup list to show it alone.
- Components with a VRML 3D model (`.wrl`, also found next to a `.step` reference) are drawn
  with it. **These models are illustration only**: the view marks them with "⚠ 3D models:
  illustration only, not simulated".
- The simulated bodies (connector pin, flange, legs, coax feed; metal cyan, dielectric violet)
  are read from the openEMS model file `<board>_sim.xml` next to the board (or
  `kicad2view.py --sim model.xml`); board copper and laminate are not repeated from it. When
  both exist, the Components section switches between "⚠ 3D model", "Simulated" and "Both".
- When exploded, what is above the board stays rigid on the top layer, what is below rigid under
  the bottom layer, and the slice within the board thickness (an edge connector's flange and
  barrel band) stretches with the layers, so legs never cut a layer. Separation 0 and Z scale 1x
  show the true geometry.

## Limits of this version

- Component heights come from the 3D model bounding box (VRML/STEP, KiCad's
  model directories and the board folder are searched); without a model file
  they are estimated from the footprint name and marked as estimated. STEP-only models are
  drawn as these boxes (no STEP tessellation yet).
- Copper text and silkscreen are not shown; mask openings are not cut.
- Large boards: the Jetson Nano baseboard (8 layers) exports to 22.6 MB and
  loads in about 3 s; a binary format comes later.
