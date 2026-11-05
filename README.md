# Z2 EMS

Electromagnetic simulation of printed circuit boards, straight from KiCad.

Z2 EMS is a distribution of [openEMS](https://openEMS.de) (the free FDTD field solver by Thorsten Liebig and
contributors) with the parts a PCB designer needs on top:

- **KiCad importer** (`CSXCAD.KiCadImport`): copper of every layer, plated holes, board outline and stackup are
  read straight from a `.kicad_pcb` file; two-pad R/L/C footprints become lumped elements. Fab stackup presets
  (`CSXCAD.FabStackups`, e.g. JLCPCB JLC04161H-7628) and laminate dispersion / copper roughness models
  (`CSXCAD.PCBMaterials`).
- **Board to S-parameters in one command** (`twin/kicad_sim.py`): lumped ports at the connector pads, an
  automatic mesh for PCB geometry (thirds rule on strip edges, diagonal tracks, via antipads, graded cells),
  one run per port, Touchstone output, Z0 and effective permittivity of a line, a mesh DRC that finds copper
  shorted or broken by the staircase mesh before a run.
- **CUDA engine** (`--engine=cuda`): the FDTD update, the absorbing boundaries, excitation, lumped RLC and
  dispersive materials on an NVIDIA GPU, with a memory guard that catches corrupted GPU memory instead of
  returning a wrong result.
- **Conformal PEC** (contour-path FDTD after Yu and Mittra, 2001): a Yee edge partly covered by metal
  keeps the free fraction of its length, so metal edges that do not fall on a mesh line are not snapped
  to it. On the 50 Ω test coupon this moves Z0 from 50.37 Ω (staircase) to 49.50 Ω, against 49.35 Ω
  from the closed-form formula. On by default in the twin, `Z2_CONFORMAL=0` turns it off.
- **Viewer and local server** (`gui/`): a browser 3D view of the board (exploded layers, components) that
  can also start runs and show S-parameters, field maps and far field.
- **KiCad plugin**: [z2labs/z2ems-kicad](https://github.com/z2labs/z2ems-kicad) drives all of the above from
  the PCB editor (ports, stackup, mesh preview, simulate, field overlays).

## Repository layout

```
openEMS/        submodule: z2labs/z2ems-openEMS  (openEMS + CUDA engine, conformal PEC)
CSXCAD/         submodule: z2labs/z2ems-CSXCAD   (CSXCAD + KiCad importer, fab stackups, PCB materials)
fparser, QCSXCAD, AppCSXCAD, CTB   submodules: upstream (thliebig)
twin/           KiCad board -> openEMS model -> S-parameters (see twin/README.md)
gui/            3D viewer and local simulation server (see gui/README.md)
examples/       test coupons and a worked antenna example (see examples/README.md)
scripts/        dependency installer, Python build, Windows bundling
doc-src/        openEMS documentation sources
```

## Install

**Windows:** download the engine package from the
[releases](https://github.com/z2labs/z2ems/releases) and unzip it, e.g. to `%LOCALAPPDATA%\Z2EMS\engine` (where
the KiCad plugin looks first). It contains the openEMS executables (CPU and CUDA engine), a Python with the
CSXCAD/openEMS bindings and all packages the twin needs, `twin/`, `gui/` and `examples/`. The CUDA engine needs
an NVIDIA driver; without one the CPU engine is used.

**Linux, from source:**

```bash
git clone --recursive https://github.com/z2labs/z2ems.git
cd z2ems
./scripts/install_deps.sh          # prints the package install command for your distribution
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/opt/z2ems -DBUILD_APPCSXCAD=NO \
         -DWITH_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES="61;75;86;89"   # leave out both for a CPU-only build
make -j$(nproc)
cd .. && ./scripts/build_python.sh --cpp-install-dir $HOME/opt/z2ems   # CSXCAD + openEMS Python bindings
pip install numpy h5py matplotlib scipy shapely
export Z2_OPENEMS=$HOME/opt/z2ems/bin/openEMS
```

The CUDA engine needs a CUDA 12.x toolkit for Pascal GPUs (GTX 10xx, `61`); CUDA 13 dropped them.

## First run

```bash
cd twin
python kicad_sim.py ../examples/coupon_msl50.kicad_pcb --port J1 --port J2 --fmax 6e9            # CPU
python kicad_sim.py ../examples/coupon_msl50.kicad_pcb --port J1 --port J2 --fmax 6e9 --variant cuda
```

Expected: Z0 ≈ 49.5 Ω, εeff ≈ 3.35, |S11| below −28 dB up to 6 GHz. `twin/README.md` has the options, the
model rules, the environment switches and the validation cases.

## Relation to openEMS

The engine is openEMS with additions; models, scripts and post-processing written for openEMS run unchanged,
and the openEMS documentation at [docs.openEMS.de](https://docs.openems.de) applies.

## License

GPL-3.0, see [LICENSE](LICENSE); CSXCAD, QCSXCAD and fparser are LGPL-3.0, CTB is BSD-style. openEMS and CSXCAD
are © Thorsten Liebig and contributors; the Z2 EMS additions are © Zoltan Doczi / Z2Labs. Contributions follow [AI_POLICY.md](AI_POLICY.md) (disclosure and DCO sign-off).
