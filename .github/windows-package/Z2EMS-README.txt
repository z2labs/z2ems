Z2 EMS engine for Windows
=========================

Electromagnetic simulation of KiCad boards: openEMS (CPU and CUDA engines) with
a KiCad importer, automatic PCB mesh, lumped ports and S-parameter output.
Nothing to install: unzip and run. https://github.com/z2labs/z2ems


Where to put it
---------------

The KiCad plugin (https://github.com/z2labs/z2ems-kicad) looks for the engine in

    %LOCALAPPDATA%\Z2EMS\engine      (recommended)
    C:\Z2EMS
    %USERPROFILE%\z2ems

i.e. the folder that holds bin\, python\ and twin\. Anywhere else: set the
environment variable Z2EMS_ENGINE to that folder, or write
{"engine": "D:/path/to/Z2EMS"} into %USERPROFILE%\.z2ems.json.


Contents
--------

    bin\               openEMS.exe (CPU engines and --engine=cuda), AppCSXCAD.exe, DLLs,
                       the Octave/Matlab interface and Python wheels of openEMS
    python\            a private Python with CSXCAD, openEMS, numpy, h5py,
                       matplotlib, scipy and shapely
    twin\              KiCad board -> S-parameters (kicad_sim.py), see twin\README.md
    gui\               3D viewer and local simulation server, see gui\README.md
    examples\          test coupons, see examples\README.md
    z2ems-sim.cmd      runs twin\kicad_sim.py with the private Python
    z2ems-server.cmd   runs the local viewer/simulation server


Command line
------------

    z2ems-sim.cmd examples\coupon_msl50.kicad_pcb --port J1 --port J2 --fmax 6e9
    z2ems-sim.cmd board.kicad_pcb --port J1 --port J2 --variant cuda
    z2ems-server.cmd C:\my\boards --open

The coupon gives Z0 about 49.5 ohm and an effective permittivity of about 3.35.


GPU
---

--variant cuda uses the CUDA engine: an NVIDIA GPU of the Pascal generation
(GTX 10xx) or newer with a current driver. --variant auto takes the GPU when
there is one. The engine checks its own GPU memory during the run and repeats a
run that a GPU fault spoiled (Z2_GPU_RETRIES, default 2).


License
-------

GPL-3.0, see LICENSE.txt. openEMS and CSXCAD are (C) Thorsten Liebig and
contributors, https://openEMS.de; documentation: https://docs.openems.de
