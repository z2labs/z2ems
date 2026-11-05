# Bundle a Windows build of openEMS into a self-contained folder (and zip): the executables and DLLs of the install
# prefix, every DLL they import from the vcpkg bin folder (found with tools/pe_deps.py), and the C++ runtime of the
# toolset that built them (msvcp140.dll, vcruntime140*.dll) as app-local DLLs.
#
# The CRT copy matters: a binary built with MSVC 14.40 or later (VS 2022 17.10+) needs an msvcp140.dll of at least
# that version; a machine with an older VC++ Redistributable faults at the first std::mutex lock, which in openEMS is
# the multithreaded engine's barrier (single-threaded engines keep working, so it looks like a CPU-specific crash).
# App-local DLLs next to the executable take precedence over System32, so the bundle runs without installing anything.
#
#   powershell -File scripts\win_bundle.ps1 -Prefix <install prefix> -VcpkgBin <vcpkg>\installed\x64-windows-release\bin -Out <dist folder>
#   optional: -RedistDir <...\VC\Redist\MSVC\14.44.35112\x64\Microsoft.VC143.CRT>   (default: $env:VCToolsRedistDir from a VS developer prompt)
param(
    [Parameter(Mandatory = $true)][string]$Prefix,
    [Parameter(Mandatory = $true)][string]$VcpkgBin,
    [Parameter(Mandatory = $true)][string]$Out,
    [string]$RedistDir = '',
    [switch]$NoZip
)
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$pedeps = Join-Path (Split-Path -Parent $here) 'tools\pe_deps.py'
$bin = Join-Path $Prefix 'bin'
if (-not (Test-Path $bin)) { throw "no bin folder under $Prefix" }
New-Item -ItemType Directory -Force $Out | Out-Null
Copy-Item (Join-Path $bin '*') $Out -Force
# dependencies of every executable and DLL in bin
$needed = @{}
foreach ($f in Get-ChildItem $bin -Include *.exe, *.dll -Recurse) {
    foreach ($l in (python $pedeps $f.FullName $bin $VcpkgBin)) {
        $p = ($l -split ' ', 2)[1]
        if ($p -ne 'SYSTEM' -and $p.StartsWith($VcpkgBin, [System.StringComparison]::OrdinalIgnoreCase)) { $needed[$p] = 1 }
    }
}
foreach ($p in $needed.Keys) { Copy-Item $p $Out -Force }
Write-Host ("{0} vcpkg DLLs copied" -f $needed.Count)
# C++ runtime of the building toolset
if (-not $RedistDir) {
    if ($env:VCToolsRedistDir) { $RedistDir = Join-Path $env:VCToolsRedistDir 'x64\Microsoft.VC143.CRT' }
}
if ($RedistDir -and (Test-Path $RedistDir)) {
    foreach ($n in 'msvcp140.dll', 'msvcp140_1.dll', 'msvcp140_2.dll', 'msvcp140_atomic_wait.dll', 'msvcp140_codecvt_ids.dll', 'vcruntime140.dll', 'vcruntime140_1.dll', 'vcruntime140_threads.dll', 'concrt140.dll') {
        $src = Join-Path $RedistDir $n
        if (Test-Path $src) { Copy-Item $src $Out -Force }
    }
    Write-Host ("CRT {0} from {1}" -f (Get-Item (Join-Path $Out 'msvcp140.dll')).VersionInfo.FileVersion, $RedistDir)
} else {
    Write-Warning "no VC++ redist folder found (run from a VS developer prompt or pass -RedistDir); the bundle needs msvcp140.dll >= the build toolset on the target"
}
# manifest
$lines = foreach ($f in Get-ChildItem $Out -File | Sort-Object Name) {
    $v = $f.VersionInfo.FileVersion
    '{0,-45} {1,12:N0}  {2}' -f $f.Name, $f.Length, $v
}
$lines | Set-Content (Join-Path $Out 'BUNDLE_MANIFEST.txt')
if (-not $NoZip) {
    $zip = "$Out.zip"
    if (Test-Path $zip) { Remove-Item $zip -Force }
    Compress-Archive -Path (Join-Path $Out '*') -DestinationPath $zip
    Write-Host ("{0:N0} MB  {1}" -f ((Get-Item $zip).Length / 1MB), $zip)
}
