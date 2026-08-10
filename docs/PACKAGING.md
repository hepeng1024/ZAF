# ZAF Desktop Packaging Notes

This document is for developers who build and publish ZAF desktop packages.
Ordinary users should follow the installation instructions in the main
[`README.md`](../README.md).

PyInstaller is not a cross-compiler. Build each package on its target operating
system:

```text
Build on Windows → Windows package
Build on Linux   → Linux package
Build on macOS   → macOS package
```

Do not commit generated `build/`, `dist/`, release archives, `*.spec`, or cache
directories. Publish the final archives as GitHub Release assets.

## Required Repository Inputs

Keep these files in the repository for packaging:

```text
ZAF.py
ZAF_gui.py
ZAF_instrument_settings.txt
requirements.txt
environment.yml
assets/
packaging/
.github/workflows/build-macos.yml
```

The operating-system icons are:

```text
Linux:   assets/ZAF.png
Windows: assets/ZAF.ico
macOS:   assets/ZAF.icns
```

`ZAF_gui.py` is the only PyInstaller entry point. `ZAF.py` remains a backend
module imported by the GUI and must not be packaged as a second application.

## Instrument Settings In Packages

PyInstaller bundles an internal copy of `ZAF_instrument_settings.txt` so the
template is always part of the application data. The runtime location is:

```text
Linux/Windows: beside the ZAF executable
macOS:         ~/Library/Application Support/ZAF/ZAF_instrument_settings.txt
```

On macOS, ZAF creates the per-user file from the bundled template on first
launch and preserves it across app upgrades. This is required because an
unsigned or non-notarized downloaded app may run from a temporary Gatekeeper
App Translocation path and cannot see a sibling file beside the original
`ZAF.app`. The macOS release therefore does not include a top-level settings
copy beside the app; the authoritative template remains bundled internally.

On Linux and Windows, the runtime reads the editable external copy beside the
executable. If the active file is missing or invalid, the GUI warns the user
and falls back to its built-in defaults.

## Linux x86_64 Package

The repeatable Linux workflow is designed for native Ubuntu/WSL Ubuntu x86_64
using the `zaf` Conda environment.

Create or update the environment:

```bash
conda env create -f environment.yml
conda activate zaf
```

For an existing environment:

```bash
conda env update -f environment.yml --prune
conda activate zaf
```

Build, verify, and archive the application:

```bash
python packaging/package_linux.py
```

The script verifies that it is using the `zaf` environment's Python and
PyInstaller, runs the source tests, builds a `--onedir`/`--windowed` app,
copies the launchers and editable settings file, checks all bundled ELF files
with `ldd`, runs the noninteractive bundle self-test, performs a GUI smoke test
when a display is available, and verifies a freshly extracted archive.

Generated outputs are separated under:

```text
packaging/build/work/
packaging/build/spec/
packaging/build/dist/
packaging/release/ZAF-Linux-x86_64/
packaging/release/ZAF-Linux-x86_64.tar.gz
```

The release has one top-level folder containing approximately:

```text
ZAF-Linux-x86_64/
├── ZAF
├── _internal/
├── assets/ZAF.png
├── ZAF_instrument_settings.txt
├── run_ZAF.sh
├── install_launcher.sh
├── ZAF.desktop
└── README_LINUX.txt
```

Verify an existing release directory manually with:

```bash
bash packaging/verify_linux_bundle.sh \
  packaging/release/ZAF-Linux-x86_64
```

Set `ZAF_REQUIRE_GUI_SMOKE=1` when a graphical display is available and the
verification must open the bundled GUI:

```bash
ZAF_REQUIRE_GUI_SMOKE=1 bash packaging/verify_linux_bundle.sh \
  packaging/release/ZAF-Linux-x86_64
```

Linux portability is limited by the build host's ABI. A package built on one
Ubuntu version may require the same or a newer compatible glibc on another
machine. Test the archive on at least one clean machine outside the build
environment before publishing it broadly.

## macOS Apple-Silicon Package

The macOS ARM64 package can be built with GitHub Actions or manually on an
Apple-silicon Mac. Do not attempt this build from Linux or Windows.

### GitHub Actions Build

The workflow is:

```text
.github/workflows/build-macos.yml
```

It runs manually through `workflow_dispatch` and automatically for tags that
match `v*`. On GitHub:

```text
Repository → Actions → Build macOS ZAF → Run workflow
```

After a successful run, download the `ZAF-macOS-arm64` artifact. It contains:

```text
ZAF-macOS-arm64.zip
└── ZAF-macOS-arm64/
    ├── ZAF.app
    └── README_MACOS.txt
```

The workflow uses a native `macos-15` ARM64 runner and Python 3.11. It checks
the architecture, repository inputs, source imports, bundled resources,
PyInstaller output, Mach-O slices, app metadata, ad-hoc code signature, and
the extracted archive. Calculation tests run when tracked `tests/test_*.py`
files are present; otherwise the workflow reports a warning and continues.

### Manual macOS Build

On an Apple-silicon Mac:

```bash
git clone https://github.com/hepeng1024/zone_axis_finder.git
cd zone_axis_finder
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r packaging/requirements-macos-build.txt
bash packaging/package_macos.sh
```

To validate repository inputs without building:

```bash
bash packaging/package_macos.sh --check-only
```

The output archive is:

```text
packaging/release/macos/ZAF-macOS-arm64.zip
```

The current app is ad-hoc signed for testing and institutional sharing, but it
is not Developer ID signed or notarized. A public polished release would need
an Apple Developer identity, hardened-runtime signing, and notarization.

## Windows x86_64 Package

The repeatable Windows workflow must run on native 64-bit Windows. PyInstaller
is not a cross-compiler, so do not run this build in WSL or with a WSL/Linux
Python interpreter. The scripts explicitly use the native Conda environment
named `zaf`; they do not depend on the environment that launched VS Code or the
current PowerShell session.

Create the environment if it does not exist:

```powershell
conda env create -f environment.yml
conda run -n zaf python -m pip install PyInstaller
```

For an existing environment, update it and make sure PyInstaller is installed
inside `zaf`, not in `base` or globally:

```powershell
conda env update -n zaf -f environment.yml --prune
conda run -n zaf python -m pip install --upgrade PyInstaller
```

Verify the interpreter before packaging:

```powershell
conda run -n zaf python -c "import os, sys; print(sys.executable); print(sys.platform); print(os.environ.get('CONDA_DEFAULT_ENV'))"
conda run -n zaf python -m PyInstaller --version
```

The first command must report a native Windows interpreter under the `zaf`
environment, `win32`, and `zaf`. A parent shell may still report that `base` is
active; that is harmless because every Python command in the workflow is
executed through `conda run -n zaf`.

Build, verify, smoke-test, and archive the application from PowerShell in the
repository root:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File .\packaging\windows\package_windows.ps1
```

The script validates the source inputs, icon, runtime assets, imports, and
instrument settings before building. It runs the available project tests and
the source self-tests, then invokes the environment-specific PyInstaller as
equivalent to:

```powershell
conda run -n zaf python -m PyInstaller `
  --noconfirm `
  --clean `
  --onedir `
  --windowed `
  --noupx `
  --name ZAF `
  --icon <repository>\assets\ZAF.ico `
  --add-data "<repository>\assets;assets" `
  --add-data "<repository>\ZAF_instrument_settings.txt;." `
  ZAF_gui.py
```

`ZAF_gui.py` is the only entry script. PyInstaller discovers the imported
`ZAF.py` backend and uses its maintained hooks for NumPy, SciPy, Pillow,
Matplotlib, and Tkinter. Do not add broad hidden imports or manually copied DLLs
unless a clean build and bundle test demonstrate a specific missing component.

The two settings-file copies serve different purposes:

```text
ZAF-Windows-x86_64/ZAF_instrument_settings.txt
    Editable file read by the frozen GUI at startup.

ZAF-Windows-x86_64/_internal/ZAF_instrument_settings.txt
    Bundled template used by the noninteractive bundle self-test.
```

PyInstaller creates the internal copy through `--add-data`. The packaging
script then copies the repository template beside `ZAF.exe` as an ordinary,
writable file. Both copies must initially match, but the external copy remains
independently editable by the user.

Generated outputs are separated under:

```text
build_work/windows/work/
build_work/windows/spec/
build_work/windows/staging/
dist/windows/ZAF/
release/ZAF-Windows-x86_64/
release/ZAF-Windows-x86_64.zip
```

Only those ZAF-specific generated paths are cleaned. The release ZIP has one
top-level folder containing approximately:

```text
ZAF-Windows-x86_64/
|-- ZAF.exe
|-- ZAF_instrument_settings.txt
|-- README_WINDOWS.txt
`-- _internal/
    |-- ZAF_instrument_settings.txt
    |-- assets/
    `-- bundled Python libraries, Tcl/Tk files, .pyd modules, and DLLs
```

Distribute the complete folder or ZIP, never `ZAF.exe` alone. The executable
requires `_internal`, and the editable settings file must remain beside it.

The packaging script calls `packaging/windows/verify_windows_bundle.ps1`
before creating the archive, after creating it, and again after extracting it
to a disposable path containing spaces. It checks the PE32+ x86_64 GUI
executable and embedded icon, assets, backend module, scientific libraries,
Tcl/Tk runtime, internal and external settings copies, development-file
exclusions, complete ZIP contents, `--bundle-self-test`, and
`--gui-smoke-test`. It also exercises edited, missing, and invalid external
settings in the disposable extracted copy and restores the default afterward.

Verify an existing release directory and archive manually with:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File .\packaging\windows\verify_windows_bundle.ps1 `
  -BundlePath .\release\ZAF-Windows-x86_64 `
  -ArchivePath .\release\ZAF-Windows-x86_64.zip
```

Before publishing, manually open the freshly extracted application and test
the FCC, BCC, and HCP modes, analysis and preview images, both simulators,
saved output, window resizing, and edited instrument defaults. Test on another
Windows computer when practical. The executable is currently unsigned, so a
SmartScreen reputation warning is possible; do not disable Windows security
protections, and distinguish such a warning from an actual malware detection.

## Release Checklist

1. Update any intended version/tag and user-facing release notes.
2. Confirm the source GUI starts and the calculation tests pass when present.
3. Build Windows on Windows, Linux on Linux, and macOS on Apple silicon.
4. Test FCC, BCC, and HCP landing images and analysis flows.
5. Test the fitted pattern, zone-axis map, sample-rotation simulator, tilt
   simulator, image downloads, and editable instrument defaults.
6. Test each archive after extraction, preferably on a clean machine.
7. Upload the operating-system-specific archives to GitHub Releases.
8. Do not commit generated build directories, release archives, or spec files.
