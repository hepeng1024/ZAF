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
template is always part of the application data. Each release also needs an
editable copy beside the application:

```text
Linux/Windows: beside the ZAF executable
macOS:         beside ZAF.app in the extracted release folder
```

The runtime reads the editable external copy. If it is missing or invalid, the
GUI warns the user and falls back to its built-in defaults.

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
    └── ZAF_instrument_settings.txt
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

ZAF does not yet have a repeatable Windows packaging helper equivalent to the
Linux and macOS scripts. Until one is added, build natively on Windows with the
same `--onedir` layout.

From PowerShell in the repository root:

```powershell
conda env create -f environment.yml
conda activate zaf

python -m PyInstaller `
  --onedir `
  --windowed `
  --clean `
  --noconfirm `
  --name ZAF `
  --icon "assets\ZAF.ico" `
  --add-data "assets;assets" `
  --add-data "ZAF_instrument_settings.txt;." `
  --collect-data matplotlib `
  --hidden-import matplotlib.backends.backend_tkagg `
  --hidden-import PIL._tkinter_finder `
  --exclude-module pytest `
  --exclude-module tests `
  ZAF_gui.py

Copy-Item ZAF_instrument_settings.txt dist\ZAF\ZAF_instrument_settings.txt
```

Test `dist\ZAF\ZAF.exe` on Windows, including all assets, analysis, simulators,
downloads, and the invalid/missing-settings fallback. Distribute the complete
`dist\ZAF\` folder, not only `ZAF.exe`.

Create a release ZIP after testing:

```powershell
Compress-Archive -Path dist\ZAF -DestinationPath ZAF-Windows-x86_64.zip
```

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
