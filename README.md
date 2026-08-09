# ZAF

ZAF is a Python/Tkinter zone-axis finder for FCC, BCC, and HCP TEM
diffraction patterns. It detects diffraction spots, matches them against
structure-specific analytic reference patterns, predicts reachable zone axes
for a double-tilt holder, and includes sample-rotation, tilt, pole-figure,
crystal-lattice, reciprocal-lattice, and diffraction simulators.

The start screen provides separate FCC, BCC, and HCP analysis modes. FCC and
BCC use their respective systematic reflection conditions. HCP uses a
hexagonal direct/reciprocal basis, its two-atom basis reflection condition, and
a configurable c/a ratio.

The easiest way to use ZAF is to download the desktop package for your
operating system from GitHub Releases. Running from source is also supported
and is the best option for users who want to inspect the code or receive the
latest changes with `git pull`.

## Download Desktop Packages

Open the [ZAF Releases page](https://github.com/hepeng1024/zone_axis_finder/releases)
and download the asset for your operating system. Release filenames may
include a version number. Desktop packages include the Python runtime and do
not require a separate Python or Conda installation. The available operating
systems may vary between releases.

Keep the complete extracted folder together. ZAF uses support libraries,
images, and `ZAF_instrument_settings.txt` stored beside the application.

### Windows

1. Download the ZAF Windows `.zip` asset from the Releases page.
2. Before extracting, right-click the ZIP file and select **Properties →
   General → Unblock → Apply** when the Unblock option is available.
3. Extract the complete ZIP file.
4. Open the extracted ZAF folder and double-click `ZAF.exe`.

Do not copy or share only `ZAF.exe`; the packaged support files must remain
beside it.

### Linux

1. Download the ZAF Linux `.tar.gz` asset. The unversioned filename is
   `ZAF-Linux-x86_64.tar.gz`.
2. Extract and enter the application directory:

```bash
tar -xzf ZAF-Linux-x86_64.tar.gz
cd ZAF-Linux-x86_64
```

3. Start ZAF:

```bash
./run_ZAF.sh
```

To add ZAF to the current user's application menu without `sudo`, run:

```bash
./install_launcher.sh
```

If the downloaded archive includes a version in its filename, substitute that
exact filename in the extraction command.

### macOS Apple Silicon

1. Download the ZAF macOS ARM64 `.zip` asset, normally
   `ZAF-macOS-arm64.zip`.
2. Unzip the downloaded file and open its `ZAF-macOS-arm64` folder.
3. Keep `ZAF.app` and `ZAF_instrument_settings.txt` together.
4. Double-click `ZAF.app`.

The current macOS package is intended for Apple-silicon/ARM64 Macs. Because it
is not Developer ID signed or notarized, macOS may block the first launch.
Right-click `ZAF.app`, choose **Open**, and confirm when prompted.

## Run From Source

Install [Git](https://git-scm.com/downloads) and
[Anaconda](https://www.anaconda.com/download) or
[Miniconda](https://docs.conda.io/projects/miniconda/). Then clone the
repository:

```bash
git clone https://github.com/hepeng1024/zone_axis_finder.git
cd zone_axis_finder
```

Create and activate the supplied Conda environment:

```bash
conda env create -f environment.yml
conda activate zaf
```

Start the graphical application:

```bash
python ZAF_gui.py
```

The source environment contains Tkinter, NumPy, SciPy, Pillow, Matplotlib, and
the packaging tools used by this project. If the environment already exists,
update it instead of creating it again.

## Update A Source Checkout

```bash
cd zone_axis_finder
git pull
conda env update -f environment.yml --prune
conda activate zaf
python ZAF_gui.py
```

## Quick Use

1. Choose FCC, BCC, or HCP on the landing page.
2. Select an experimental diffraction image.
3. Enter the holder alpha and beta angles at which the image was recorded.
4. Optionally enter a known current zone axis; leave it blank for automatic
   identification.
5. Select the target zone-axis families.
6. Click **Run Analysis**.
7. Review the fitted diffraction pattern, zone-axis map, reachable targets,
   and sample simulators. Use each image tab's **Download** button to save a
   result.

Use **Crystal Selection** in the analysis window to return to the FCC/BCC/HCP
landing page.

## TEM Instrument Defaults

ZAF reads `ZAF_instrument_settings.txt` once when the GUI starts. Edit the five
documented numeric values in that file to set the startup alpha tilt limits,
beta tilt limits, and image-to-holder rotation for a particular TEM. Restart
ZAF after saving changes. Holder order intentionally remains `xy` and is not
configured by this file.

For a packaged release, keep the editable settings file beside the `ZAF`
executable on Linux/Windows or beside `ZAF.app` on macOS. When running from
source, it stays in the repository root beside `ZAF_gui.py`.

If the settings file is missing, unreadable, or invalid, ZAF displays a warning
and uses the built-in defaults: alpha −35° to 35°, beta −20° to 20°, and
image-to-holder 90°.

## Crystal-System Notes

For HCP, the indexing panel includes:

- four-index Miller–Bravais directions by default, with a toggle for the
  internal three-index representation;
- a c/a entry, whose default is the ideal value of approximately 1.633;
- three-index input such as `[1 0 0]`, or four-index Miller–Bravais input such
  as `[2-1-10]` when four-index mode is enabled.

FCC and BCC use the 26 primitive three-index families with nonnegative h, k,
l and h+k+l ≤ 8. HCP uses the nonduplicated four-index catalog `<0001>`,
`<2-1-10>`, `<10-10>`, `<10-11>`, `<10-12>`, `<11-23>`, `<21-30>`, and
`<40-43>`.

## Command Line

Advanced users running from source can call the calculation backend directly:

```bash
python ZAF.py IMAGE --alpha 0 --beta 0 --crystal-structure FCC
python ZAF.py IMAGE --alpha 0 --beta 0 --crystal-structure BCC
python ZAF.py IMAGE --alpha 0 --beta 0 --crystal-structure HCP \
  --hcp-c-over-a 1.633 --current-zone "2-1-10"
```

Run `python ZAF.py --help` for all matching, calibration, holder, map, and
export options.

## Scale-Bar Calibration

ZAF can detect the bright horizontal scale-bar line and measure its pixel
length. The printed reciprocal-space value, for example 5 nm⁻¹, still needs to
be entered because that physical value is not reliably encoded in the line
itself. Supplying an expected lattice parameter lets ZAF use this calibration
to help rank otherwise ambiguous candidate zones.

## Local Processing And Privacy

ZAF processes diffraction images locally on the user's computer. It does not
upload selected images to a server. Saved result images are written only to
the location selected by the user.

## Developer Documentation

Application-building instructions have been moved out of this user README:

- [Desktop packaging notes](docs/PACKAGING.md)

Ordinary users do not need the packaging tools or developer instructions.

## Main Files

- `ZAF_gui.py`: graphical interface, landing page, and simulators.
- `ZAF.py`: indexing, crystal geometry, matching, plotting, and tilt math.
- `ZAF_instrument_settings.txt`: editable startup TEM tilt/calibration values.
- `environment.yml`: Conda environment definition (`zaf`).
- `requirements.txt`: pinned Python runtime dependencies.
- `assets/`: application icons and GUI artwork.
- `packaging/`: developer packaging and bundle-verification scripts.
- `docs/`: developer documentation.
- `tests/`: local regression and crystal-structure tests when included in the
  checkout.
