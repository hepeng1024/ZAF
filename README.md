# ZAF

ZAF is a Python/Tkinter zone-axis finder for TEM diffraction patterns. It can
analyze its built-in FCC, BCC, and HCP structures or a crystal supplied as a
CIF file. It detects diffraction spots, matches them against calculated
reference patterns, predicts reachable zone axes for a double-tilt holder,
and includes sample-rotation, tilt, pole-figure, crystal-lattice,
reciprocal-lattice, and diffraction simulators.

The start screen provides separate FCC, BCC, and HCP analysis modes plus a
section for your own CIF structures. FCC and BCC use their respective
systematic reflection conditions. HCP uses a hexagonal direct/reciprocal
basis, its two-atom basis reflection condition, and a configurable c/a ratio.

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

Keep the complete extracted folder together. ZAF uses bundled support
libraries and images. The editable instrument-settings location depends on
the operating system, as described under **TEM Instrument Defaults** below.

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
3. Read `README_MACOS.txt` for the first-launch and instrument-settings notes.
4. Double-click `ZAF.app`.

The current macOS package is intended for Apple-silicon/ARM64 Macs. Because it
is not Developer ID signed or notarized, macOS may block the first launch.
Right-click `ZAF.app`, choose **Open**, and confirm when prompted.

On first launch, the app creates its editable instrument settings at
`~/Library/Application Support/ZAF/ZAF_instrument_settings.txt`. This per-user
location is necessary because macOS may run downloaded applications through
App Translocation. The macOS ZIP intentionally does not place another settings
file beside `ZAF.app`, because edits to that copy would be unreliable and
confusing.

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

The source environment contains Tkinter, NumPy, SciPy, Pillow, Matplotlib,
pymatgen for CIF import, and the packaging tools used by this project. If the
environment already exists, update it instead of creating it again.

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

Use **Crystal Selection** in the analysis window to return to the landing
page.

### Your crystal structures

On the landing page, use **Your crystal structures (CIF)** → **Import CIF...**
to select a `.cif` file. ZAF saves its own copy and opens analysis with that
structure selected. For a later session, choose it in the saved-crystal list
and click **Analyze selected**; the original CIF need not be selected again.
**Remove from list** deletes only ZAF's saved copy, not your original file.

The library is private to the current computer user. Its copied CIF files are
in `~/.local/share/ZAF/crystals/` on Linux (or under `$XDG_DATA_HOME/ZAF` if
set), `~/Library/Application Support/ZAF/crystals/` on macOS, and
`%APPDATA%\ZAF\crystals\` on Windows. It is separate from the editable TEM
instrument-settings file.

ZAF reads the full CIF cell without converting it to a primitive cell, so
displayed indices refer to the supplied cell. Hexagonal CIF structures use
the HCP-style four-index direction notation by default, with a three-index
toggle; other custom structures use three-index directions. Automatic search
tests directions with |u|+|v|+|w| ≤ 8 for cubic CIFs, ≤ 6 for tetragonal,
orthorhombic, and trigonal CIFs, or ≤ 4 for monoclinic and triclinic CIFs.
Hexagonal CIFs use the eight HCP-style candidate families. To test another
direction, enter it as the **Known current zone**. This is a finite-candidate,
kinematic geometry match, not a full dynamical TEM simulation. Atom-site
occupancies and atomic numbers are used to screen systematic absences through
an approximate structure-factor phase sum, not to score spot intensities;
experimental intensities can differ substantially. A CIF must contain a
readable cell and occupied atom sites. The detailed tilt simulator remains
available for the built-in structures; the sample-rotation map can be used
with custom structures.

### Multiphase diffraction

If a diffraction image contains more than one phase, select **two or more**
built-in or saved-CIF phases in the landing page's **Multiphase diffraction**
section. Enter the shared image and holder angles, then click **Run Multiphase Analysis**.
ZAF fits a separate reciprocal lattice to each phase and shows a combined
image with a different colored ring for each phase, along with individual
fitted-pattern previews and downloadable images. Rings can overlap where both
phases explain a spot. The results also list reachable target axes for each
phase. A match is supporting evidence, not proof that a phase is present or a
quantitative phase-fraction measurement.

Under **Selected phases**, leave a zone blank to search automatically or enter
a known zone for that phase. The optional lattice parameter `a` is in nm. If
the image has a printed scale bar, enter its labeled value in `1/nm` in the
**Advanced** tab; ZAF measures the line length in pixels. A known lattice
parameter plus that scale can resolve half-spacing aliases caused by spots
from another phase. Without a printed bar, a known lattice value can still
help compare relative spacings when an orientation relationship connects two
phases with known cell dimensions.

An optional orientation relationship accepts **one or two specific, signed
parallel direction pairs** between Phase 1 and Phase 2. Relationship 1 is
required when this option is enabled; Relationship 2 can be blank for both
phases. One pair can constrain its corresponding zone when one fitted phase is
on that specified axis, but it cannot determine rotation around that axis or
predict unrelated zones. ZAF therefore performs one bounded image fit rather
than searching indefinitely over that free rotation. Two pairs define a full
three-dimensional variant.

The relationship does not make every member of two families parallel. For the
D0₂₄ Ni₃Ti example, one complete variant is FCC `[1 1 0]` ∥ Ni₃Ti
`[2 -1 -1 0]` and FCC `[1 -1 1]` ∥ Ni₃Ti `[0 0 0 1]`. With the supplied CIF,
that variant maps FCC `[1 0 1]` close to Ni₃Ti `[2 0 -2 3]`, and FCC `[1 0 0]`
close to Ni₃Ti `[4 -4 0 3]`. Different signed or permuted FCC variants lead
to different predictions. ZAF reports a warning if two entered pairs are
geometrically inconsistent or if the diffraction evidence conflicts with an
orientation-relationship prediction. The multiphase page does not enable the
single-phase sample simulator; inspect each phase's fitted overlay and target
table.

## TEM Instrument Defaults

ZAF reads `ZAF_instrument_settings.txt` once when the GUI starts. Edit the five
documented numeric values in that file to set the startup alpha tilt limits,
beta tilt limits, and image-to-holder rotation for a particular TEM. Restart
ZAF after saving changes. Holder order intentionally remains `xy` and is not
configured by this file.

For packaged Linux and Windows releases, edit the settings file beside the
`ZAF` executable. When running from source, edit the repository-root copy
beside `ZAF_gui.py`.

For packaged macOS releases, edit:

```text
~/Library/Application Support/ZAF/ZAF_instrument_settings.txt
```

ZAF creates this per-user file from its bundled template on first launch and
does not overwrite later user changes. In Finder, choose **Go → Go to Folder**
and enter `~/Library/Application Support/ZAF` to open it. This location remains
stable when Gatekeeper launches the downloaded app from an App Translocation
path or when `ZAF.app` is moved.

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
the location selected by the user. Imported CIFs are also processed locally;
ZAF retains a per-user copy only when you add one to its crystal list.

## Developer Documentation

Application-building instructions have been moved out of this user README:

- [Desktop packaging notes](docs/PACKAGING.md)

Ordinary users do not need the packaging tools or developer instructions.

## Main Files

- `ZAF_gui.py`: graphical interface, landing page, and simulators.
- `ZAF.py`: indexing, crystal geometry, matching, plotting, and tilt math.
- `ZAF_crystals.py`: CIF parsing and the per-user crystal library.
- `ZAF_multiphase.py`: separate fits, shared overlay, and multiphase predictions.
- `ZAF_orientation.py`: signed-variant orientation-relationship geometry.
- `ZAF_instrument_settings.txt`: editable startup TEM tilt/calibration values.
- `environment.yml`: Conda environment definition (`zaf`).
- `requirements.txt`: pinned Python runtime dependencies.
- `assets/`: application icons and GUI artwork.
- `packaging/`: developer packaging and bundle-verification scripts.
- `docs/`: developer documentation.
- `tests/`: local regression and crystal-structure tests when included in the
  checkout.
