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

## Installation

Install [Anaconda](https://www.anaconda.com/download) or
[Miniconda](https://docs.conda.io/projects/miniconda/), then clone or download
the project:

```bash
git clone https://github.com/hepeng1024/ZAF
cd ZAF
conda env create -f environment.yml
conda activate zaf
```

To update an existing environment:

```bash
conda env update -f environment.yml --prune
```

## Run the GUI

```bash
python ZAF_gui.py
```

Select FCC, BCC, or HCP on the landing page, choose a diffraction image, enter
the current holder alpha/beta angles, and click **Run Analysis**. Use
**Crystal Selection** in the analysis window to return to the landing page.

For HCP, the indexing panel includes:

- four-index Miller–Bravais directions by default, with a toggle for the
  internal three-index representation;
- a c/a entry (the default is the ideal value, about 1.633);
- three-index input such as `[1 0 0]`, or four-index Miller–Bravais input
  such as `[2-1-10]` when four-index mode is enabled.

FCC and BCC use the 26 primitive three-index families with nonnegative h, k,
l and h+k+l ≤ 8. HCP instead uses this nonduplicated four-index catalog:
`<0001>`, `<2-1-10>`, `<10-10>`, `<10-11>`, `<10-12>`, `<11-23>`,
`<21-30>`, and `<40-43>`.

## Command line

The core program can also be run directly:

```bash
python ZAF.py IMAGE --alpha 0 --beta 0 --crystal-structure FCC
python ZAF.py IMAGE --alpha 0 --beta 0 --crystal-structure BCC
python ZAF.py IMAGE --alpha 0 --beta 0 --crystal-structure HCP \
  --hcp-c-over-a 1.633 --current-zone "2-1-10"
```

Run `python ZAF.py --help` for all matching, calibration, holder, map, and
export options.

## Scale-bar calibration

ZAF can detect the bright horizontal scale-bar line and measure its pixel
length. The printed reciprocal-space value (for example, 5 nm⁻¹) still needs
to be entered because that physical value is not reliably encoded in the
line itself. Supplying an expected lattice parameter lets ZAF use this
calibration to help rank otherwise ambiguous candidate zones.

## Project files

- `ZAF_gui.py`: graphical interface, landing page, and simulators.
- `ZAF.py`: indexing, crystal geometry, matching, plotting, and tilt math.
- `environment.yml`: Conda environment definition (`zaf`).
- `assets/`: GUI artwork, including the FCC/BCC/HCP landing images.
- `tests/`: regression and crystal-structure tests.

Keep the `assets/` folder beside `ZAF_gui.py`. If Tkinter does not open,
recreate or update the `zaf` environment from `environment.yml`.
