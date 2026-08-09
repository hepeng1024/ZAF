#!/usr/bin/env python3
"""Build, verify, and archive the Linux ZAF application bundle."""

from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import NoReturn


APP_NAME = "ZAF"
ARCHIVE_STEM = "ZAF-Linux-x86_64"
REQUIRED_PROJECT_PATHS = (
    "ZAF.py",
    "ZAF_gui.py",
    "ZAF_instrument_settings.txt",
    "assets",
    "assets/ZAF.png",
)
REQUIRED_IMPORTS = (
    "tkinter",
    "numpy",
    "scipy",
    "PIL",
    "matplotlib",
    "ZAF",
    "ZAF_gui",
    "PyInstaller",
)


def fail(message: str) -> NoReturn:
    raise SystemExit(f"Linux packaging failed: {message}")


def run(command: list[str], *, cwd: Path, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(shlex_quote(part) for part in command), flush=True)
    try:
        subprocess.run(command, cwd=cwd, env=env, check=True)
    except subprocess.CalledProcessError as exc:
        fail(f"command exited with status {exc.returncode}: {command[0]}")


def shlex_quote(value: str) -> str:
    import shlex

    return shlex.quote(value)


def ensure_descendant(path: Path, parent: Path) -> None:
    try:
        path.resolve().relative_to(parent.resolve())
    except ValueError:
        fail(f"refusing to clean path outside {parent}: {path}")


def remove_generated_path(path: Path, build_root: Path) -> None:
    ensure_descendant(path, build_root)
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def verify_environment(repo_root: Path) -> str:
    if not sys.platform.startswith("linux"):
        fail(f"this workflow only supports Linux; detected {sys.platform}")
    if platform.machine() not in {"x86_64", "amd64"}:
        fail(f"this release name expects x86_64; detected {platform.machine()}")

    missing_paths = [
        relative
        for relative in REQUIRED_PROJECT_PATHS
        if not (repo_root / relative).exists()
    ]
    if missing_paths:
        fail("required project paths are missing: " + ", ".join(missing_paths))

    if Path(sys.prefix).name.casefold() != "zaf":
        fail(
            "run this script with the zaf Conda environment's Python "
            f"(current interpreter: {sys.executable}, prefix: {sys.prefix})"
        )

    sys.path.insert(0, str(repo_root))
    import_errors: list[str] = []
    for module_name in REQUIRED_IMPORTS:
        try:
            importlib.import_module(module_name)
        except Exception as exc:
            import_errors.append(f"{module_name}: {exc}")
    if import_errors:
        fail("required imports failed:\n  " + "\n  ".join(import_errors))

    try:
        import tkinter as tk

        tk.Tcl().eval("info patchlevel")
    except Exception as exc:
        fail(f"the Tcl/Tk runtime is unavailable: {exc}")

    try:
        return importlib.metadata.version("pyinstaller")
    except importlib.metadata.PackageNotFoundError:
        fail("PyInstaller is not installed in the zaf environment")


def display_is_available(repo_root: Path) -> bool:
    if not os.environ.get("DISPLAY"):
        return False
    probe = [
        sys.executable,
        "-c",
        "import tkinter as tk; root=tk.Tk(); root.update(); root.destroy()",
    ]
    result = subprocess.run(
        probe,
        cwd=repo_root,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        timeout=10,
        check=False,
    )
    return result.returncode == 0


def write_linux_readme(
    destination: Path,
    *,
    pyinstaller_version: str,
    glibc_version: str,
) -> None:
    text = f"""ZAF — Linux x86_64 release
================================

Run ZAF
-------

1. Extract the complete archive.
2. Open a terminal in this directory.
3. Run: ./run_ZAF.sh

Paths containing spaces are supported. Keep ZAF, _internal, assets,
ZAF_instrument_settings.txt, and the launcher scripts together in this
directory.

Instrument defaults
-------------------

Edit ZAF_instrument_settings.txt in this directory to set the startup alpha
minimum/maximum, beta minimum/maximum, and image-to-holder rotation for your
TEM. Keep the five documented setting names unchanged and restart ZAF after
editing. If the file is missing or invalid, ZAF warns at startup and uses its
built-in defaults. Holder order remains xy and is not read from this file.

Application-menu launcher (optional)
------------------------------------

Run: ./install_launcher.sh

This installs a user-local desktop entry under
~/.local/share/applications and requires no sudo. The installer prints the
exact removal command when it finishes. Run it again if this extracted folder
is moved, because desktop entries use an absolute launcher path.

Compatibility
-------------

This bundle was built for Linux x86_64 on {platform.platform()} with
PyInstaller {pyinstaller_version}. The build host reports {glibc_version}.
It was verified on that host only. Linux binaries are not universally portable:
another machine may need a compatible Linux ABI and the same or a newer glibc,
plus a working graphical desktop/X11 environment. Wayland-only systems may
need XWayland for Tk.

Troubleshooting
---------------

Launch ./run_ZAF.sh from a terminal to retain useful error output. If the GUI
does not open, check that the archive was extracted on a Linux filesystem with
executable permissions preserved and that a graphical display is available.
"""
    destination.write_text(text, encoding="utf-8")


def collect_symlinks(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): os.readlink(path)
        for path in root.rglob("*")
        if path.is_symlink()
    }


def create_archive(release_root: Path, release_dir: Path, archive_path: Path) -> None:
    with tarfile.open(archive_path, "w:gz", format=tarfile.PAX_FORMAT) as archive:
        archive.add(release_dir, arcname=release_dir.name, recursive=True)

    original_symlinks = collect_symlinks(release_dir)
    with tempfile.TemporaryDirectory(prefix="zaf archive verify ") as temp_name:
        temp_root = Path(temp_name)
        with tarfile.open(archive_path, "r:gz") as archive:
            members = archive.getmembers()
            top_levels = {Path(member.name).parts[0] for member in members if member.name}
            if top_levels != {release_dir.name}:
                fail(f"archive has unexpected top-level entries: {sorted(top_levels)}")
            archive.extractall(temp_root, filter="data")
        extracted = temp_root / release_dir.name
        extracted_symlinks = collect_symlinks(extracted)
        if extracted_symlinks != original_symlinks:
            fail("archive did not preserve the release's symbolic links")

        verify_script = release_root.parent / "verify_linux_bundle.sh"
        verification_env = os.environ.copy()
        verification_env["ZAF_REQUIRE_GUI_SMOKE"] = "0"
        run(["bash", str(verify_script), str(extracted)], cwd=release_root.parent, env=verification_env)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-tests",
        action="store_true",
        help="skip the pre-build unittest suite (intended only for iterative local builds)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    packaging_dir = Path(__file__).resolve().parent
    repo_root = packaging_dir.parent
    build_root = packaging_dir / "build"
    work_dir = build_root / "work"
    spec_dir = build_root / "spec"
    pyinstaller_dist = build_root / "dist"
    release_root = packaging_dir / "release"
    release_dir = release_root / ARCHIVE_STEM
    archive_path = release_root / f"{ARCHIVE_STEM}.tar.gz"

    pyinstaller_version = verify_environment(repo_root)
    print(f"Using Python: {sys.executable}")
    print(f"Python version: {platform.python_version()}")
    print(f"Conda prefix: {sys.prefix}")
    print(f"PyInstaller version: {pyinstaller_version}")

    if not args.skip_tests:
        tests_dir = repo_root / "tests"
        if not tests_dir.is_dir():
            fail(f"test directory is missing: {tests_dir}")
        run(
            [sys.executable, "-m", "unittest", "discover", "-s", str(tests_dir), "-v"],
            cwd=repo_root,
        )

    for generated in (work_dir, spec_dir, pyinstaller_dist, release_dir, archive_path):
        remove_generated_path(generated, packaging_dir)
    for directory in (work_dir, spec_dir, pyinstaller_dist, release_root):
        directory.mkdir(parents=True, exist_ok=True)

    asset_data = f"{repo_root / 'assets'}{os.pathsep}assets"
    settings_data = (
        f"{repo_root / 'ZAF_instrument_settings.txt'}{os.pathsep}."
    )
    pyinstaller_command = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--onedir",
        "--windowed",
        "--clean",
        "--noconfirm",
        "--name",
        APP_NAME,
        "--workpath",
        str(work_dir),
        "--specpath",
        str(spec_dir),
        "--distpath",
        str(pyinstaller_dist),
        "--add-data",
        asset_data,
        "--add-data",
        settings_data,
        "--collect-data",
        "matplotlib",
        "--hidden-import",
        "matplotlib.backends.backend_tkagg",
        "--hidden-import",
        "PIL._tkinter_finder",
        "--exclude-module",
        "pytest",
        "--exclude-module",
        "tests",
        str(repo_root / "ZAF_gui.py"),
    ]
    run(pyinstaller_command, cwd=repo_root)

    built_app = pyinstaller_dist / APP_NAME
    if not (built_app / APP_NAME).is_file():
        fail(f"PyInstaller output is incomplete: {built_app}")
    shutil.copytree(built_app, release_dir, symlinks=True)

    top_level_assets = release_dir / "assets"
    top_level_assets.mkdir()
    shutil.copy2(repo_root / "assets" / "ZAF.png", top_level_assets / "ZAF.png")
    shutil.copy2(
        repo_root / "ZAF_instrument_settings.txt",
        release_dir / "ZAF_instrument_settings.txt",
    )
    for filename in ("run_ZAF.sh", "install_launcher.sh", "ZAF.desktop"):
        shutil.copy2(packaging_dir / filename, release_dir / filename)
    write_linux_readme(
        release_dir / "README_LINUX.txt",
        pyinstaller_version=pyinstaller_version,
        glibc_version=" ".join(platform.libc_ver()).strip() or "an unknown glibc version",
    )
    for filename in (APP_NAME, "run_ZAF.sh", "install_launcher.sh", "ZAF.desktop"):
        (release_dir / filename).chmod(0o755)

    gui_available = display_is_available(repo_root)
    verify_env = os.environ.copy()
    verify_env["ZAF_REQUIRE_GUI_SMOKE"] = "1" if gui_available else "0"
    verify_script = packaging_dir / "verify_linux_bundle.sh"
    run(["bash", str(verify_script), str(release_dir)], cwd=repo_root, env=verify_env)

    create_archive(release_root, release_dir, archive_path)
    print(f"GUI display available for smoke test: {'yes' if gui_available else 'no'}")
    print(f"Release directory: {release_dir}")
    print(f"Archive: {archive_path} ({archive_path.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
