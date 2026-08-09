#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
DEFAULT_RELEASE="$SCRIPT_DIR/release/ZAF-Linux-x86_64"
RELEASE_DIR="${1:-$DEFAULT_RELEASE}"

fail() {
    printf 'VERIFY FAILED: %s\n' "$*" >&2
    exit 1
}

pass() {
    printf 'VERIFY OK: %s\n' "$*"
}

[[ -d "$RELEASE_DIR" ]] || fail "release directory not found: $RELEASE_DIR"
RELEASE_DIR="$(cd -- "$RELEASE_DIR" && pwd -P)"

for executable in ZAF run_ZAF.sh install_launcher.sh; do
    [[ -f "$RELEASE_DIR/$executable" ]] || fail "missing $executable"
    [[ -x "$RELEASE_DIR/$executable" ]] || fail "$executable is not executable"
done
pass "main executable and launcher scripts are executable"

[[ -r "$RELEASE_DIR/assets/ZAF.png" ]] || fail "top-level Linux icon is missing or unreadable"
pass "desktop icon is present"

SETTINGS_FILE="ZAF_instrument_settings.txt"
[[ -r "$RELEASE_DIR/$SETTINGS_FILE" ]] \
    || fail "editable instrument settings file is missing or unreadable"
[[ -r "$RELEASE_DIR/_internal/$SETTINGS_FILE" ]] \
    || fail "PyInstaller instrument settings template is missing or unreadable"
pass "editable and bundled instrument settings files are present"

if [[ -f "$RELEASE_DIR/ZAF.desktop" ]]; then
    grep -Eq '^Exec=.*run_ZAF\.sh' "$RELEASE_DIR/ZAF.desktop" \
        || fail "ZAF.desktop does not call run_ZAF.sh"
    if grep -Eq '^Exec=[^#]*([[:space:]/]|^)ZAF([[:space:]]|$)' "$RELEASE_DIR/ZAF.desktop" \
        && ! grep -q 'run_ZAF.sh' "$RELEASE_DIR/ZAF.desktop"; then
        fail "ZAF.desktop bypasses run_ZAF.sh"
    fi
    pass "ZAF.desktop routes startup through run_ZAF.sh"
fi

RUNTIME_ASSET_DIR="$RELEASE_DIR/_internal/assets"
[[ -d "$RUNTIME_ASSET_DIR" ]] || fail "runtime asset directory missing: $RUNTIME_ASSET_DIR"
runtime_assets=(
    "BCC.png"
    "FCC.png"
    "HCP.png"
    "Negative alpha tilt arrow.png"
    "Negative beta tilt arrow.png"
    "Positive alpha tilt arrow.png"
    "Positive beta tilt arrow.png"
    "Sample holder object.png"
    "TEM lamella object.png"
    "Schematic of the double-tilt holder and the zone-axis of the sample.png"
    "Schematic of the effect of sample rotation.png"
)
for asset in "${runtime_assets[@]}"; do
    [[ -r "$RUNTIME_ASSET_DIR/$asset" ]] || fail "missing runtime asset: $asset"
done
pass "all ${#runtime_assets[@]} GUI runtime PNG assets are present"

if find "$RELEASE_DIR" -type d \
    \( -name tests -o -name .git -o -name __pycache__ -o -name .pytest_cache \) \
    -print -quit | grep -q .; then
    fail "a development/test/cache directory is present"
fi
if find "$RELEASE_DIR" -type f \( -name '*.spec' -o -name '*.pyc' \) \
    -print -quit | grep -q .; then
    fail "a generated spec or Python cache file is present"
fi
pass "tests, VCS metadata, and development caches are absent"

self_test_output="$("$RELEASE_DIR/run_ZAF.sh" --bundle-self-test 2>&1)" \
    || fail "bundled import/data self-test failed: $self_test_output"
printf '%s\n' "$self_test_output"
grep -q '^ZAF bundle self-test: OK$' <<< "$self_test_output" \
    || fail "bundle self-test did not report success"
for component in 'ZAF backend:' 'NumPy:' 'SciPy:' 'Pillow:' 'Matplotlib:' 'Tcl/Tk:'; do
    grep -q "$component" <<< "$self_test_output" \
        || fail "bundle self-test did not confirm $component"
done
pass "scientific, imaging, plotting, Tkinter, and local backend imports work"

launcher_test_root="$(mktemp -d -t zaf-launcher-verify-XXXXXX)"
cleanup_launcher_test() {
    rm -rf -- "$launcher_test_root"
}
trap cleanup_launcher_test EXIT
installer_output="$(XDG_DATA_HOME="$launcher_test_root/share" "$RELEASE_DIR/install_launcher.sh" 2>&1)" \
    || fail "user-local launcher installer failed: $installer_output"
installed_desktop="$launcher_test_root/share/applications/ZAF.desktop"
[[ -x "$installed_desktop" ]] || fail "launcher installer did not create an executable desktop entry"
grep -Fq "run_ZAF.sh" "$installed_desktop" \
    || fail "installed desktop entry does not call run_ZAF.sh"
grep -Fq "Icon=$RELEASE_DIR/assets/ZAF.png" "$installed_desktop" \
    || fail "installed desktop entry does not use the release's PNG icon"
pass "no-sudo launcher installer creates an absolute user-local desktop entry"

ldd_checked=0
runtime_library_path="$RELEASE_DIR/_internal"
while IFS= read -r -d '' candidate; do
    file_description="$(LC_ALL=C file -b -- "$candidate")"
    if [[ "$file_description" != ELF* ]]; then
        continue
    fi
    if [[ "$file_description" != *"dynamically linked"* && "$file_description" != *"shared object"* ]]; then
        continue
    fi
    # PyInstaller's bootloader prepends _internal to LD_LIBRARY_PATH before it
    # loads Python. Use that same runtime context when inspecting extensions;
    # several Conda-built modules retain RPATHs relative to their old prefix.
    ldd_output="$(LD_LIBRARY_PATH="$runtime_library_path" ldd "$candidate" 2>&1)" \
        || fail "ldd could not inspect ${candidate#"$RELEASE_DIR/"}: $ldd_output"
    if grep -q 'not found' <<< "$ldd_output"; then
        printf '%s\n' "$ldd_output" >&2
        fail "unresolved library dependency in ${candidate#"$RELEASE_DIR/"}"
    fi
    ((ldd_checked += 1))
done < <(find "$RELEASE_DIR" -type f -print0)
((ldd_checked > 0)) || fail "no dynamically linked ELF files were found"
pass "ldd found no unresolved runtime libraries across $ldd_checked ELF files"

if [[ "${ZAF_REQUIRE_GUI_SMOKE:-0}" == "1" ]]; then
    command -v timeout >/dev/null 2>&1 || fail "timeout is needed for the GUI smoke test"
    gui_output="$(timeout 20 "$RELEASE_DIR/run_ZAF.sh" --gui-smoke-test 2>&1)" \
        || fail "GUI smoke test failed: $gui_output"
    printf '%s\n' "$gui_output"
    grep -q '^ZAF GUI smoke test: OK$' <<< "$gui_output" \
        || fail "GUI smoke test did not report success"
    pass "bundled GUI opened and closed cleanly"
else
    printf 'VERIFY SKIP: GUI smoke test (set ZAF_REQUIRE_GUI_SMOKE=1 when a display is available)\n'
fi

pass "bundle verification completed for $RELEASE_DIR"
