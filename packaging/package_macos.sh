#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd -P)"
APP_NAME="ZAF"
ENTRY_SCRIPT="$REPO_ROOT/ZAF_gui.py"
BACKEND_MODULE="$REPO_ROOT/ZAF.py"
SETTINGS_FILE="$REPO_ROOT/ZAF_instrument_settings.txt"
MACOS_README="$SCRIPT_DIR/README_MACOS.txt"
ASSET_DIR="$REPO_ROOT/assets"
ICON_FILE="$ASSET_DIR/ZAF.icns"
BUNDLE_IDENTIFIER="edu.umich.hepeng.ZAF"
BUILD_ROOT="$SCRIPT_DIR/build/macos"
WORK_DIR="$BUILD_ROOT/work"
SPEC_DIR="$BUILD_ROOT/spec"
DIST_DIR="$BUILD_ROOT/dist"
RELEASE_ROOT="$SCRIPT_DIR/release/macos"
RELEASE_DIR="$RELEASE_ROOT/${APP_NAME}-macOS-arm64"
APP_BUNDLE="$DIST_DIR/${APP_NAME}.app"
ARCHIVE_PATH="$RELEASE_ROOT/${APP_NAME}-macOS-arm64.zip"

fail() {
    printf 'macOS packaging failed: %s\n' "$*" >&2
    exit 1
}

require_file() {
    [[ -f "$1" ]] || fail "required file is missing: $1"
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || fail "required command is unavailable: $1"
}

check_repository_inputs() {
    require_file "$ENTRY_SCRIPT"
    require_file "$BACKEND_MODULE"
    require_file "$SETTINGS_FILE"
    require_file "$MACOS_README"
    require_file "$ICON_FILE"
    require_file "$ASSET_DIR/ZAF.png"
    require_file "$ASSET_DIR/FCC.png"
    require_file "$ASSET_DIR/BCC.png"
    require_file "$ASSET_DIR/HCP.png"
    [[ -d "$ASSET_DIR" ]] || fail "asset directory is missing: $ASSET_DIR"
}

clean_generated_directory() {
    local target="$1"
    if [[ "$target" != "$BUILD_ROOT" && "$target" != "$RELEASE_ROOT" ]]; then
        fail "refusing to clean unexpected path: $target"
    fi
    rm -rf -- "$target"
}

check_repository_inputs

if [[ "${1:-}" == "--check-only" ]]; then
    [[ $# -eq 1 ]] || fail "--check-only does not accept additional arguments"
    printf 'macOS packaging inputs: OK\n'
    printf '  repository: %s\n' "$REPO_ROOT"
    printf '  entry point: %s\n' "$ENTRY_SCRIPT"
    printf '  instrument settings: %s\n' "$SETTINGS_FILE"
    printf '  macOS readme: %s\n' "$MACOS_README"
    printf '  icon: %s\n' "$ICON_FILE"
    exit 0
fi
[[ $# -eq 0 ]] || fail "unknown argument: $1"

[[ "$(uname -s)" == "Darwin" ]] \
    || fail "this script must run natively on macOS; detected $(uname -s)"
[[ "$(uname -m)" == "arm64" ]] \
    || fail "an Apple-silicon ARM64 runner is required; detected $(uname -m)"

for command_name in python ditto file lipo plutil codesign shasum cmp; do
    require_command "$command_name"
done

python_arch="$(python -c 'import platform; print(platform.machine())')"
[[ "$python_arch" == "arm64" ]] \
    || fail "Python must run natively as ARM64; detected $python_arch"

(
    cd "$REPO_ROOT"
    python -c 'import tkinter as tk; import numpy, scipy, PIL, matplotlib, ZAF, ZAF_gui; tk.Tcl().eval("info patchlevel")'
    python -m PyInstaller --version
)

clean_generated_directory "$BUILD_ROOT"
clean_generated_directory "$RELEASE_ROOT"
mkdir -p -- "$WORK_DIR" "$SPEC_DIR" "$DIST_DIR" "$RELEASE_ROOT"

(
    cd "$REPO_ROOT"
    python -m PyInstaller \
        --onedir \
        --windowed \
        --clean \
        --noconfirm \
        --name "$APP_NAME" \
        --target-arch arm64 \
        --icon "$ICON_FILE" \
        --osx-bundle-identifier "$BUNDLE_IDENTIFIER" \
        --workpath "$WORK_DIR" \
        --specpath "$SPEC_DIR" \
        --distpath "$DIST_DIR" \
        --add-data "$ASSET_DIR:assets" \
        --add-data "$SETTINGS_FILE:." \
        --collect-data matplotlib \
        --hidden-import matplotlib.backends.backend_tkagg \
        --hidden-import PIL._tkinter_finder \
        --exclude-module pytest \
        --exclude-module tests \
        "$ENTRY_SCRIPT"
)

MAIN_EXECUTABLE="$APP_BUNDLE/Contents/MacOS/$APP_NAME"
INFO_PLIST="$APP_BUNDLE/Contents/Info.plist"
[[ -d "$APP_BUNDLE" ]] || fail "PyInstaller did not create $APP_BUNDLE"
[[ -x "$MAIN_EXECUTABLE" ]] || fail "app executable is missing or not executable"
require_file "$INFO_PLIST"

plutil -lint "$INFO_PLIST"
actual_identifier="$(/usr/libexec/PlistBuddy -c 'Print :CFBundleIdentifier' "$INFO_PLIST")"
[[ "$actual_identifier" == "$BUNDLE_IDENTIFIER" ]] \
    || fail "unexpected bundle identifier: $actual_identifier"

main_architectures="$(lipo -archs "$MAIN_EXECUTABLE")"
[[ "$main_architectures" == "arm64" ]] \
    || fail "main executable is not ARM64-only: $main_architectures"

mach_o_count=0
while IFS= read -r -d '' candidate; do
    description="$(file -b "$candidate")"
    [[ "$description" == *"Mach-O"* ]] || continue
    architectures="$(lipo -archs "$candidate")" \
        || fail "could not inspect architecture: $candidate"
    case " $architectures " in
        *" arm64 "*) ;;
        *) fail "Mach-O file lacks ARM64 support: $candidate ($architectures)" ;;
    esac
    if [[ " $architectures " == *" x86_64 "* ]]; then
        fail "ARM64 release contains an x86_64 slice: $candidate ($architectures)"
    fi
    ((mach_o_count += 1))
done < <(find "$APP_BUNDLE" -type f -print0)
((mach_o_count > 0)) || fail "no Mach-O files were found in the app bundle"
printf 'Verified %d ARM64-only Mach-O files.\n' "$mach_o_count"

for source_asset in "$ASSET_DIR"/*.png; do
    asset_name="$(basename "$source_asset")"
    bundled_matches="$(find "$APP_BUNDLE" -path "*/assets/$asset_name" -print)"
    [[ -n "$bundled_matches" ]] || fail "bundled PNG asset is missing: $asset_name"
done

bundled_settings="$(find "$APP_BUNDLE" -name 'ZAF_instrument_settings.txt' -print -quit)"
[[ -n "$bundled_settings" && -r "$bundled_settings" ]] \
    || fail "bundled instrument settings template is missing"
cmp -s "$SETTINGS_FILE" "$bundled_settings" \
    || fail "bundled instrument settings template differs from the repository file"

if [[ -n "$(find "$APP_BUNDLE" -type d \( -name tests -o -name .git -o -name __pycache__ -o -name .pytest_cache \) -print)" ]]; then
    fail "the app bundle contains repository tests, VCS metadata, or development caches"
fi

"$MAIN_EXECUTABLE" --bundle-self-test \
    || fail "the frozen application's headless bundle self-test failed"
codesign --verify --deep --strict --verbose=2 "$APP_BUNDLE"

mkdir -p -- "$RELEASE_DIR"
ditto "$APP_BUNDLE" "$RELEASE_DIR/${APP_NAME}.app"
cp "$MACOS_README" "$RELEASE_DIR/README_MACOS.txt"
ditto -c -k --sequesterRsrc --keepParent "$RELEASE_DIR" "$ARCHIVE_PATH"
[[ -s "$ARCHIVE_PATH" ]] || fail "release ZIP was not created"

VERIFY_DIR="$(mktemp -d -t zaf-macos-archive-verify)"
cleanup_verify_dir() {
    rm -rf -- "$VERIFY_DIR"
}
trap cleanup_verify_dir EXIT
ditto -x -k "$ARCHIVE_PATH" "$VERIFY_DIR"
EXTRACTED_RELEASE="$VERIFY_DIR/${APP_NAME}-macOS-arm64"
EXTRACTED_APP="$EXTRACTED_RELEASE/${APP_NAME}.app"
[[ -x "$EXTRACTED_APP/Contents/MacOS/$APP_NAME" ]] \
    || fail "archive did not preserve the app executable"
[[ -r "$EXTRACTED_RELEASE/README_MACOS.txt" ]] \
    || fail "archive did not preserve README_MACOS.txt"
cmp -s "$MACOS_README" "$EXTRACTED_RELEASE/README_MACOS.txt" \
    || fail "archive changed README_MACOS.txt"
[[ ! -e "$EXTRACTED_RELEASE/ZAF_instrument_settings.txt" ]] \
    || fail "archive contains a misleading settings file beside ZAF.app"
[[ "$(lipo -archs "$EXTRACTED_APP/Contents/MacOS/$APP_NAME")" == "arm64" ]] \
    || fail "archive did not preserve the ARM64 executable"
codesign --verify --deep --strict --verbose=2 "$EXTRACTED_APP"

if ! diff -u \
    <(cd "$APP_BUNDLE" && find . -type l -exec sh -c 'for item do printf "%s -> %s\n" "$item" "$(readlink "$item")"; done' sh {} + | sort) \
    <(cd "$EXTRACTED_APP" && find . -type l -exec sh -c 'for item do printf "%s -> %s\n" "$item" "$(readlink "$item")"; done' sh {} + | sort); then
    fail "archive did not preserve the app bundle's symbolic links"
fi

printf 'macOS ARM64 package: OK\n'
printf '  app: %s\n' "$APP_BUNDLE"
printf '  archive: %s\n' "$ARCHIVE_PATH"
printf '  SHA-256: '
shasum -a 256 "$ARCHIVE_PATH" | awk '{print $1}'
