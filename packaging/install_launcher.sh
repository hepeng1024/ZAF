#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
RUNNER="$SCRIPT_DIR/run_ZAF.sh"
ICON="$SCRIPT_DIR/assets/ZAF.png"
APPLICATIONS_DIR="${XDG_DATA_HOME:-$HOME/.local/share}/applications"
DESKTOP_FILE="$APPLICATIONS_DIR/ZAF.desktop"

if [[ ! -x "$RUNNER" ]]; then
    printf 'Launcher installation failed: %s is missing or not executable.\n' "$RUNNER" >&2
    exit 1
fi
if [[ ! -r "$ICON" ]]; then
    printf 'Launcher installation failed: icon is missing or unreadable: %s\n' "$ICON" >&2
    exit 1
fi
if [[ "$RUNNER" == *$'\n'* || "$RUNNER" == *$'\r'* || "$ICON" == *$'\n'* || "$ICON" == *$'\r'* ]]; then
    printf 'Launcher installation failed: the installation path contains a newline.\n' >&2
    exit 1
fi

desktop_exec_quote() {
    local value="$1"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    value="${value//\`/\\\`}"
    value="${value//\$/\\\$}"
    value="${value//%/%%}"
    printf '"%s"' "$value"
}

mkdir -p -- "$APPLICATIONS_DIR"
QUOTED_RUNNER="$(desktop_exec_quote "$RUNNER")"
{
    printf '%s\n' '[Desktop Entry]'
    printf '%s\n' 'Version=1.0'
    printf '%s\n' 'Type=Application'
    printf '%s\n' 'Name=ZAF'
    printf '%s\n' 'Comment=TEM Zone-Axis Finder'
    printf 'Exec=%s\n' "$QUOTED_RUNNER"
    printf 'Icon=%s\n' "$ICON"
    printf '%s\n' 'Terminal=false'
    printf '%s\n' 'Categories=Science;Education;'
    printf '%s\n' 'StartupNotify=true'
} > "$DESKTOP_FILE"
chmod 755 -- "$DESKTOP_FILE"

if command -v desktop-file-validate >/dev/null 2>&1; then
    desktop-file-validate "$DESKTOP_FILE"
fi
if command -v update-desktop-database >/dev/null 2>&1; then
    update-desktop-database "$APPLICATIONS_DIR"
fi

printf 'Installed the ZAF application-menu launcher at:\n  %s\n' "$DESKTOP_FILE"
printf 'It launches:\n  %s\n' "$RUNNER"
printf 'To remove it later, run:\n  rm -- %q\n' "$DESKTOP_FILE"
printf 'Then refresh the menu, if available, with:\n  update-desktop-database %q\n' "$APPLICATIONS_DIR"
