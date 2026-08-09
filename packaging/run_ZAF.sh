#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
EXECUTABLE="$SCRIPT_DIR/ZAF"

if [[ ! -x "$EXECUTABLE" ]]; then
    printf 'ZAF launcher error: executable not found or not executable: %s\n' "$EXECUTABLE" >&2
    exit 1
fi

exec "$EXECUTABLE" "$@"
