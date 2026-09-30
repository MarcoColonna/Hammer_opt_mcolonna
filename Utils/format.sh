#!/usr/bin/env bash
# Runs clang-format -i on all Hammer source files:
#   - include/Hammer/**/*.hh and *.hh.in
#   - src/**/*.cc
#   - Tests/**/Test*.cc  (GTest subdirectory excluded)
#
# The .clang-format config at the repository root is picked up automatically.
#
# Usage: format.sh [SOURCE_DIR]
#   SOURCE_DIR  defaults to the directory containing this script's parent

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="${1:-$(cd "$SCRIPT_DIR/.." && pwd)}"

# ---------------------------------------------------------------------------
# Tool detection (prefer xcrun on macOS)
# ---------------------------------------------------------------------------
CLANG_FORMAT=""
if command -v clang-format &>/dev/null; then
    CLANG_FORMAT="clang-format"
elif xcrun --find clang-format &>/dev/null 2>&1; then
    CLANG_FORMAT="$(xcrun --find clang-format)"
else
    echo "ERROR: clang-format not found. Install Xcode command-line tools or LLVM." >&2
    exit 1
fi

echo "=== Hammer clang-format ==="
echo "  Source dir    : $SOURCE_DIR"
echo "  clang-format  : $CLANG_FORMAT ($($CLANG_FORMAT --version))"
echo

# ---------------------------------------------------------------------------
# Collect files
# ---------------------------------------------------------------------------
mapfile -d '' HEADERS < <(
    find "$SOURCE_DIR/include/Hammer" \( -name "*.hh" -o -name "*.hh.in" \) -print0 | sort -z
)

mapfile -d '' SOURCES < <(
    find "$SOURCE_DIR/src" -name "*.cc" -print0 | sort -z
)

mapfile -d '' TESTS < <(
    find "$SOURCE_DIR/Tests" -path "*/GTest" -prune -o -name "Test*.cc" -print0 | sort -z
)

TOTAL=$(( ${#HEADERS[@]} + ${#SOURCES[@]} + ${#TESTS[@]} ))
echo "  Files found   : $TOTAL  (${#HEADERS[@]} headers, ${#SOURCES[@]} sources, ${#TESTS[@]} tests)"
echo

# ---------------------------------------------------------------------------
# Format
# ---------------------------------------------------------------------------
format_group() {
    local label="$1"; shift
    local files=("$@")
    if [[ ${#files[@]} -eq 0 ]]; then return; fi
    echo "--- $label (${#files[@]} files) ---"
    for f in "${files[@]}"; do
        printf "  %s\n" "${f#$SOURCE_DIR/}"
        "$CLANG_FORMAT" -i "$f"
    done
    echo
}

format_group "Headers  (include/Hammer)" "${HEADERS[@]}"
format_group "Sources  (src)"            "${SOURCES[@]}"
format_group "Tests    (Tests)"          "${TESTS[@]}"

echo "=== Done: $TOTAL files formatted ==="
