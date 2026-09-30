#!/usr/bin/env bash
# Runs all CTest tests and example binaries instrumented for LLVM coverage,
# merges the .profraw files, and generates an HTML + text report.
#
# Requires a build configured with -DCOVERAGE=ON -DENABLE_TESTS=ON -WITH_EXAMPLES=ON.
#
# IMPORTANT — demo examples require a prior install:
#   cmake --install <BUILD_DIR>
#
# The demo binaries look for data files (HepMC event files, YAML cards) in
# the examples/data/ subdirectory of their working directory.  Those files
# are only placed there by `cmake --install`; they are NOT present in the
# build tree.  Demos are therefore run from the installed examples directory
# ($CMAKE_INSTALL_PREFIX/share/Hammer/examples) using the instrumented binary
# from the build tree.  Any demo whose install directory cannot be found is
# skipped with a warning.
#
# Usage: coverage.sh [BUILD_DIR] [SOURCE_DIR] [REPORT_DIR]
#   BUILD_DIR   defaults to ~/BuildTests/Hammer/Profile
#   SOURCE_DIR  defaults to the directory containing this script's parent
#   REPORT_DIR  defaults to BUILD_DIR/coverage_report

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="${3:-$(cd "$SCRIPT_DIR/.." && pwd)}"
BUILD_DIR="${1:-$HOME/BuildTests/Hammer/Profile}"
REPORT_DIR="${2:-$BUILD_DIR/coverage_report}"
PROFRAW_DIR="$REPORT_DIR/profraw"

BIN_DIR="$BUILD_DIR/bin"

# Derive the installed examples directory from CMakeCache.txt.
# cmake --install places demo data files there; demos must run from that dir.
CMAKE_CACHE="$BUILD_DIR/CMakeCache.txt"
INSTALL_PREFIX=""
INSTALL_DATADIR=""
if [[ -f "$CMAKE_CACHE" ]]; then
    INSTALL_PREFIX="$(grep -m1 '^CMAKE_INSTALL_PREFIX:' "$CMAKE_CACHE" | cut -d= -f2)"
    INSTALL_DATADIR="$(grep -m1 '^CMAKE_INSTALL_DATADIR:' "$CMAKE_CACHE" | cut -d= -f2)"
fi
# GNUInstallDirs defaults CMAKE_INSTALL_DATADIR to "share" when left empty
INSTALL_DATADIR="${INSTALL_DATADIR:-share}"
EXAMPLES_INSTALL_DIR="${INSTALL_PREFIX:+$INSTALL_PREFIX/$INSTALL_DATADIR/Hammer/examples}"

# ---------------------------------------------------------------------------
# Tool detection (prefer xcrun on macOS)
# ---------------------------------------------------------------------------
find_tool() {
    local name="$1"
    if command -v "$name" &>/dev/null; then
        echo "$name"
    elif xcrun --find "$name" &>/dev/null 2>&1; then
        xcrun --find "$name"
    else
        echo ""
    fi
}

LLVM_PROFDATA="$(find_tool llvm-profdata)"
LLVM_COV="$(find_tool llvm-cov)"

if [[ -z "$LLVM_PROFDATA" || -z "$LLVM_COV" ]]; then
    echo "ERROR: llvm-profdata or llvm-cov not found. Install Xcode command-line tools." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
if [[ ! -d "$BUILD_DIR" ]]; then
    echo "ERROR: build directory not found: $BUILD_DIR" >&2
    echo "Configure with -DCOVERAGE=ON -DENABLE_TESTS=ON first." >&2
    exit 1
fi

if [[ ! -d "$BIN_DIR" ]]; then
    echo "ERROR: no bin/ directory in $BUILD_DIR — build the project first." >&2
    exit 1
fi

mkdir -p "$PROFRAW_DIR"

echo "=== Hammer coverage run ==="
echo "  Build dir   : $BUILD_DIR"
echo "  Source dir  : $SOURCE_DIR"
echo "  Report dir  : $REPORT_DIR"
echo "  Examples dir: ${EXAMPLES_INSTALL_DIR:-<not found — demos will be skipped>}"
echo

# ---------------------------------------------------------------------------
# Discover tests via ctest --show-only=json-v1
# The JSON format has properties as a list of {name, value} objects.
# ---------------------------------------------------------------------------
CTEST_JSON="$REPORT_DIR/ctest_list.json"
if ! ctest --test-dir "$BUILD_DIR" --show-only=json-v1 2>/dev/null > "$CTEST_JSON"; then
    echo "WARNING: ctest --show-only=json-v1 failed; falling back to name-based discovery." >&2
    CTEST_JSON=""
fi

# ---------------------------------------------------------------------------
# Run a single binary and capture its .profraw
# ---------------------------------------------------------------------------
PROFRAW_FILES=()
ALL_OBJECTS=()
PASS=0
FAIL=0

run_binary() {
    local label="$1"
    local binary="$2"
    local workdir="$3"

    if [[ ! -x "$binary" ]]; then
        echo "  [SKIP] $label — not found or not executable: $binary"
        return
    fi

    local safe_label="${label//\//_}"
    local profraw="$PROFRAW_DIR/${safe_label}.profraw"
    local log="$PROFRAW_DIR/${safe_label}.log"

    printf "  %-55s " "$label"

    if ( cd "$workdir" && env LLVM_PROFILE_FILE="$profraw" "$binary" ) > "$log" 2>&1; then
        echo "PASS"
        PASS=$(( PASS + 1 ))
    else
        local rc=$?
        echo "FAIL (exit $rc)"
        FAIL=$(( FAIL + 1 ))
    fi

    if [[ -f "$profraw" ]]; then
        PROFRAW_FILES+=("$profraw")
        ALL_OBJECTS+=("$binary")
    else
        echo "    WARNING: no .profraw written for $label" >&2
    fi
}

# ---------------------------------------------------------------------------
# 1. CTest-registered tests
# ---------------------------------------------------------------------------
echo "--- Tests (via ctest) ---"

if [[ -n "$CTEST_JSON" && -f "$CTEST_JSON" ]]; then
    # properties in json-v1 is a list of {name, value} objects, not a dict
    while IFS='|' read -r tname texe twd; do
        run_binary "$tname" "$texe" "$twd"
    done < <(python3 - "$CTEST_JSON" "$BIN_DIR" <<'PYEOF'
import json, sys, os

data = json.load(open(sys.argv[1]))
bin_dir = sys.argv[2]

for t in data.get("tests", []):
    name = t["name"]
    cmd  = t.get("command", [])
    props = {p["name"]: p["value"] for p in t.get("properties", [])}
    wd = props.get("WORKING_DIRECTORY", bin_dir)
    if not wd:
        wd = bin_dir
    exe = cmd[0] if cmd else os.path.join(bin_dir, name)
    if not os.path.isabs(exe):
        exe = os.path.join(bin_dir, exe)
    print(name + "|" + exe + "|" + wd)
PYEOF
)
else
    # Fallback: find Test* executables in bin/
    while IFS= read -r -d '' binary; do
        name="$(basename "$binary")"
        run_binary "$name" "$binary" "$BIN_DIR"
    done < <(find "$BIN_DIR" -maxdepth 1 -type f -perm +111 -name 'Test*' -print0 | sort -z)
fi

# ---------------------------------------------------------------------------
# 2. Example (demo) binaries
#
# Demos are run from the installed examples directory so that they can find
# their data files (HepMC files, YAML cards).  The instrumented binary from
# the build tree is used so that coverage is captured.  If the install
# directory is missing, all demos are skipped.
# ---------------------------------------------------------------------------
echo
echo "--- Examples (demo*) ---"

if [[ -z "$EXAMPLES_INSTALL_DIR" || ! -d "$EXAMPLES_INSTALL_DIR" ]]; then
    echo "  [SKIP ALL] installed examples directory not found."
    echo "             Run 'cmake --install $BUILD_DIR' first, then re-run this script."
else
    while IFS= read -r -d '' binary; do
        name="$(basename "$binary")"
        run_binary "example/$name" "$binary" "$EXAMPLES_INSTALL_DIR"
    done < <(find "$BIN_DIR" -maxdepth 1 -type f -perm +111 -name 'demo*' -print0 2>/dev/null | sort -z)
fi

# ---------------------------------------------------------------------------
# Summary of runs
# ---------------------------------------------------------------------------
echo
echo "=== Run summary: $PASS passed, $FAIL failed ==="
echo "    profraw files collected: ${#PROFRAW_FILES[@]}"
echo

if [[ ${#PROFRAW_FILES[@]} -eq 0 ]]; then
    echo "ERROR: no .profraw files generated. Was the build configured with -DCOVERAGE=ON?" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Merge .profraw → .profdata
# ---------------------------------------------------------------------------
PROFDATA="$REPORT_DIR/merged.profdata"
echo "Merging ${#PROFRAW_FILES[@]} profraw files → $PROFDATA"
"$LLVM_PROFDATA" merge -sparse "${PROFRAW_FILES[@]}" -o "$PROFDATA"

# ---------------------------------------------------------------------------
# Select objects for llvm-cov.
#
# In a shared library build, the same header-defined functions (inline,
# templates, exception types) get compiled into every test/demo binary AND
# into the dylibs.  Passing all 170 binaries as objects causes llvm-cov to
# encounter the same function with different coverage-mapping hashes across
# binaries → "mismatched data" warning for each duplicate.
#
# Passing only the dylibs avoids this: each function exists once per dylib,
# and the profraw files collected from running every test already contain
# accumulated coverage for all dylib code (the dylibs are loaded into each
# test process at run time).
#
# For static builds there are no dylibs, so fall back to the test/demo
# binaries collected during the run (original behaviour).
# ---------------------------------------------------------------------------
LIB_DIR="$BUILD_DIR/lib"
UNIQUE_OBJECTS=()

if [[ -d "$LIB_DIR" ]]; then
    while IFS= read -r -d '' dylib; do
        UNIQUE_OBJECTS+=("$dylib")
    done < <(find "$LIB_DIR" -maxdepth 1 -type f -name 'lib*.dylib' -print0 | sort -z)
fi

if [[ ${#UNIQUE_OBJECTS[@]} -gt 0 ]]; then
    echo "Shared build detected — using ${#UNIQUE_OBJECTS[@]} dylibs as coverage objects"
    for obj in "${UNIQUE_OBJECTS[@]}"; do echo "  $(basename "$obj")"; done
else
    echo "Static build detected — using test/demo binaries as coverage objects"
    declare -A SEEN_OBJECTS
    for obj in "${ALL_OBJECTS[@]}"; do
        if [[ -z "${SEEN_OBJECTS[$obj]+set}" ]]; then
            SEEN_OBJECTS[$obj]=1
            UNIQUE_OBJECTS+=("$obj")
        fi
    done
    echo "  ${#UNIQUE_OBJECTS[@]} unique binaries"
fi
echo

# Source filter: restrict report to Hammer source, exclude third-party
IGNORE_REGEX='(Tests/GTest|external|yaml-cpp|flatbuffers|gtest)'

run_llvm_cov() {
    local subcmd="$1"; shift
    local first_obj="${UNIQUE_OBJECTS[0]}"
    local extra_objects=()
    for obj in "${UNIQUE_OBJECTS[@]:1}"; do
        extra_objects+=(-object "$obj")
    done
    "$LLVM_COV" "$subcmd" \
        "$first_obj" \
        "${extra_objects[@]+"${extra_objects[@]}"}" \
        -instr-profile="$PROFDATA" \
        -ignore-filename-regex="$IGNORE_REGEX" \
        "$@"
}

# ---------------------------------------------------------------------------
# Text report (stdout + file)
# ---------------------------------------------------------------------------
TEXT_REPORT="$REPORT_DIR/coverage.txt"
echo "Generating text report  → $TEXT_REPORT"
run_llvm_cov report --use-color=false | tee "$TEXT_REPORT"

# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
HTML_DIR="$REPORT_DIR/html"
echo
echo "Generating HTML report  → $HTML_DIR/index.html"
run_llvm_cov show \
    -format=html \
    -output-dir="$HTML_DIR" \
    -show-line-counts-or-regions \
    -show-expansions

# ---------------------------------------------------------------------------
# Export JSON (machine-readable)
# ---------------------------------------------------------------------------
JSON_REPORT="$REPORT_DIR/coverage.json"
echo "Generating JSON export  → $JSON_REPORT"
run_llvm_cov export -format=text > "$JSON_REPORT"

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo
echo "=== Coverage report complete ==="
echo "  Text  : $TEXT_REPORT"
echo "  HTML  : $HTML_DIR/index.html"
echo "  JSON  : $JSON_REPORT"
echo "  Data  : $PROFDATA"
