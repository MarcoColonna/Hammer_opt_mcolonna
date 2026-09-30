#!/usr/bin/env python3
"""
bump_version.py — prepare a new HAMMER release

Usage (run from any directory):
    python Utils/bump_version.py <new_version> <new_zenodo_id>

Arguments
---------
new_version   : dot-separated X.Y.Z (exactly three numeric components)
new_zenodo_id : Zenodo record identifier in any of these forms:
                  12345678                 bare record ID
                  zenodo.12345678          prefixed
                  10.5281/zenodo.12345678  full DOI

Actions
-------
1. Bumps HAMMER_VERSION_MAJOR / _MINOR / _PATCH in CMakeLists.txt.
2. Updates "ver. X.Y.Z" in README.md and
   "HAMMER RELEASE VERSION X.Y.Z" in the plain-text README.
3. Replaces the version-specific Zenodo DOI badge (second badge) in README.md.
4. Prepends a dated stub entry to CHANGELOG.
5. Updates the final year in all "Copyright (C) YYYY - YYYY The HAMMER
   Collaboration" headers across the source tree to the current calendar year.

The script must be run from within the Hammer source tree (or pass --root to
specify the repository root explicitly).
"""

import argparse
import re
from datetime import date
from pathlib import Path


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def find_repo_root(start: Path) -> Path:
    """Walk upward from *start* until we find CMakeLists.txt + CHANGELOG."""
    for parent in [start, *start.parents]:
        if (parent / "CMakeLists.txt").exists() and (parent / "CHANGELOG").exists():
            return parent
    raise SystemExit(
        "error: could not locate the Hammer repository root "
        "(expected CMakeLists.txt and CHANGELOG in the same directory)."
    )


def parse_version(version_str: str) -> tuple[int, int, int]:
    """Return (major, minor, patch) from a 'X.Y.Z' string."""
    parts = version_str.strip().split(".")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        raise SystemExit(
            f"error: version must be X.Y.Z with three numeric components, got {version_str!r}"
        )
    return int(parts[0]), int(parts[1]), int(parts[2])


def parse_zenodo_id(raw: str) -> str:
    """Return bare Zenodo record ID (digits only) from any supported input form."""
    # strip leading/trailing whitespace
    raw = raw.strip()
    # accept: 12345678 | zenodo.12345678 | 10.5281/zenodo.12345678
    m = re.search(r"zenodo\.(\d+)", raw)
    if m:
        return m.group(1)
    if re.fullmatch(r"\d+", raw):
        return raw
    raise SystemExit(
        f"error: cannot parse Zenodo ID from {raw!r}. "
        "Expected a bare record ID (digits), 'zenodo.<id>', or '10.5281/zenodo.<id>'."
    )


def read_file(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_file(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


# ---------------------------------------------------------------------------
# per-file update functions
# ---------------------------------------------------------------------------

def update_cmake(path: Path, major: int, minor: int, patch: int) -> None:
    text = read_file(path)

    replacements = [
        (r'(set\(HAMMER_VERSION_MAJOR\s+")\d+(")', rf'\g<1>{major}\2'),
        (r'(set\(HAMMER_VERSION_MINOR\s+")\d+(")', rf'\g<1>{minor}\2'),
        (r'(set\(HAMMER_VERSION_PATCH\s+")\d+(")', rf'\g<1>{patch}\2'),
    ]

    for pattern, repl in replacements:
        new_text, n = re.subn(pattern, repl, text)
        if n == 0:
            print(f"  warning: pattern not found in {path.name}: {pattern}")
        text = new_text

    write_file(path, text)
    print(f"  updated {path.name}: HAMMER_VERSION_MAJOR/MINOR/PATCH → {major}.{minor}.{patch}")


def update_readme_md(path: Path, major: int, minor: int, patch: int, zenodo_id: str) -> None:
    version_str = f"{major}.{minor}.{patch}"
    text = read_file(path)

    # 1. Version in the title heading
    new_text, n = re.subn(
        r'(HAMMER\s+-+\s+Helicity Amplitude Module for Matrix Element Reweighting,\s+ver\.)\s+[\d.]+',
        rf'\1 {version_str}',
        text,
    )
    if n == 0:
        print(f"  warning: title version pattern not found in {path.name}")
    text = new_text

    # 2. Version-specific DOI badge: the second [![DOI]...] line.
    #    We match any zenodo record ID on that line and replace it.
    #    The concept DOI (first badge) is left untouched.
    badge_pattern = re.compile(
        r'(\[!\[DOI\]\(https://zenodo\.org/badge/DOI/10\.5281/zenodo\.)(\d+)(\.svg\)\]\(https://doi\.org/10\.5281/zenodo\.)(\d+)(\))'
    )
    matches = list(badge_pattern.finditer(text))
    if len(matches) < 2:
        print(f"  warning: fewer than two DOI badges found in {path.name}; version-specific badge not updated")
    else:
        # Replace only the second badge (index 1)
        m = matches[1]
        replacement = f"{m.group(1)}{zenodo_id}{m.group(3)}{zenodo_id}{m.group(5)}"
        text = text[: m.start()] + replacement + text[m.end() :]
        print(f"  updated {path.name}: version-specific DOI → zenodo.{zenodo_id}")

    write_file(path, text)
    print(f"  updated {path.name}: title version → {version_str}")


def update_readme_plain(path: Path, major: int, minor: int, patch: int) -> None:
    version_str = f"{major}.{minor}.{patch}"
    text = read_file(path)

    new_text, n = re.subn(
        r'(HAMMER RELEASE VERSION)\s+[\d.]+',
        rf'\1 {version_str}',
        text,
    )
    if n == 0:
        print(f"  warning: 'HAMMER RELEASE VERSION' line not found in {path.name}")
    else:
        write_file(path, new_text)
        print(f"  updated {path.name}: HAMMER RELEASE VERSION → {version_str}")


def update_copyright_years(root: Path, new_year: int) -> None:
    """Replace the final year in all HAMMER copyright headers under *root*."""
    extensions = {".cc", ".hh", ".fhh", ".py", ".pyx", ".pxd", ".in"}
    plain_names = {"README", "COPYING"}

    # Match: "Copyright (C) YYYY - YYYY The HAMMER Collaboration" (any case)
    pattern = re.compile(
        r'(Copyright\s+\(C\)\s+\d{4}\s+-\s+)(\d{4})(\s+The\s+HAMMER\s+Collaboration)',
        re.IGNORECASE,
    )

    updated = 0
    skipped = 0

    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix not in extensions and path.name not in plain_names:
            continue
        # Skip vendored third-party trees
        if any(part in {"GTest", "flatbuffers", "pybind11"} for part in path.parts):
            continue

        text = path.read_text(encoding="utf-8", errors="replace")
        new_text, n = re.subn(
            pattern,
            lambda m: f"{m.group(1)}{new_year}{m.group(3)}",
            text,
        )
        if n > 0:
            path.write_text(new_text, encoding="utf-8")
            updated += 1
        else:
            skipped += 1

    print(f"  updated copyright year → {new_year} in {updated} file(s) ({skipped} file(s) had no match)")


def prepend_changelog(path: Path, major: int, minor: int, patch: int) -> None:
    version_str = f"{major}.{minor}.{patch}"
    today = date.today().strftime("%Y-%m-%d")

    stub = (
        f"{today} ---   Version {version_str}\n"
        f"     Additions:\n"
        f"      - \n"
        f"     Changes:\n"
        f"      - Bumped to version {version_str}, added version-specific doi\n"
        f"     Fixes:\n"
        f"      - \n"
        f"\n"
    )

    existing = read_file(path)
    write_file(path, stub + existing)
    print(f"  updated {path.name}: prepended stub entry for version {version_str}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bump the HAMMER version and DOI for a new release.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("version", help="New version string, e.g. 1.5.0")
    parser.add_argument("zenodo_id", help="New Zenodo record ID (bare, prefixed, or full DOI)")
    parser.add_argument(
        "--root",
        metavar="DIR",
        help="Path to the Hammer repository root (auto-detected if omitted)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be changed without modifying any files",
    )
    args = parser.parse_args()

    major, minor, patch = parse_version(args.version)
    zenodo_id = parse_zenodo_id(args.zenodo_id)

    root = Path(args.root).resolve() if args.root else find_repo_root(Path.cwd())
    print(f"Repository root: {root}")
    print(f"New version    : {major}.{minor}.{patch}")
    print(f"New Zenodo DOI : 10.5281/zenodo.{zenodo_id}")

    if args.dry_run:
        print("\n[dry-run] no files will be modified.")
        return

    current_year = date.today().year
    print()
    update_cmake(root / "CMakeLists.txt", major, minor, patch)
    update_readme_md(root / "README.md", major, minor, patch, zenodo_id)
    update_readme_plain(root / "README", major, minor, patch)
    prepend_changelog(root / "CHANGELOG", major, minor, patch)
    update_copyright_years(root, current_year)

    print(
        f"\nDone. Review the changes, fill in the CHANGELOG stub, "
        f"then commit and tag the release."
    )


if __name__ == "__main__":
    main()
